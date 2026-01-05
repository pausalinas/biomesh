#!/usr/bin/env bash
set -o pipefail

SCRIPT_NAME=$(basename "$0")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CONFIG="${SCRIPT_DIR}/config/default.conf"
VERSION="1.0.0"

print_banner() {
    cat <<'EOF'
========================
   BioMesh CLI Runner
========================
EOF
}

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [options] <input.pdb>

Options:
  -i, --input FILE           Input PDB file (positional also accepted)
  -o, --output BASENAME      Output file basename (default: mesh)
                             Generates: <basename>_occupied.msh and <basename>_empty.msh
  -v, --voxel-size VALUE     Voxel size in Angstroms (default: 1.0)
  -p, --padding VALUE        Bounding box padding in Angstroms (default: 2.0)
  -f, --format FORMAT        Output format (vtk|gid|json|txt) (default: gid)
      --filter TYPE          Filter type (none|all|protein-only|no-water|custom)
      --occupied BOOL        Generate occupied mesh (true|false) (default: true)
      --empty BOOL           Generate empty mesh (true|false) (default: true)
  -c, --config FILE          INI-style config file
      --generate-config      Print default config and exit
      --batch LIST           Comma-separated list of PDB files to process
  -V, --verbose              Verbose output
  -q, --quiet                Minimal output
      --version              Show version
  -h, --help                 Show this help message

Examples:
  # Generate both occupied and empty meshes
  ${SCRIPT_NAME} protein.pdb

  # Custom output basename
  ${SCRIPT_NAME} protein.pdb -o my_protein
  # → my_protein_occupied.msh and my_protein_empty.msh

  # Generate only occupied mesh
  ${SCRIPT_NAME} protein.pdb --empty false

  # Generate only empty mesh
  ${SCRIPT_NAME} protein.pdb --occupied false

  # Use config file
  ${SCRIPT_NAME} --config config/default.conf

  # Batch processing
  ${SCRIPT_NAME} --batch file1.pdb,file2.pdb -p 1.5
EOF
}

die() { echo "Error: $1" >&2; exit 1; }

trim() { echo "$1" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'; }

declare -A CONFIG_MAP
parse_config_file() {
    local file="$1"
    [[ -f "$file" ]] || die "Config file not found: $file"
    local section=""
    while IFS= read -r line || [[ -n "$line" ]]; do
        line="${line%%#*}"
        line="${line%%;*}"
        line=$(trim "$line")
        [[ -z "$line" ]] && continue
        if [[ "$line" =~ ^\[(.+)\]$ ]]; then
            section="${BASH_REMATCH[1]}"
            continue
        fi
        if [[ "$line" =~ ^([^=]+)=(.*)$ ]]; then
            local key
            key=$(trim "${BASH_REMATCH[1]}")
            local value
            value=$(trim "${BASH_REMATCH[2]}")
            CONFIG_MAP["${section}.${key}"]="$value"
        fi
    done <"$file"
}

set_from_config() {
    local key="$1" var_name="$2"
    if [[ -n "${CONFIG_MAP[$key]:-}" ]]; then
        printf -v "$var_name" "%s" "${CONFIG_MAP[$key]}"
    fi
}

log() { $quiet || echo "$@"; }
logv() { $verbose && ! $quiet && echo "$@"; }

generate_default_config() {
    cat <<'EOF'
[input]
# Input PDB file
file = input.pdb

[filter]
type = protein-only
keep_proteins = true
keep_water = false
keep_nucleic_acids = false
keep_ions = false
keep_lipids = true
keep_ligands = true

[voxelization]
voxel_size = 1.0
padding = 2.0

[output]
# Base name for output files (without extension)
# Will generate: <basename>_occupied.msh and <basename>_empty.msh
basename = mesh
format = gid

# Control which meshes to generate
generate_occupied = true
generate_empty = true

[options]
verbose = true
EOF
}

is_number() { [[ "$1" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; }

allowed_filter() {
    case "$1" in
    none|all|protein-only|no-water|custom) return 0 ;;
    *) return 1 ;;
    esac
}

allowed_format() {
    case "$1" in
    vtk|gid|json|txt) return 0 ;;
    *) return 1 ;;
    esac
}

normalize_bool() {
    case "$(echo "$1" | tr '[:upper:]' '[:lower:]')" in
    1|true|yes|on) echo "true" ;;
    0|false|no|off|"") echo "false" ;;
    *) echo "$1" ;;
    esac
}

append_batch_items() {
    local value="$1"
    IFS=',' read -r -a parts <<<"$value"
    for part in "${parts[@]}"; do
        part=$(trim "$part")
        [[ -z "$part" ]] && continue
        cli_batch+=("$part")
    done
}

find_mesh_executables() {
    local build_dir="${SCRIPT_DIR}/build"
    OCCUPIED_CMD=""
    EMPTY_CMD=""
    
    # Find occupied voxel mesh generator
    if [[ -x "${build_dir}/occupied_voxel_to_gid" ]]; then
        OCCUPIED_CMD="${build_dir}/occupied_voxel_to_gid"
        logv "Found occupied mesh generator: $OCCUPIED_CMD"
    else
        logv "Warning: occupied_voxel_to_gid not found or not executable"
    fi
    
    # Find empty voxel mesh generator
    if [[ -x "${build_dir}/empty_voxel_to_gid" ]]; then
        EMPTY_CMD="${build_dir}/empty_voxel_to_gid"
        logv "Found empty mesh generator: $EMPTY_CMD"
    else
        logv "Warning: empty_voxel_to_gid not found or not executable"
    fi
    
    # Fail if both are missing
    if [[ -z "$OCCUPIED_CMD" && -z "$EMPTY_CMD" ]]; then
        die "No mesh generator executables found. Build the project in ${build_dir} first."
    fi
}

write_summary() {
    local mesh_type="$1"
    local input_file="$2"
    local output_file="$3"
    local summary_file="${output_file}.summary"
    
    mkdir -p "$(dirname "$summary_file")" 2>/dev/null || true
    cat >"$summary_file" <<EOF
# BioMesh dual mesh generation summary
Mesh type: $mesh_type
Input file: $input_file
Output file: $output_file
Output format: $output_format
Voxel size: $voxel_size
Padding: $padding
Filter: $filter_type
Generation timestamp: $(date -u +"%Y-%m-%d %H:%M:%S UTC")
EOF
    logv "Wrote summary to $summary_file"
}

get_file_extension() {
    case "$output_format" in
    gid) echo "msh" ;;
    vtk) echo "vtk" ;;
    json) echo "json" ;;
    txt) echo "txt" ;;
    *) echo "msh" ;;  # Default to msh
    esac
}

run_dual_mesh_generation() {
    local input_file="$1"
    local run_index="$2"
    local run_total="$3"
    
    # Determine base filename for output and file extension
    local base="$output_basename"
    local ext
    ext=$(get_file_extension)
    local stem=""
    if (( run_total > 1 )); then
        stem=$(basename "${input_file%.*}")
    fi
    
    # Track generated files for summary
    local generated_files=()
    
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log "Processing: $input_file"
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    # Generate occupied mesh if enabled
    if [[ "$generate_occupied" == "true" ]]; then
        if [[ -n "$OCCUPIED_CMD" ]]; then
            local occupied_output
            if (( run_total > 1 )); then
                occupied_output="${base}_${stem}_occupied.${ext}"
            else
                occupied_output="${base}_occupied.${ext}"
            fi
            
            log ""
            log "═══ Generating Occupied Mesh ═══"
            logv "Running: $OCCUPIED_CMD $input_file $voxel_size $occupied_output $padding"
            
            if "$OCCUPIED_CMD" "$input_file" "$voxel_size" "$occupied_output" "$padding"; then
                log "✓ Occupied mesh generated: $occupied_output"
                write_summary "occupied" "$input_file" "$occupied_output"
                generated_files+=("$occupied_output")
            else
                die "Occupied mesh generation failed for: $input_file"
            fi
        else
            log "⚠ Warning: occupied_voxel_to_gid not available, skipping occupied mesh"
        fi
    else
        logv "Skipping occupied mesh generation (disabled)"
    fi
    
    # Generate empty mesh if enabled
    if [[ "$generate_empty" == "true" ]]; then
        if [[ -n "$EMPTY_CMD" ]]; then
            local empty_output
            if (( run_total > 1 )); then
                empty_output="${base}_${stem}_empty.${ext}"
            else
                empty_output="${base}_empty.${ext}"
            fi
            
            log ""
            log "═══ Generating Empty Mesh ═══"
            logv "Running: $EMPTY_CMD $input_file $voxel_size $empty_output $padding"
            
            if "$EMPTY_CMD" "$input_file" "$voxel_size" "$empty_output" "$padding"; then
                log "✓ Empty mesh generated: $empty_output"
                write_summary "empty" "$input_file" "$empty_output"
                generated_files+=("$empty_output")
            else
                die "Empty mesh generation failed for: $input_file"
            fi
        else
            log "⚠ Warning: empty_voxel_to_gid not available, skipping empty mesh"
        fi
    else
        logv "Skipping empty mesh generation (disabled)"
    fi
    
    log ""
    log "✓ Completed processing: $input_file"
    if (( ${#generated_files[@]} > 0 )); then
        log "  Generated files:"
        for f in "${generated_files[@]}"; do
            log "    - $f"
        done
    fi
    log ""
}

# Defaults
input_file=""
output_basename=""  # Will be set to "mesh" if still empty after config parsing
voxel_size="1.0"
padding="2.0"
filter_type="none"
output_format="gid"
generate_occupied="true"
generate_empty="true"
config_file=""
verbose=false
quiet=false
batch_inputs=()
cli_input=""
cli_output=""
cli_voxel=""
cli_padding=""
cli_filter=""
cli_format=""
cli_occupied=""
cli_empty=""
cli_batch=()
cli_verbose_set=false
cli_quiet_set=false

[[ -f "$DEFAULT_CONFIG" ]] && config_file="$DEFAULT_CONFIG"

require_value() {
    local option="$1"
    local next_index="$2"
    local total="$3"
    (( next_index < total )) || die "Option ${option} requires a value."
}

args=("$@")
idx=0
while [[ $idx -lt ${#args[@]} ]]; do
    arg="${args[$idx]}"
    case "$arg" in
    -i|--input)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_input="${args[$idx]}"
        ;;
    --input=*)
        cli_input="${arg#*=}"
        ;;
    -o|--output)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_output="${args[$idx]}"
        ;;
    --output=*)
        cli_output="${arg#*=}"
        ;;
    -v|--voxel-size)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_voxel="${args[$idx]}"
        ;;
    --voxel-size=*)
        cli_voxel="${arg#*=}"
        ;;
    -p|--padding)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_padding="${args[$idx]}"
        ;;
    --padding=*)
        cli_padding="${arg#*=}"
        ;;
    -f|--format)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_format="${args[$idx]}"
        ;;
    --format=*)
        cli_format="${arg#*=}"
        ;;
    --filter)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_filter="${args[$idx]}"
        ;;
    --filter=*)
        cli_filter="${arg#*=}"
        ;;
    --occupied)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_occupied="${args[$idx]}"
        ;;
    --occupied=*)
        cli_occupied="${arg#*=}"
        ;;
    --empty)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); cli_empty="${args[$idx]}"
        ;;
    --empty=*)
        cli_empty="${arg#*=}"
        ;;
    -c|--config)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); config_file="${args[$idx]}"
        ;;
    --config=*)
        config_file="${arg#*=}"
        ;;
    --generate-config)
        generate_default_config
        exit 0
        ;;
    --batch)
        require_value "$arg" $((idx + 1)) ${#args[@]}
        ((idx++)); append_batch_items "${args[$idx]}"
        ;;
    --batch=*)
        append_batch_items "${arg#*=}"
        ;;
    -V|--verbose)
        cli_verbose_set=true
        verbose=true
        ;;
    -q|--quiet)
        cli_quiet_set=true
        quiet=true
        ;;
    --version)
        echo "${SCRIPT_NAME} version ${VERSION}"
        exit 0
        ;;
    -h|--help)
        usage
        exit 0
        ;;
    --)
        ((idx++))
        for ((; idx < ${#args[@]}; idx++)); do
            cli_batch+=("${args[$idx]}")
        done
        break
        ;;
    -*)
        die "Unknown option: $arg"
        ;;
    *)
        if [[ -z "$cli_input" ]]; then
            cli_input="$arg"
        else
            append_batch_items "$arg"
        fi
        ;;
    esac
    ((idx++))
done

if [[ -n "$config_file" ]]; then
    parse_config_file "$config_file"
    set_from_config "input.file" input_file
    set_from_config "filter.type" filter_type
    set_from_config "voxelization.voxel_size" voxel_size
    set_from_config "voxelization.padding" padding
    
    # Support both old (output.file) and new (output.basename) config format
    set_from_config "output.basename" output_basename
    if [[ -z "$output_basename" ]]; then
        # Fallback for backward compatibility with old output.file format
        old_output_file=""
        set_from_config "output.file" old_output_file
        if [[ -n "$old_output_file" ]]; then
            # Extract basename without extension, handling paths with multiple dots
            # e.g., /path/to/file.v1.0.msh -> /path/to/file.v1.0
            base_with_path="${old_output_file%.*}"
            # If the original had a directory component, preserve it
            if [[ "$old_output_file" == */* ]]; then
                output_basename="$base_with_path"
            else
                # For simple filenames, just remove the extension
                output_basename="$base_with_path"
            fi
        fi
    fi
    
    set_from_config "output.format" output_format
    set_from_config "output.generate_occupied" generate_occupied
    set_from_config "output.generate_empty" generate_empty
    set_from_config "options.verbose" verbose
    set_from_config "options.quiet" quiet
    
    # Normalize boolean values
    verbose=$(normalize_bool "$verbose")
    quiet=$(normalize_bool "$quiet")
    generate_occupied=$(normalize_bool "$generate_occupied")
    generate_empty=$(normalize_bool "$generate_empty")
fi

[[ -n "$cli_input" ]] && input_file="$cli_input"
[[ ${#cli_batch[@]} -gt 0 ]] && batch_inputs=("${cli_batch[@]}")
[[ -n "$cli_output" ]] && output_basename="$cli_output"
[[ -n "$cli_voxel" ]] && voxel_size="$cli_voxel"
[[ -n "$cli_padding" ]] && padding="$cli_padding"
[[ -n "$cli_filter" ]] && filter_type="$cli_filter"
[[ -n "$cli_format" ]] && output_format="$cli_format"
[[ -n "$cli_occupied" ]] && generate_occupied=$(normalize_bool "$cli_occupied")
[[ -n "$cli_empty" ]] && generate_empty=$(normalize_bool "$cli_empty")
$cli_verbose_set && verbose=true
$cli_quiet_set && quiet=true

# Set default output basename if not set by config or CLI
[[ -z "$output_basename" ]] && output_basename="mesh"

all_inputs=()
if [[ -n "$input_file" ]]; then
    all_inputs+=("$input_file")
fi
if [[ ${#batch_inputs[@]} -gt 0 ]]; then
    all_inputs+=("${batch_inputs[@]}")
fi

if [[ ${#all_inputs[@]} -eq 0 ]]; then
    usage
    die "No input PDB file specified."
fi

for f in "${all_inputs[@]}"; do
    [[ -f "$f" ]] || die "Input file not found: $f"
done

is_number "$voxel_size" || die "Voxel size must be numeric."
is_number "$padding" || die "Padding must be numeric."
allowed_filter "$filter_type" || die "Invalid filter type: $filter_type"
allowed_format "$output_format" || die "Invalid output format: $output_format"

# Check that at least one mesh type is enabled
if [[ "$generate_occupied" != "true" && "$generate_empty" != "true" ]]; then
    die "At least one mesh type must be enabled (--occupied true or --empty true)"
fi

print_banner

log "BioMesh Dual Mesh Generation:"
log "  Inputs             : ${all_inputs[*]}"
log "  Output basename    : $output_basename"
log "  Voxel size         : $voxel_size Angstroms"
log "  Padding            : $padding Angstroms"
log "  Filter type        : $filter_type"
log "  Format             : $output_format"
log "  Generate occupied  : $generate_occupied"
log "  Generate empty     : $generate_empty"
log "  Verbose            : $verbose"
log "  Quiet              : $quiet"
log "  Config file        : ${config_file:-<none>}"
log ""

find_mesh_executables

log "Available executables:"
[[ -n "$OCCUPIED_CMD" ]] && log "  ✓ Occupied mesh: $OCCUPIED_CMD" || log "  ✗ Occupied mesh: not found"
[[ -n "$EMPTY_CMD" ]] && log "  ✓ Empty mesh: $EMPTY_CMD" || log "  ✗ Empty mesh: not found"
log ""

for idx in "${!all_inputs[@]}"; do
    run_dual_mesh_generation "${all_inputs[$idx]}" "$((idx + 1))" "${#all_inputs[@]}"
done

log "═══════════════════════════════════════════════════════════"
log "All mesh generation tasks completed successfully."
log "═══════════════════════════════════════════════════════════"
