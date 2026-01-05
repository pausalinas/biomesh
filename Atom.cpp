#include "biomesh/Atom.hpp"

namespace biomesh {

Atom::Atom(const std::string& element) 
    : chemicalElement_(element), x_(0.0), y_(0.0), z_(0.0), atomicRadius_(0.0), atomicMass_(0.0), id_(0) {
}

Atom::Atom(const std::string& element, double radius) 
    : chemicalElement_(element), x_(0.0), y_(0.0), z_(0.0), atomicRadius_(radius), atomicMass_(0.0), id_(0) {
}

Atom::Atom(const std::string& element, double radius, double mass) 
    : chemicalElement_(element), x_(0.0), y_(0.0), z_(0.0), atomicRadius_(radius), atomicMass_(mass), id_(0) {
}

void Atom::setCoordinates(double x, double y, double z) {
    x_ = x;
    y_ = y;
    z_ = z;
}

} // namespace biomesh