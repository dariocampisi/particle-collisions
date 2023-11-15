#include "particle_type.hpp"

#include <iostream>

const char *ParticleType::GetName() const { return fName; }

double ParticleType::GetMass() const { return fMass; }

int ParticleType::GetCharge() const { return fCharge; }

void ParticleType::Print() const {
  std::cout << "Name: " << fName << "\nMass: " << fMass
            << "\nCharge: " << fCharge << '\n';
}
