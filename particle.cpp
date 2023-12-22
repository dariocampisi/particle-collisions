#include "particle.hpp"

#include <cmath>    // for M_PI
#include <cstdlib>  // for RAND_MAX
#include <iostream>

std::vector<ParticleType *> Particle::fParticleTypes;

Particle::Particle(const char *name, double px, double py, double pz)
    : fPx{px}, fPy{py}, fPz{pz} {
  if (FindParticleType(name) == -1) {
    throw std::invalid_argument("particle type not found\n");
  } else {
    fIndex = FindParticleType(name);
  }
}

int Particle::FindParticleType(const char *name) const {
  for (auto it = fParticleTypes.begin(); it != fParticleTypes.end(); ++it) {
    if ((*it)->GetName() == name) {
      return (int)(it - fParticleTypes.begin());
    }
  }
  return -1;
}

int Particle::GetIndex() const { return fIndex; }
int Particle::GetCharge() const { return fParticleTypes[fIndex]->GetCharge(); }
double Particle::GetPx() const { return fPx; }
double Particle::GetPy() const { return fPy; }
double Particle::GetPz() const { return fPz; }
double Particle::GetP() const {
  return std::sqrt(fPx * fPx + fPy * fPy + fPz * fPz);
}
double Particle::GetTrsP() const { return std::sqrt(fPx * fPx + fPy * fPy); }
double Particle::GetMass() const { return fParticleTypes[fIndex]->GetMass(); }

void Particle::SetIndex(const int index) {
  if (index >= 0 && index < (int)(fParticleTypes.size())) {
    fIndex = index;
  } else {
    std::cerr << "Index not found\n";
  }
}
void Particle::SetIndex(const char *name) {
  fIndex = FindParticleType(name);

  if (fIndex == -1) {
    std::cerr << "Particle name not found\n";
  }
}

void Particle::SetP(double px, double py, double pz) {
  fPx = px;
  fPy = py;
  fPz = pz;
}

void Particle::AddParticleType(const char *name, const double mass,
                               const int charge, const double width) {
  ParticleType *type = new ResonanceType(name, mass, charge, width);
  fParticleTypes.push_back(type);
}

void Particle::PrintParticleTypes() {
  for (auto type : fParticleTypes) {
    type->Print();
  }
}

void Particle::PrintParticleInfo() const {
  std::cout << fIndex << " " << fParticleTypes[fIndex]->GetName()
            << "Impulse (Px,Py,Pz): (" << fPx << ", " << fPy << ", " << fPz
            << ")\n";
}

double Particle::TotalEnergy() const {
  return std::sqrt((this->GetMass() * this->GetMass()) +
                   (fPx * fPx + fPy * fPy + fPz * fPz));
}

double Particle::InvariantMass(const Particle &other) const {
  return std::sqrt((this->TotalEnergy() + other.TotalEnergy()) *
                       (this->TotalEnergy() + other.TotalEnergy()) -
                   ((fPx + other.GetPx()) * (fPx + other.GetPx()) +
                    (fPy + other.GetPy()) * (fPy + other.GetPy()) +
                    (fPz + other.GetPz()) * (fPz + other.GetPz())));
}

int Particle::Decay2body(Particle &dau1, Particle &dau2) const {
  if (this->GetMass() == 0.0) {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }

  double massMot = this->GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if (fIndex > -1) {  // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * std::rand() * invnum - 1.0;
      x2 = 2.0 * std::rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = std::sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;

    massMot += fParticleTypes[fIndex]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    printf(
        "Decayment cannot be preformed because mass is too low in this "
        "channel\n");
    return 2;
  }

  double pout =
      std::sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = std::rand() * norm;
  double theta = std::rand() * norm * 0.5 - M_PI / 2.;
  dau1.SetP(pout * std::sin(theta) * std::cos(phi),
            pout * std::sin(theta) * std::sin(phi), pout * std::cos(theta));
  dau2.SetP(-pout * std::sin(theta) * std::cos(phi),
            -pout * std::sin(theta) * std::sin(phi), -pout * std::cos(theta));

  double energy =
      std::sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

  double bx = fPx / energy;
  double by = fPy / energy;
  double bz = fPz / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}

void Particle::Boost(double bx, double by, double bz) {
  double energy = this->TotalEnergy();

  // Boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / std::sqrt(1.0 - b2);
  double bp = bx * fPx + by * fPy + bz * fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  fPx += gamma2 * bp * bx + gamma * bx * energy;
  fPy += gamma2 * bp * by + gamma * by * energy;
  fPz += gamma2 * bp * bz + gamma * bz * energy;
}
