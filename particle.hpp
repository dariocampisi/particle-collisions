#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>

#include "particle_type.hpp"
#include "resonance_type.hpp"

class Particle {
 public:
  Particle() {}
  Particle(const char *name, double px = 0, double py = 0, double pz = 0);

  int GetIndex() const;
  int GetCharge() const;
  double GetPx() const;
  double GetPy() const;
  double GetPz() const;
  double GetP() const;
  double GetTrsP() const;
  double GetMass() const;

  void SetIndex(const int index);
  void SetIndex(const char *name);
  void SetP(double px, double py, double pz);

  static void AddParticleType(const char *name, const double mass,
                              const int charge, const double width = 0);

  static void PrintParticleTypes();
  void PrintParticleInfo() const;

  double TotalEnergy() const;
  double InvariantMass(const Particle &other) const;

  int Decay2body(Particle &dau1, Particle &dau2) const;

 private:
  static std::vector<ParticleType *> fParticleTypes;
  int fIndex;

  double fPx;
  double fPy;
  double fPz;

  int FindParticleType(const char *name) const;
  void Boost(double bx, double by, double bz);
};

#endif
