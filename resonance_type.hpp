#ifndef RESONANCE_TYPE_HPP
#define RESONANCE_TYPE_HPP

#include "particle_type.hpp"

class ResonanceType : public ParticleType {
 public:
  ResonanceType(const char *name, const double mass, const int charge,
                const double width)
      : ParticleType(name, mass, charge), fWidth{width} {}

  double GetWidth() const;

  void Print() const;

 private:
  const double fWidth;
};

#endif
