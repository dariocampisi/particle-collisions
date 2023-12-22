#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP

class ParticleType {
 public:
  ParticleType(const char *name, const double mass, const int charge)
      : fName{name}, fMass{mass}, fCharge{charge} {}

  const char *GetName() const;
  double GetMass() const;
  int GetCharge() const;

  virtual double GetWidth() const { return 0; }

  virtual void Print() const;

 private:
  const char *fName;
  const double fMass;
  const int fCharge;
};

#endif
