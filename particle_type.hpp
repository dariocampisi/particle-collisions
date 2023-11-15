#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP

class ParticleType {
 public:
  ParticleType(const char *name, const double mass, const int charge)
      : fName{name}, fMass{mass}, fCharge{charge} {}

  // i getters che ritornano by value non devono ritornare const, perché se si
  // ritorna by value viene ritornata una copia dell'oggetto originale, che
  // quindi non può essere modificato in ogni caso
  // GetName() ritorna invece un puntatore, in questo caso quindi bisogna
  // ritornare const
  const char *GetName() const;
  double GetMass() const;
  int GetCharge() const;

  virtual double GetWidth() const { return 0; }

  // Print() sarà ridefinita nella in ResonanceType, derivata da ParticleType
  // se un metodo viene ridefinito in una classe derivata è necessario
  // dichiararlo virtual nella classe base
  virtual void Print() const;

 private:
  const char *fName;
  const double fMass;
  const int fCharge;
};

#endif
