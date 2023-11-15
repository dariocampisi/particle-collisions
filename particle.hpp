#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>

#include "particle_type.hpp"
#include "resonance_type.hpp"

class Particle {
 public:
  Particle() {}  // default constructor, meglio definirlo esplicitamente
  Particle(const char *name, double px = 0, double py = 0, double pz = 0);

  int GetIndex() const;
  int GetPx() const;
  int GetPy() const;
  int GetPz() const;
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

  // dalla spiegazione dell'Arcelli
  // Il metodo simiula la cinematica del decadimento di una particella p1 in
  // altre due, p2 e p3 quindi, p1.(p2,p3) fa decadere p1 in p2 e p3 dopo la
  // chiamata p2 e p3 contengono gli impulsi finali delle figlie del decadimento
  // se tutto è andato a buon fine la funzione ritorna 0, altrimenti ritorna un
  // numero > 0
  int Decay2body(Particle &dau1, Particle &dau2) const;

 private:
  static std::vector<ParticleType *> fParticleTypes;
  int fIndex;

  double fPx;
  double fPy;
  double fPz;

  // deve essere privato, è scritto nella consegna
  int FindParticleType(const char *name) const;

  // serve a Decay2body()
  void Boost(double bx, double by, double bz);
};

#endif
