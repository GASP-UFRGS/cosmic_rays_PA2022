#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"

class B1Hits: public G4VHit {

public:
  B1Hits();
  ~B1Hits();

  void print();
  void set_partdef (const G4String particle_name);

  const G4String getParticleInTarget() {
    return fParticleInTarget;
  }
private:

 G4String fParticleInTarget;

};

typedef G4THitsCollection<B1Hits> B1HitsCollection;
