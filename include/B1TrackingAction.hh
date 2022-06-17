

#ifndef TrackingAction_h
#define TrackingAction_h

#include "G4UserTrackingAction.hh"
#include "G4LogicalVolume.hh"
#include <map>

class G4Region;
class G4ParticleDefinition;
class B1DetectorConstruction;


class B1TrackingAction : public G4UserTrackingAction
{
private:
  //aqui ficam as  variáveis da classe
  B1DetectorConstruction* fworld;
  G4Region* fworldRegion;
  G4LogicalVolume* worldVolume;
  const G4ParticleDefinition* particle_def;
  std::map<const G4ParticleDefinition*, int> fNParticleOutsideWorld;
  std::map<const G4ParticleDefinition*, int> fNParticleInWorld;

public:
 //aqui ficam as funçes que a classe terá
    B1TrackingAction(B1DetectorConstruction* world); //construtor
    ~B1TrackingAction(); //destrutor
    virtual void PreUserTrackingAction(const G4Track* track);

    std::map<const G4ParticleDefinition*, int>& GetNParticlesCreatedOutsideWorld()
    {
        return fNParticleOutsideWorld;
    }

    std::map<const G4ParticleDefinition*, int>& GetNParticlesCreatedInWorld()
    {
        return fNParticleInWorld;
    }
    void clearParticles();
};

#endif
