
#include "B1TrackingAction.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4Region.hh"
#include "B1DetectorConstruction.hh"

using namespace std;

//implementação do construtor
B1TrackingAction::B1TrackingAction(B1DetectorConstruction* world) {
  fworld = world;
  fworldRegion = 0;
   worldVolume = 0;
}

//implementação do destrutor
B1TrackingAction::~B1TrackingAction() {
  fworld = 0;
  fworldRegion = 0;
   worldVolume = 0;
}

//implementação da função(instância da classe B1TrackingAction) PreUserTrackingAction
void B1TrackingAction::PreUserTrackingAction(const G4Track* track) {
  //pega o nome da partícula relacionado ao seu TrackID.
  /*
  particle_def = track->GetParticleDefinition();
  G4int particle_track = track->GetTrackID();
  //const G4String& particle_name =	particle_def->GetParticleName();
  //se a variável fworldRegio não possui a regiao do mundo, fazer isso
  if(fworldRegion == 0) {
    worldVolume= fworld->GetScoringVolume();
    fworldRegion = worldVolume->GetRegion();
  }


  const G4ThreeVector& position = track->GetPosition();
  int N =  fworldRegion->GetNumberOfRootVolumes();
  std::vector<G4LogicalVolume*>::iterator it_logicalVolumeInRegion = fworldRegion->GetRootLogicalVolumeIterator();
  bool inside_world = false;
  //checar se está dentro do mundo
  for(int i = 0; i < N ; i++, it_logicalVolumeInRegion++)
        {
            EInside test_status = (*it_logicalVolumeInRegion)->GetSolid()->Inside(position) ;
            if(test_status == kInside)
            {
                inside_world = true;

                break;
            }
        }
        //coloca eles no std::map
        if (inside_world) {
          fNParticleInWorld[particle_def]++;
        } else {
          fNParticleOutsideWorld[particle_def]++;
        }
 */
}
void B1TrackingAction::clearParticles() {
  /*
  fNParticleInWorld.clear();
  fNParticleOutsideWorld.clear();
  */

}
