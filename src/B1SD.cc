#include "B1SD.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

B1SD::B1SD(G4String SDname): G4VSensitiveDetector(SDname),
  hitCollection(nullptr), HCID(-1) {
  //cria a hit collection
  G4cout << "Criando Hit Collection com nome: " << SDname <<G4endl;
  collectionName.insert(SDname);
  track_id = 0;
  sdname = SDname;
}

B1SD::~B1SD() {

}

G4bool B1SD::ProcessHits(G4Step* step, G4TouchableHistory* ROhist) {
  G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
  //get name, momentum and energy
  const G4String particle_name = step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  G4ThreeVector momentum = step->GetTrack()->GetMomentumDirection();
  G4double energy = step->GetTrack()->GetTotalEnergy();
  //Get filename from RunAction class
  B1RunAction* runAction = (B1RunAction*) G4RunManager::GetRunManager()->GetUserRunAction();

  G4String filename = runAction->Get_filename();
  std::ofstream data("data.txt",std::ios_base::app);

  data << particle_name << " "  << " " << energy  << " " << momentum.getX() << " " << momentum.getY() << " " << momentum.getZ() << "\n";
//  G4cout << filename;
  //G4cout << particle_name << " "  << " " << energy  << " " << momentum.getX() << " " << momentum.getY() << " " << momentum.getZ() << " " << G4endl;
  data.close();

}

void B1SD::Initialize(G4HCofThisEvent* HCE) {
  hitCollection = new B1HitsCollection(SensitiveDetectorName,collectionName[0]);
  if (HCID < 0) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, hitCollection);
}

void B1SD::EndOfEvent(G4HCofThisEvent* HCE) {
  track_id = 0;
}
