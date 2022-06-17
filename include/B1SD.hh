#ifndef B1SD_h
#define B1SD_h 1

#include "G4VSensitiveDetector.hh"
#include "B1Hits.hh"
#include "B1RunAction.hh"
#include <vector>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class B1Hits;
class B1RunAction;
using namespace std;
class B1SD: public G4VSensitiveDetector
{
public:
  //construtor
  B1SD(G4String SDname);
  //destrutor
  ~B1SD();
  //função chamda a cada hit
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* ROhist);

  void Initialize(G4HCofThisEvent* HCE);

  //função chamada ao final do evento
  void EndOfEvent(G4HCofThisEvent* HCE);

private:
  B1HitsCollection* hitCollection;
  G4int HCID;
  std::vector<G4int> tracks;
  G4int track_id;
  G4String sdname;
};

#endif
