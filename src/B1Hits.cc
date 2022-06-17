#include "B1Hits.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


B1Hits::B1Hits() {

}
B1Hits::~B1Hits() {

}

 void B1Hits::set_partdef(const G4String particle_name) {
      fParticleInTarget = particle_name;
}
