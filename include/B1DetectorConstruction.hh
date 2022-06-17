//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1DetectorConstruction.hh
/// \brief Definition of the B1DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include <vector>

class B1SD;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4PVPlacement;
class G4Tubs;
class G4Material;

using namespace std;
/// Detector construction class to define materials and geometry.

class B1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B1DetectorConstruction();
    virtual ~B1DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    G4LogicalVolume* GetScoringDetector() const { return fScoringDetector; }
    G4int return_num_detec() {return number_detectors; }

    G4double get_delta(){return delta;}
    G4double get_size() {return size;}
    G4double GetDetecPos(G4int num_detec);
    G4double get_position(){return posz;}
    std::vector<G4double>& GetPositions() {
      return positions;}
    G4double density_function(G4double h);
    G4double get_mass_overburden(G4double h);
/*
    struct detector_air {
      G4String name;
      G4double position;
      G4LogicalVolume* detector = new G4LogicalVolume(detector, world_mat, name);
      G4PVPlacement* place = new G4PVPlacement(0,
                        G4ThreeVector(0,0,(position)),
                        detector,
                          name,
                         fScoringVolume,
                         false,
                         1,
                         checkOverlaps);
    };
*/
  protected:
    G4LogicalVolume*  fScoringVolume;
    G4LogicalVolume* fScoringDetector;
    G4LogicalVolume* logicWorld;
    std::vector<G4double> positions;
    std::vector<G4LogicalVolume* > logicDetector;
    G4int number_detectors;
    G4double delta;
    std::vector<G4double> a_atm = {-186.5562*(g/cm2), -94.919*(g/cm2), 0.61289*(g/cm2), 0.0*(g/cm2), 0.01128292*(g/cm2)};
    std::vector<G4double> b_atm = {1222.6562*(g/cm2), 1144.9069*(g/cm2), 1305.5948*(g/cm2), 540.1778*(g/cm2), 1.0*(g/cm2)};
    std::vector<G4double> c_atm = {994186.38*cm, 878153.55*cm, 636143.04*cm, 772170.16*cm, 1.0e9*cm};
    std::vector<G4double> h_atm = {0*km, 4.*km, 10.*km, 40.*km, 100*km};
    G4double size;
    G4Tubs* detector;
    G4double posz;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
