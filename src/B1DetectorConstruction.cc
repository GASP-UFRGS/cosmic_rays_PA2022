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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "B1SD.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Track.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Element.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //making air materiais
  //Nitrogen
  G4double a_N = 14.01*g/mole;
  G4String name, symbol;
  G4double z;
  G4Element* ele_N = new G4Element(name="Nitrogen", symbol="N", z=7., a_N);
  //Oxygen
  G4double a_o = 16.00*g/mole;
  G4Element* ele_O = new G4Element(name="Oxygen", symbol="O", z=8., a_o);
  //Argon
  G4double a_Ar = 39.948*g/mole;
  G4Element* ele_Ar = new G4Element(name="Argon", symbol="Ar", z=18., a_Ar);

  G4Material* air_material = nist->FindOrBuildMaterial("G4_AIR");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  //caracteristicas do cilíndro -> world
  G4double raio_i = 0;
  G4double raio_e = 20*m; // RAIO EXTERNO CILINDRO
  //G4double h_detector = 40*m;
  G4double h=40*m; // ALTURA CILINDRO
  G4double theta_0 = 0.*deg;
  G4double theta_f = 360.*deg;
  //making world of vacuum
  posz = h/2;

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* Oxygen = nist->FindOrBuildMaterial("G4_O");
  G4Material* Nitrogen = nist->FindOrBuildMaterial("G4_N");
  G4Material* Argon = nist->FindOrBuildMaterial("G4_Ar");

  G4Tubs* cilindro = new G4Tubs("world", raio_i, raio_e, h/2, theta_0, theta_f);

    logicWorld =
    new G4LogicalVolume(cilindro,          //its solid
                        air_material,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  //MUNDO INVISIVEL
  /*
//Camadas

  G4int num_layers = 250;
  std::vector<G4LogicalVolume* > air_layers;
  std::vector<G4Material *> layers_material;
  layers_material.reserve(num_layers);
  air_layers.reserve(num_layers);
  G4double dh = h/num_layers;
  G4Tubs* cilindro_ar = new G4Tubs("air_layer", raio_i, raio_e, dh/2 , theta_0, theta_f);
  std::ofstream density_file("density.csv");
  std::ofstream mass_file("mass.csv");

  for(G4int i=0; i<num_layers; i++){
     //Calculates density
     G4double atm_position = (dh/2) + dh*i;
     G4double density  = density_function(atm_position);
     G4double mass_overburden = get_mass_overburden(atm_position);
     density_file << density/(g/cm3) << " " << atm_position/(km) << "\n";
     mass_file << mass_overburden/(g/cm2) << " " << atm_position/(km) << "\n";
     //make layer material
     layers_material[i] = new G4Material(name="Air" + std::to_string(i), density, 3);
     layers_material[i]->AddMaterial(Nitrogen, 78.1*perCent);
     layers_material[i]->AddMaterial(Oxygen, 21.0*perCent);
     layers_material[i]->AddMaterial(Argon, 0.9*perCent);

      air_layers[i] = new G4LogicalVolume(cilindro_ar,          //its solid
                          layers_material[i],           //its material
                        "layer" + std::to_string(i + 1));

       //calculing the position in relation to the mother volume

       G4double position = (h/2) - (dh/2) -dh*i;



        new G4PVPlacement(0,                     //no rotation
                         G4ThreeVector(0,0,position),       //at (0,0,0)
                         air_layers[i],            //its logical volume
                         ("lay" + std::to_string(i + 1)),               //its name
                         logicWorld,                     //its mother  volume
                          false,                 //no boolean operation
                          i,
                          checkOverlaps);
                                 //copynumber




  }
 density_file.close();
*/
  //
  //always return the physical World
  //
   fScoringVolume = logicWorld;

   //começar aqui os layers da atmosfera em um for semelhante ao dos detetores, porém grudados e preenchendo todo o cilindro.




/*
  G4double safe_distance = 0.02*m;
  //Coloque aqui o número de detetores
  G4int num_detector = 0;
  number_detectors = num_detector;

  G4double total_safe = h_detector + safe_distance;
  G4double delta_h = (2*h - 2*safe_distance - 4*h_detector)/(num_detector - 1);
  G4double initial_pos = -h + total_safe;
  G4double final_pos = h - total_safe;
  size = h;
  delta = delta_h;




  /*
  G4LogicalVolume* logicDetector = new G4LogicalVolume(detector,          //its solid
                     world_mat,           //its material
                   "logicDetector");


//Criando os Logical volume
/*
logicDetector.reserve(num_detector);
 for(G4int i=0; i < num_detector; i++) {
    logicDetector[i] = new G4LogicalVolume(detector,          //its solid
                       world_mat,           //its material
                     "logicDetector" + std::to_string(i + 1));

    positions.push_back(initial_pos + i*delta_h);
    std::cout << positions[i] << '\n';

    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,positions[i]),       //at (0,0,0)
                      logicDetector[i],            //its logical volume
                      ("detector" + std::to_string(i + 1)),               //its name
                      logicWorld,                     //its mother  volume
                      false,                 //no boolean operation
                      i,                     //copy number
                      checkOverlaps);        //overlaps checking

  }



//Coloca os deteroes dentro do cilindro e coloca a sua posição em um vetor com push_back()

  //fScoringDetector = logicDetector;
 //Associando os detetores à classe sensitive detector
  auto sdman = G4SDManager::GetSDMpointer();


  for(G4int i = 0; i< num_detector; i++) {
  G4String SDname = "SD" + std::to_string(i + 1);
  auto sensitive = new B1SD(SDname);
  sdman->AddNewDetector(sensitive);
  logicDetector[i]->SetSensitiveDetector(sensitive);

}
    //colocando cor vermelha para detectores

  G4VisAttributes* worldVisAtt1 = new G4VisAttributes(G4Colour(1.0,0.0,0.0));

  // DEIXA DETETORES INVISIVEIS
  worldVisAtt1->SetVisibility(true);
    // DEIXA DETETORES INVISIVEIS
  //worldVisAtt1->SetVisibility(false);

  for(G4int i = 0; i < num_detector; i++) {
  logicDetector[i]->SetVisAttributes(worldVisAtt1 );
}

// Detector on the ground
G4int number_detectors = 1;

G4Material* detec_mat = nist->FindOrBuildMaterial("G4_Galactic");
G4double detec_height = 1*cm;
G4double detec_z = (dh/2) - (detec_height/2);

  G4Tubs* detec_tub = new G4Tubs("detector", raio_i, raio_e, detec_height/2, theta_0, theta_f);

  G4LogicalVolume* logicDetector =
  new G4LogicalVolume(detec_tub,          //its solid
                      detec_mat,           //its material
                      "detec_volume");            //its name

// air_layers[0] -> First layer(On the ground)
G4VPhysicalVolume* physDetector =
  new G4PVPlacement(0,                     //no rotation
                    G4ThreeVector(0,0,detec_z),       //at (0,0,0)
                    logicDetector,            //its logical volume
                    "detector",               //its name
                    air_layers[0],                     //its mother  volume
                    false,                 //no boolean operation
                    0,                     //copy number
                    checkOverlaps);        //overlaps checking


//Applying red color
G4VisAttributes* worldVisAtt1 = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
worldVisAtt1->SetVisibility(true);
logicDetector->SetVisAttributes(worldVisAtt1);

//Set as sensitive detector

auto sdman = G4SDManager::GetSDMpointer();
G4String SDname = "SD" + std::to_string(number_detectors);
auto sensitive = new B1SD(SDname);
sdman->AddNewDetector(sensitive);
logicDetector->SetSensitiveDetector(sensitive);
*/

  return physWorld;
}

void B1DetectorConstruction::ConstructSDandField() {
  //creanting uniform magnetic field
/*
  G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(.0000005*LT,0.,0.));
  G4FieldManager* FieldMgr  =new G4FieldManager(magField);
  logicWorld->SetFieldManager(FieldMgr, true);
  */
}

G4double B1DetectorConstruction::density_function(G4double h) {
  G4int l; //layer
  G4double res;
  G4cout << h << G4endl;
  if (h >= h_atm[0] && h < h_atm[1]) {
    l = 0;
  } else if (h >= h_atm[1] && h < h_atm[2]) {
    l = 1;
  } else if (h >= h_atm[2] && h < h_atm[3]) {
    l = 2;
  } else if (h >= h_atm[3] && h < h_atm[4]) {
    l = 3;
  } else if (h >= h_atm[5]) {
    l = 4;
  }

  if (l == 4) {
     res = b_atm[4] / c_atm[4];

  } else {

    G4double div = b_atm[l]/(c_atm[l]);
    res = div*exp(-(h)/(c_atm[l]));

  }
  return res;
}

G4double B1DetectorConstruction::get_mass_overburden(G4double h) {
  //G4cout << h << G4endl;
  G4int l;
  G4double res;
  if (h >= h_atm[0] && h < h_atm[1]) {
    l = 0;
  } else if (h >= h_atm[1] && h < h_atm[2]) {
    l = 1;
  } else if (h >= h_atm[2] && h < h_atm[3]) {
    l = 2;
  } else if (h >= h_atm[3] && h < h_atm[4]) {
    l = 3;
  } else if (h >= h_atm[5]) {
    l = 4;
  }
  if (l == 4) {
    res = (a_atm[4]) - ((b_atm[4])/(c_atm[4]))*(h);
  } else {
    res = (a_atm[l]) + (b_atm[l])*exp(-(h)/(c_atm[l]));
    G4cout << "TT" << G4endl;
    //G4cout << res << G4endl;
  }
  return res;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
