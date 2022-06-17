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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"
#include "B1Hits.hh"
#include "B1DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "B1Analysis.hh"

#include <iostream>
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event* event)
{
  fEdep = 0.;



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* event)
{
  // accumulate statistics in run action
  /*
  fRunAction->AddEdep(fEdep);

// pega as collections ID's
 G4SDManager * SDman = G4SDManager::GetSDMpointer();
 G4HCofThisEvent* HCE = event->GetHCofThisEvent();

 detectorConstruction = (B1DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
 G4int num = detectorConstruction->return_num_detec();
 dec_pos = detectorConstruction->GetPositions();
  half_height = detectorConstruction->get_size();
  variation = detectorConstruction->get_delta();

  std::ofstream gammafile("gamma.txt");
  std::ofstream mu_plusfile("mu+.txt");
  std::ofstream mu_minusfile("mu-.txt");
  std::ofstream nu_e_file("nu_e.txt");
  std::ofstream nu_mu_file("nu_mu.txt");
  std::ofstream anti_nu_e_file("anti_nu_e.txt");
  std::ofstream anti_nu_mu_file("anti_nu_mu.txt");

 std::vector<G4int> col;
 std::vector<B1HitsCollection*> HitsCol;

 col.reserve(num);
 HitsCol.reserve(num);

 for(G4int i=0; i < num; i++) {
   col[i] = SDman->GetCollectionID("SD" + std::to_string(i + 1));

   if(HCE) {
     HitsCol[i] = (B1HitsCollection*)(HCE->GetHC(col[i]));
   }

   if(HitsCol[i]) {
     G4double pos_detec = SetPosition(i + 1);
     int n_hit = HitsCol[i]->entries();
     G4cout << "My detector " + std::to_string(i+1) +" has " << n_hit << "hits" << G4endl;
     B1Hits* hit = new B1Hits;
     std::map<const G4String, int> fparticles;
     for(int i1 = 0; i1 < n_hit; i1++) {
      B1Hits* hit = (*HitsCol[i])[i1];
      const G4String name = hit->getParticleInTarget();
      //G4cout << name << G4endl;
      fparticles[name]++;
      //WriteHistogram(name, 1);
    }
    PrintParticles(fparticles, gammafile, mu_plusfile, mu_minusfile, nu_e_file, nu_mu_file, anti_nu_e_file, anti_nu_mu_file, pos_detec);
    fparticles.clear();
   }
   }

   gammafile.close();
   mu_plusfile.close();
   mu_minusfile.close();
   nu_e_file.close();
   nu_mu_file.close();
   anti_nu_e_file.close();
   anti_nu_mu_file.close();
 }


void B1EventAction::PrintParticles(std::map<const G4String, int>& container, std::ofstream& gamma, std::ofstream& mu_minus, std::ofstream& mu_plus,
                                   std::ofstream& nu_e, std::ofstream& nu_mu, std::ofstream& anti_nu_e, std::ofstream& anti_nu_mu, G4double position) {
  std::map<const G4String, int>::iterator it;
//  G4cout << "Número de párticulas identificadas no detetor: " << G4endl;
   G4bool isGamma, ismu_minus, ismu_plus, isnu_e, isnu_mu, isanti_nu_e, isanti_nu_mu = false;

    for(it = container.begin() ;it != container.end(); it ++)
    {
      G4cout << "N " << it->first << " : " << it->second << G4endl;
      //FUTURO -> substituir por switch;
      if(it->first == "gamma") {
        //Escreve o valor no arquivo
        gamma << it->second << " " << position << "\n";
        isGamma = true;
      } else if (it->first == "mu+") {
        mu_plus << it->second << " " << position << "\n";
        ismu_plus= true;
      } else if (it->first == "mu-") {
        mu_minus << it->second << " " << position << "\n";
        ismu_minus = true;
      } else if (it->first == "nu_e") {
        nu_e << it->second << " " << position << "\n";
        isnu_e = true;
      } else if (it->first == "nu_mu") {
        nu_mu << it->second << " " << position << "\n";
        isnu_mu = true;
      } else if (it->first == "anti_nu_e") {
        anti_nu_e << it->second << " " << position << "\n";
        isanti_nu_e = true;
      } else if (it->first == "anti_nu_mu") {
        anti_nu_mu << it->second << " " << position << "\n";
        isanti_nu_mu = true;
      }
  }

  if (!isGamma) {
    gamma << "0 " << position << "\n";
  }
  if (!ismu_plus) {
    mu_plus << "0 " << position << "\n";
  }
  if (!ismu_minus) {
    mu_minus << "0 " << position << "\n";
  }
  if (!isnu_e) {
    nu_e << "0 " << position << "\n";
  }
  if (!isnu_mu) {
    nu_mu << "0 " << position << "\n";
  }
  if (!isanti_nu_e) {
    anti_nu_e << "0 " << position << "\n";
  }
  if (!isanti_nu_mu) {
    anti_nu_mu << "0 " << position << "\n";
  }
}

void B1EventAction::WriteHistogram(const G4String name, G4int Detec) {
  /*
   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if (name == "gamma") {
    analysisManager->FillNtupleDColumn(2, Detec);
  } else if (name == "mu+") {
    analysisManager->FillNtupleDColumn(1, Detec);
  } else if (name == "mu-") {
    analysisManager->FillNtupleDColumn(0, Detec);
  }
  analysisManager->AddNtupleRow();
*/
}

G4double B1EventAction::SetPosition(G4int detec) {
  /*
  //coloca a origem no topo da atmosfera
  G4int detec_index = detec - 1;
  G4double first_value = dec_pos[0]/1000000;
  G4cout << first_value << G4endl;
  G4double value =  dec_pos[detec_index];
  G4double value_c = value/1000000;
  G4double half_height_real = half_height/1000000;
  G4double h_real = half_height_real*2;
  G4double real_delta = variation/1000000;
  G4double first_distance = half_height_real + first_value ;
  if (detec_index == 0) {
    return first_distance;
  }
  return (first_distance + detec_index*real_delta);
  */
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
