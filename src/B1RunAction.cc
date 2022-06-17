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
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1SteppingAction.hh"
//adicionando Tracking Action
#include "B1TrackingAction.hh"
// #include "B1Run.hh"
#include "B1Analysis.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{
  // add new units for dose
  //
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{

//  createHistogram();
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values

  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

  const G4ParticleGun* particleGun = generatorAction->GetParticleGun();

  G4double primaryEnergy = particleGun->GetParticleEnergy()/GeV;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();

  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }

  // Print
  //

  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }

  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : "
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

//   WriteHistogram();
//  delete G4AnalysisManager::Instance();
  //printando partículas criadas
/*
  std::map<const G4ParticleDefinition*, int>&
    particlesCreatedInWorld = fptrackingAction->GetNParticlesCreatedInWorld();

  G4cout << "Partículas criadas dentro do mundo :" << G4endl;

  PrintParticles(particlesCreatedInWorld);

  G4cout << "_______________________" << G4endl;

  std::map<const G4ParticleDefinition*, int>&
    particlesCreatedOutsideWorld = fptrackingAction->GetNParticlesCreatedOutsideWorld();

  G4cout << "Prtículas criadas fora do mundo :" << G4endl;
  PrintParticles(particlesCreatedOutsideWorld);
*/
   //limpa dados
   fptrackingAction->clearParticles();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

void B1RunAction::createHistogram()
{
  //Avisa que está criando
  /*
  G4cout << "CRIANDO HISTOGRAMAS" << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //mostrando o tipo
  G4cout << "utilizando " << analysisManager->GetType() << "Analysis Manager" << G4endl;
  //setando o verbose
  analysisManager->SetVerboseLevel(1);

  //Abrindo arquivo de saída
  G4String filename = "air_shower";
  analysisManager->OpenFile(filename);

  analysisManager->CreateNtuple("air_shower", "physics");
  analysisManager->CreateNtupleDColumn("mu_minus");
  analysisManager->CreateNtupleDColumn("mu_plus");
  analysisManager->CreateNtupleDColumn("gamma");
  analysisManager->FinishNtuple();
  */

  //criando ntuples
  /*
  analysisManager->CreateNtuple("air_shower", "physics");
	analysisManager->CreateNtupleDColumn("proton");
	analysisManager->CreateNtupleDColumn("e_plus");
	analysisManager->CreateNtupleDColumn("e_minus");
	analysisManager->CreateNtupleDColumn("mu_minus");
	analysisManager->CreateNtupleDColumn("mu_plus");
	analysisManager->CreateNtupleDColumn("pi_plus");
	analysisManager->CreateNtupleDColumn("pi_minus");
	analysisManager->CreateNtupleDColumn("gamma");
  analysisManager->CreateNtupleDColumn("z");
	analysisManager->FinishNtuple();
*/
}

void B1RunAction::WriteHistogram()
{
  /*
  G4cout<<"ESCREVENDO"<<G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //salvando ghistogramas

  analysisManager->Write();
  analysisManager->CloseFile();
 */
}

//função para printar partículas
void B1RunAction::PrintParticles(std::map<const G4ParticleDefinition*, int>& container) {
  /*
  std::map<const G4ParticleDefinition*, int>::iterator it;
    for(it = container.begin() ;it != container.end(); it ++)
    {
      G4cout << "N " << it->first->GetParticleName() << " : " << it->second << G4endl;
    }
    */
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
