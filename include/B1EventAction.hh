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
/// \file B1EventAction.hh
/// \brief Definition of the B1EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <map>
#include <vector>
#include "G4ParticleDefinition.hh"
#include <iostream>
#include <fstream>

class B1RunAction;
class B1DetectorConstruction;
/// Event action class
///
using namespace std;

class B1EventAction : public G4UserEventAction
{
  public:
    B1EventAction(B1RunAction* runAction);
    virtual ~B1EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    void PrintParticles(std::map<const G4String, int>& container, std::ofstream& gamma, std::ofstream& mu_minus, std::ofstream& mu_plus,
                                       std::ofstream& nu_e, std::ofstream& nu_mu, std::ofstream& anti_nu_e, std::ofstream& anti_nu_mu, G4double position);
    void WriteHistogram(const G4String name, G4int Detec);
    void AddEdep(G4double edep) { fEdep += edep; };
    G4double SetPosition(G4int detec);

  private:
    B1RunAction* fRunAction;
    G4double     fEdep;
    std::vector<G4double> dec_pos;
    B1DetectorConstruction* detectorConstruction;
    G4double half_height;
    G4double variation;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
