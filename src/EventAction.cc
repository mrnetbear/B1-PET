//ТУТ ЗАПОЛНЯЮТСЯ ГИСТОГРАММЫ ПО ИТОГАМ ОДНОГО (?) забега


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
/// \file B1/src/EventAction.cc
/// \brief Implementation of the B1::EventAction class

#define Pi 3.1415926535

#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "G4AnalysisManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


#include <iostream>
#include <fstream>

using namespace std;

const G4int nOfVolumes = 64;
const G4double rOfCircle = 63.75*mm;

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  
  for(G4int i  = 0; i < nOfVolumes; i++){
    fEdep[i] = 0.;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{
  // accumulate statistics in run action
  //fRunAction->AddEdep(fEdep);
  auto analysisManager = G4AnalysisManager::Instance();
  
  //
  // calculation of synogram coordinates
  //
  G4int n1 = -1, n2 = 0;

  // searching of worked volumes
  for (G4int i = 0; i < nOfVolumes; i++){
    if (fEdep[i] > 500.0*keV){
      if (n1 > 0)
        n2 = i;
      else
        n1 = i;
    }
  }
  
  // calculation of angle & radius
  G4double theta = Pi / nOfVolumes * (n1 + n2);
  G4double rForSynogram = rOfCircle * cos(theta - 2 * Pi * n1 / nOfVolumes);
  theta = theta * 180.0 / Pi - 90.0;

  //
  // file forming
  //

  fstream ofile;
  if (n1*n2){
    std::string clipboard = to_string(rForSynogram) + "\t" + to_string(theta);
    ofile.open("raw_synogram.txt", ios::app);

    
    if(ofile.is_open()){
      //ofile << "Hello world!" << endl;
      //ofile << n1 << "; " << n2 << endl;
      ofile << clipboard << std::endl;
    }

    ofile.close();
  }

  //fill hists
  analysisManager->FillH1(0, fEdep[0]);
  analysisManager->FillH1(1, fEdep[25]);

  //fill ntuple
  analysisManager->FillNtupleDColumn(0, fEdep[0]);
  analysisManager->FillNtupleDColumn(1, fEdep[25]);
  analysisManager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
