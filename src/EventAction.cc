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

const G4int nOfVolumes = 32;
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
  G4bool b1 = false, b2 = false;

  // searching of worked volumes
  for (G4int i = 0; i < nOfVolumes; i++){
    if (fEdep[i] > 500.0*keV){
      if (b1){
        n2 = i;
        b2 = true;
      }
      else{
        n1 = i;
        b1 = true;
      }
    }
  }
  
  // calculation of angle & radius
  G4double cot_theta = 0;
  G4double theta = 0;
  G4double rForSynogram = 0;
  G4double r1 = 0, r2 = 0;
  G4double theta1 = 0, theta2 = 0;
  G4double phi1 = 2 * M_PI / nOfVolumes * n1;
  G4double phi2 = 2 * M_PI / nOfVolumes * n2;

  theta = M_PI / nOfVolumes * (n2 + n1);

  rForSynogram = rOfCircle * cos((phi2 - phi1) / 2);
  if (theta > M_PI){
    theta = theta - M_PI;
    rForSynogram *= -1;
  }

  /*if (abs(sin(phi2) - sin(phi1)) < 1e-6) {
    theta = M_PI / 2;
    rForSynogram = rOfCircle * cos (phi1);
  }
  else if (n2 - n1 == 32){
    theta = phi1 + M_PI / 2;
    rForSynogram = 0;
  }
  else if (abs(cos(phi2) - cos(phi1)) < 1e-6){
    if (phi1 < M_PI / 2) theta = 0;
    else theta = M_PI;
    rForSynogram = (rOfCircle * cos (phi1 - theta));
  }
  else {
    G4double tanTheta = (cos(phi2) - cos(phi1))/(sin(phi1) - sin(phi2));
    theta  = atan(tanTheta);
    if (tanTheta < 0)
      theta += M_PI;
    rForSynogram = (rOfCircle * cos (phi1 - theta));
  }*/
  /*if (!(sin(phi1) - sin(phi2))) cos(phi1) >= cos(phi2) ? theta : theta += M_PI;
  else{
    cot_theta = (sin(phi2) - sin(phi1)) / (cos(phi1) - cos(phi2));
    if (!(cot_theta + tan(phi1))) theta = phi1;
    else cot_theta < 0 ? theta = M_PI + atan(1 / cot_theta) : theta = atan(1 / cot_theta);
    abs(cos(phi1 - theta)) != 1 ? rForSynogram = rOfCircle * abs(cos(phi1 - theta)) : rForSynogram;
  }*/
  /*if (n2 - n1 <= nOfVolumes / 2)
    rForSynogram = rOfCircle * cos( M_PI * (G4double) (n2 - n1) / (G4double) nOfVolumes);
  else {
    //theta <= M_PI ? theta += M_PI : theta -= M_PI;
    rForSynogram = rOfCircle * cos(M_PI * (1.0 - (G4double) (n2 - n1) / (G4double) nOfVolumes)); 
  }
  theta = 2 * M_PI / nOfVolumes * n1 - rForSynogram / rOfCircle;*/
  /*theta1 = theta - M_PI / 2;
  theta2 = theta + M_PI / 2;
  r1 = rForSynogram;
  r2 = rForSynogram;
  if (theta1 > M_PI / 2)
    r1 *= -1;
  if (theta2 > M_PI / 2)
    r2 *= -1;

  theta1 *= 180.0 / M_PI;
  theta2 *= 180.0 / M_PI;*/

/*G4double x1 = rOfCircle * cos(phi1), 
  y1 = rOfCircle * sin(phi1), 
  x2 = rOfCircle * cos(phi2), 
  y2 = rOfCircle * cos(phi2);

  G4double k0 = 0; 
  if (x1 - x2) k0 = (y1 - y2) / (x1 - x2);
  G4double b0 = 0;
  G4double k = 0;
  if (x1 - x2) {
    b0 = y1 - k0 * x1;
    if (!k0){
      theta = M_PI / 2;
      rForSynogram = b0;
    }
    else{
      k = - 1 / k0;
      G4double x = b0 / (k - k0);
      G4double y = k * x;
      rForSynogram = sqrt (x * x + y * y);
      theta = atan (k);
      if (k < 0)
        theta += M_PI;
      //if (b0/k0 < 0)
        //rForSynogram *= -1;
    }
  }
  else {
    rForSynogram = x1;
    if (x1 < 0)
     theta = - M_PI;
  }*/

  theta = theta * 180.0 / M_PI;

  //
  // file forming
  //

  fstream ofile;
  if (b1 && b2){
    std::string clipboard1 = to_string(rForSynogram) + "\t" + to_string(theta);
    ofile.open("dat.dat", ios::app);

    
    if(ofile.is_open()){
      //ofile << "Hello world!" << endl;
      //ofile << n1 << "; " << n2 << endl;
      ofile << clipboard1 << std::endl;
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
