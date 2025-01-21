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
/// \file B1/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#define Pi 3.1415926535

#include <time.h>
#include <math.h>
#include <random>

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GenericMessenger.hh"
#include "G4RandomDirection.hh"

const G4double innerRadius = 20*mm;
const G4double randRMin = 0.*mm;
const G4double randRMax = 5.*mm;
const G4double randAMin = 0.*mm;
const G4double randAMax = 2 * Pi *mm;


using namespace std;

double randRadius(){
  srand(time(NULL));
  return (double) rand() / (double) RAND_MAX * (randRMax - randRMin) + randRMin;
}

double randAngle(){
  srand(time(NULL));
  return (double) rand() / (double) RAND_MAX * (randAMax - randAMin) + randAMin;
}

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(): G4VUserPrimaryGeneratorAction()
{
  // define new partical gun
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  fParticleGun1 = new G4ParticleGun(n_particle);

  // define temporary particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="gamma");

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun1->SetParticleDefinition(particle);

  //
  // define particle kinematic
  //
  
  // define generator position

  G4double x0 = 0*mm;
  G4double y0 = 20*mm;
  G4double z0 = 0*mm;

  

  /*G4int sourcePos = 9; 

  G4double angle = 2 * Pi / 10 * sourcePos;

  G4double x0 = 0;
  G4double y0 = innerRadius * sin(angle);
  G4double z0 = innerRadius * cos(angle);*/


  /*G4double r = 0, phi = 0;

  r = 5 - 5 * G4UniformRand();
  phi = 2 * Pi * G4UniformRand();

  G4double x0 = 0;
  G4double y0 = r * sin(phi) + 20;
  G4double z0 = r * cos(phi) + 20;
  cout << "(" << x0 << ", " << y0 << ", " << z0 << ")" << endl;*/

  // set particle properties
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->SetParticleEnergy(0.511*MeV);

  fParticleGun1->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun1->SetParticleEnergy(0.511*MeV);
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of each event
  //
  // define particle momentum vector

  //G4Random::setTheSeed(time(NULL));

  //random point in circle

  /*G4double r = 0, phi = 0;

  r = 5 - 5 * G4UniformRand();
  phi = 2 * Pi * G4UniformRand();

  G4double sPosY = 20, sPosZ = 20;

  G4double x0 = 0;
  G4double y0 = r * sin(phi) + sPosY;
  G4double z0 = r * cos(phi) + sPosZ;*/
  //cout << "(" << x0 << ", " << y0 << ", " << z0 << ")" << endl;

  G4double x = 0;
  G4double angle = 2 * Pi * G4UniformRand();
  G4double angle1 = 0;

  if (angle >= Pi)
    angle1 = angle - Pi;
  else
    angle1 = angle + Pi;

  G4ThreeVector momentumUnitVector = G4ThreeVector(0., sin(angle), cos(angle));
  G4ThreeVector momentumUnitVector1 = G4ThreeVector(0., sin(angle1), cos(angle1));

  //G4ThreeVector momentumUnitVector = G4ThreeVector(0., 0., 1);
  //G4ThreeVector momentumUnitVector1 = G4ThreeVector(0., 0., -1);


  //fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->SetParticleMomentumDirection(momentumUnitVector);
  fParticleGun->GeneratePrimaryVertex(anEvent);

  //fParticleGun1->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun1->SetParticleMomentumDirection(momentumUnitVector1);
  fParticleGun1->GeneratePrimaryVertex(anEvent);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}


