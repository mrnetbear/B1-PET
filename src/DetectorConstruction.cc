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
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#define Pi 3.1415926535

#include <vector>
#include <string>
#include <math.h>

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4SubtractionSolid.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"

const G4double rOfCircle = 63.75*mm;
const G4int nOfDetectors = 32;

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  for (G4int i = 0; i < nOfDetectors; i++){
    fScoringVolume[i] = nullptr;
  }

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 60*cm, env_sizeZ = 60*cm;
  //G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;
  G4bool fCheckOverlaps = true;

  G4double a;  // mass of a mole; 
  G4double density; 
  G4double z;  // z=mean number of protons; 
  G4int natoms;
   
   
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);
  G4Element* C = new G4Element("Carboneum", "C", z=6 , a=12.01*g/mole);
  G4Element* O = new G4Element("Oxygen" , "O", z=8 , a=16.00*g/mole);

  //Gd3Al2Ga3O12 6.63 г/см3
  
  G4Element* Gd = new G4Element("Gadoliniumd", "Gd", z=64 , a=157.25*g/mole);
  G4Element* Al = new G4Element("Aluminium", "Al", z=13 , a=26.98*g/mole); 
  G4Element* Ga = new G4Element("Galium", "Ga", z=31 , a=69.72*g/mole);  

  G4double nGd=3, nAl=2, nGa=3, nO=12; 
   
  G4Material* GAGG = new G4Material("GAGG", density= 6.63*g/cm3, 4);
   
  GAGG->AddElement(Gd, natoms=nGd);
  GAGG->AddElement(Al, natoms=nAl);
  GAGG->AddElement(Ga, natoms=nGa);  
  GAGG->AddElement(O, natoms=nO);

  // Scintillator GAGG
  G4double photonEnergy[] = {2.757*eV, 8.819*eV};

  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

  //Lu2SiO5 (LSO) 7.40 г/см3
  
  G4Element* Lu = new G4Element("Lutetium", "Lu", z=71 , a=174.96*g/mole);
  G4Element* Si = new G4Element("Silicium", "Si", z=14 , a=28.085*g/mole);

  G4double nlLu=2, nlSi=1, nlO=5;

  G4Material* LSO = new G4Material("LSO", density= 7.4*g/cm3, 3);

  LSO->AddElement(Lu, natoms=nlLu);
  LSO->AddElement(Si, natoms=nlSi);
  LSO->AddElement(O, natoms=nlO);

//Bi4Ge3O12 (BGO) 7.13 г/см3
  
  G4Element* Bi = new G4Element("Bismuthum", "Bi", z=83 , a=208.98*g/mole);
  G4Element* Ge = new G4Element("Germаnium", "Ge", z=32 , a=72.63*g/mole);

  G4double nbBi=4, nbGe=3, nbO=12;

  G4Material* BGO = new G4Material("BGO", density= 7.13*g/cm3, 3);

  BGO->AddElement(Bi, natoms=nbBi);
  BGO->AddElement(Ge, natoms=nbGe);
  BGO->AddElement(O, natoms=nbO);

//
// GAGG
  

  G4double refractiveIndex1[] = {1.9, 1.9};

  G4double absorption[] = {1.3*cm,  1.3*cm};

  G4double scintilFast[] = {0.19, 0.28};

  G4double scintilSlow[] = {0.93, 1.00};

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1, nEntries);
  myMPT1->AddProperty("ABSLENGTH",    photonEnergy, absorption, nEntries);
  myMPT1->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintilFast, nEntries);
  myMPT1->AddProperty("SCINTILLATIONCOMPONENT2", photonEnergy, scintilSlow, nEntries);
        
  myMPT1->AddConstProperty("SCINTILLATIONYIELD",0./MeV); 
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 92*ns);
  myMPT1->AddConstProperty("SCINTILLATIONTIMECONSTANT2",630*ns);
  myMPT1->AddConstProperty("SCINTILLATIONYIELD1",0.18);

  GAGG->SetMaterialPropertiesTable(myMPT1); 

//
// LSO

  G4double refractiveIndexLu[] = {1.82, 1.82};
  G4double absorptionLu[] = {1.14*cm, 1.14*cm};
  G4double scintilFastLu[] = {0.19, 0.28};
  
  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();

  myMPT2->AddProperty("RINDEX",       photonEnergy, refractiveIndexLu, nEntries);
  myMPT2->AddProperty("ABSLENGTH",    photonEnergy, absorptionLu, nEntries);
  myMPT2->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintilFastLu, nEntries);
        
  myMPT2->AddConstProperty("SCINTILLATIONYIELD",0./MeV); 
  myMPT2->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT2->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 40*ns);
  myMPT2->AddConstProperty("SCINTILLATIONYIELD1",0.18);

  LSO->SetMaterialPropertiesTable(myMPT2); 

//
// BGO

  G4double refractiveIndexBi[] = {2.15, 2.15};
  G4double absorptionBi[] = {1.118*cm, 1.118*cm};
  G4double scintilFastBi[] = {0.19, 0.28};
  
  G4MaterialPropertiesTable* myMPT3 = new G4MaterialPropertiesTable();

  myMPT3->AddProperty("RINDEX",       photonEnergy, refractiveIndexBi, nEntries);
  myMPT3->AddProperty("ABSLENGTH",    photonEnergy, absorptionBi, nEntries);
  myMPT3->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintilFastBi, nEntries);
        
  myMPT3->AddConstProperty("SCINTILLATIONYIELD",0./MeV); 
  myMPT3->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT3->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 300*ns);
  myMPT3->AddConstProperty("SCINTILLATIONYIELD1",0.18);

  BGO->SetMaterialPropertiesTable(myMPT3); 
        
//
// SiPM features description
//
        
  G4Material* SiPM_mat = nist->FindOrBuildMaterial("G4_Si");
              
  G4double refractiveIndexSiPM[] = {1.4626, 1.4626};


  G4double absorptionSiPM[] = {100.*nm,  100.*nm};
  G4double efficiencySiPM[] = {0.23, 0.23};


  G4MaterialPropertiesTable* myMPT4 = new G4MaterialPropertiesTable();

  myMPT4->AddProperty("RINDEX", photonEnergy, refractiveIndexSiPM,nEntries, true);
  myMPT4->AddProperty("ABSLENGTH", photonEnergy, absorptionSiPM, nEntries, true);    
  myMPT4->AddProperty("EFFICIENCY",photonEnergy, efficiencySiPM, nEntries, true);

  SiPM_mat->SetMaterialPropertiesTable(myMPT4);


  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  //world_mat->SetMaterialPropertiesTable(myMPT2);

  auto solidWorld = new G4Box("World",                           // its name
    0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking

 
  //
  //Shape GAGG crystals
  //

  G4double shapeGAGG_dx = 0.15*cm;
  G4double shapeGAGG_dy = 0.15*cm;
  G4double shapeGAGG_dz = 1.0*cm;
  
  std::string shape;
  G4Box* solidShapeSc[nOfDetectors];
  G4LogicalVolume* logicShapeSc[nOfDetectors];
  G4VPhysicalVolume* physVolumeSc[nOfDetectors];

//creation of crystals

  for (int i = 0; i < nOfDetectors; i++){
    shape = "Shape" + std::to_string(i);
    G4Box* solidShape = new G4Box(shape, shapeGAGG_dx, shapeGAGG_dy, shapeGAGG_dz);
    G4LogicalVolume* logicShape = new G4LogicalVolume(solidShape, BGO, shape);
    double angle = 2 * Pi * i /nOfDetectors;
    double posX = 0;
    double posY = sin(angle) * rOfCircle;
    double posZ = cos(angle) * rOfCircle;
    G4RotationMatrix* xRot = new G4RotationMatrix;
    xRot -> rotateX(angle);
    G4VPhysicalVolume* physVolume = new G4PVPlacement(xRot, G4ThreeVector(posX, posY, posZ), logicShape, shape, logicWorld, false, 0, checkOverlaps);
    solidShapeSc[i] = solidShape;
    logicShapeSc[i] = logicShape;
    physVolumeSc[i] = physVolume;
  }

  //
  // Set GAGG crystal No.0 as scoring volume
  //
  fScoringVolume[0] = logicShapeSc[0];

  for (G4int i = 0; i < nOfDetectors; i++){
    fScoringVolume[i] = logicShapeSc[i];
  }

  G4VisAttributes* CrystallVisAtt= new G4VisAttributes(G4Colour(255.0,255.0,0.0));
    CrystallVisAtt->SetVisibility(true);
    for (int  i = 0; i < nOfDetectors; i++){
      logicShapeSc[i]->SetVisAttributes(CrystallVisAtt);
    }

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
