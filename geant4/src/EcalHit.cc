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
/// \file EcalHit.cc
/// \brief Implementation of the B5::EcalHit class

#include "EcalHit.hh"
#include "DetectorConstruction.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<EcalHit> *EcalHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EcalHit::EcalHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EcalHit::EcalHit(G4int barID,G4int layerID)
: fBarID(barID), fLayerID(layerID)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EcalHit::~EcalHit()
{}

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

G4int EcalHit::GetPDG() const
{
  if (fPDs.size() == 0) {
    return -2; // if no particle interaction
  }
  bool there_is_muon = false;
  bool there_is_electron = false;
  for (auto pd : fPDs) {
    if (pd != nullptr) {
      if (pd->GetParticleDefinitionID() == 407 || pd->GetParticleDefinitionID() == 406) {
        there_is_muon = true;
      } else if (pd->GetParticleDefinitionID() == 324 || pd->GetParticleDefinitionID() == 323) {
        there_is_electron = true;
      }
    }
  }
  if (there_is_muon && there_is_electron) {
    return -5; // special code
  } else if (there_is_muon) {
    return 407;
  } else if (there_is_electron) {
    return 324;
  } else {
    return fPDs[0]->GetParticleDefinitionID();
  }
  // G4int tmp = -1;
  // if (fPD != nullptr) {
    // G4cout << "ooo " << fPD->GetParticleDefinitionID() << " " << fPD->GetParticleName() << G4endl;
    // G4cout << fPD->GetParticleName() << G4endl;
    // return fPD->GetParticleDefinitionID();
  // } else {
    // return -2;
  // }
  // if (fPD->GetParticleName() != "e-")  G4cout << fPD->GetParticleName() << G4endl;
  // if (tmp == 407) { // mu-
  //   G4cout << "here" << G4endl;
  //   return 1;
  // } else if (tmp == 324){ //e-
  //   return 2;
  // } else if (tmp == 323) { //e+
  //   return 3;
  // } else if (tmp == 344) { // gamma
  //   return 4;
  // } else if (tmp == 428) { //pi0
  //   return 5;
  // } else if (tmp == 426) { // pi+
  //   return 6;
  // } else if (tmp == 427) { // pi-
  //   return 7;
  // } else if (tmp == 406) { // mu+
  //   return 8;
  // } else { // anything else
  //   G4cout << "UNKOWN PARTICLE " << fPD->GetParticleDefinitionID() << " " << fPD->GetParticleName() << G4endl;
  //   return 9;
  // }
}

    // //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    // EcalHit::EcalHit(const EcalHit &right)
    // : G4VHit(),
    //   fBarID(right.fBarID),
    //   fLayerID(right.fLayerID),
    //   fEdep(right.fEdep),
    //   fPos(right.fPos),
    //   fRot(right.fRot)
    // {}

    // //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    // const EcalHit& EcalHit::operator=(
    //         const EcalHit &right)
    // {
    //   fBarID = right.fBarID;
    //   fLayerID = right.fLayerID;
    //   fEdep = right.fEdep;
    //   fPos = right.fPos;
    //   fRot = right.fRot;
    //   return *this;
    // }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool EcalHit::operator==(const EcalHit &right) const
{
  return ( fBarID==right.fBarID && fLayerID==right.fLayerID );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EcalHit::Draw()
{
  auto visManager = G4VVisManager::GetConcreteInstance();
  if (! visManager || (fEdep==0.)) return;

  // Draw a calorimeter cell with depth propotional to the energy deposition
  G4Transform3D trans(fRot.inverse(),fPos);
  G4VisAttributes attribs;
  G4Colour colour(1.,0.,0.);
  attribs.SetColour(colour);
  attribs.SetForceSolid(true);
  G4Box box("dummy",15.*cm,15.*cm,1.*m*fEdep/(0.1*GeV));
  visManager->Draw(box,attribs,trans);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String,G4AttDef>* EcalHit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("EcalHit",isNew);

  if (isNew) {
    (*store)["HitType"]
      = G4AttDef("HitType","Hit Type","Physics","","G4String");

    (*store)["Bar"]
      = G4AttDef("Bar","Bar ID","Physics","","G4int");

    (*store)["Row"]
      = G4AttDef("Row","Row ID","Physics","","G4int");

    (*store)["Energy"]
      = G4AttDef("Energy","Energy Deposited","Physics","G4BestUnit",
                 "G4double");

    (*store)["Pos"]
      = G4AttDef("Pos", "Position", "Physics","G4BestUnit",
                 "G4ThreeVector");
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* EcalHit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;

  values
    ->push_back(G4AttValue("HitType","EcalHit",""));
  values
    ->push_back(G4AttValue("Bar",G4UIcommand::ConvertToString(fBarID),
                           ""));
  values
    ->push_back(G4AttValue("Row",G4UIcommand::ConvertToString(fLayerID),""));
  values
    ->push_back(G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),""));
  values
    ->push_back(G4AttValue("Pos",G4BestUnit(fPos,"Length"),""));

  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EcalHit::Print()
{
  G4cout << "  Cell[" << fLayerID << ", " << fBarID << "] "
    << fEdep/MeV << " (MeV) " << fPos << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......