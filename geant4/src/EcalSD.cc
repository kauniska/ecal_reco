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
/// \file EcalSD.cc
/// \brief Implementation of the EcalSD class

#include "EcalSD.hh"
#include "EcalHit.hh"
#include "Constants.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EcalSD::EcalSD(G4String name)
: G4VSensitiveDetector(name)
{
  collectionName.insert("EcalCalorimeterLayer");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EcalSD::~EcalSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EcalSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new EcalHitsCollection(SensitiveDetectorName, collectionName[0]);
  if (fHCID<0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);

  // fill calorimeter hits with zero energy deposition
  for (auto layer = 0; layer < kNofEcalLayers; layer ++) {
    for (auto bar = 0; bar < kNofEcalBars; bar++) {
      fHitsCollection->insert(new EcalHit());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool EcalSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;

  const G4VTouchable* touchable = step->GetPreStepPoint()->GetTouchable();
  auto physical = touchable->GetVolume();
  auto copyNo = physical->GetCopyNo();
  G4int layerNo = touchable->GetCopyNumber(2);
  G4int barNo = touchable->GetCopyNumber(3);
  G4int hitID = kNofEcalLayers * layerNo + barNo;
  EcalHit* hit = (*fHitsCollection)[hitID];
  const G4ParticleDefinition* pd = step->GetTrack()->GetDefinition();
  // G4cout << "ooo " << touchable->GetVolume()->GetName() << G4endl;

  // check if it is first touch
  if (hit->GetBarID()<0) {
    hit->SetBarID(barNo);
    hit->SetLayerID(layerNo);
    auto depth = touchable->GetHistory()->GetDepth();
    auto transform = touchable->GetHistory()->GetTransform(depth-2);
    transform.Invert();
    hit->SetRot(transform.NetRotation());
    hit->SetPos(transform.NetTranslation());
  }
  // add energy deposition
  hit->AddEdep(edep);
  if (step->GetTrack()->GetParticleDefinition()->GetParticleName() == "e-") {
    hit->AddElectronEdep(edep);
  } else if (step->GetTrack()->GetParticleDefinition()->GetParticleName() == "mu-") {
    hit->AddMuonEdep(edep);
  }
  hit->AddNHit(1);
  hit->SetPD(pd);
  hit->SetCopyNo(copyNo);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......