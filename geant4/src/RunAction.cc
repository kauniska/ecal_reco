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
/// \file electromagnetic/TestEm5/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "Run.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction *det, EventAction *eventAction, PrimaryGeneratorAction *kin)
    : G4UserRunAction(), fDetector(det), fPrimary(kin), fRun(0), fEventAction(eventAction), fFileName("data")
{
  // Book predefined histograms
  // fHistoManager = new HistoManager(fEventAction);
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetFileName("data");
  analysisManager->SetVerboseLevel(1);
  // analysisManager->SetNtupleMerging(true);

  // create ROOT tree
  analysisManager->CreateNtuple("events", "recorded info per event");
  analysisManager->CreateNtupleDColumn(0, "E");                                           // 0
  analysisManager->CreateNtupleDColumn(0, "Edep");                                        // 1
  analysisManager->CreateNtupleIColumn(0, "pdg", fEventAction->GetEcalPDG());                      // 2
  analysisManager->CreateNtupleDColumn(0, "EcalEdep", fEventAction->GetEcalEdep()); // 3
  analysisManager->CreateNtupleIColumn(0, "layerID", fEventAction->GetEcalLayers()); // 4
  // analysisManager->CreateNtupleIColumn(0, "barID", fEventAction->GetEcalBars()); // 5
  analysisManager->CreateNtupleIColumn(0, "barID", fEventAction->GetEcalCopyNo()); // 5
  analysisManager->CreateNtupleIColumn(0, "Nhits", fEventAction->GetEcalHits()); // 6
  // analysisManager->CreateNtupleIColumn(0, "copyNo", fEventAction->GetEcalCopyNo()); // 7

  analysisManager->FinishNtuple(0);

  analysisManager->SetActivation(false); // enable inactivation of histograms

  // Define histograms start values
  const G4int kMaxHisto = 6;
  const G4String id[] = {"dummy", "Ecal_tot", "Evis_tot", "Etot_profile", "Evis_rofile", "Evis_scint"};
  const G4String title[] =
      {
          "dummy",              // 0
          "total Etot in Ecal", // 1
          "total Evis in Ecal", // 2
          "Etot profile",       // 3
          "Evis profile",       // 4
          "Evis per scint"      // 5
      };

  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated
  // as we have not yet set nbins, vmin, vmax
  for (G4int k = 0; k < kMaxHisto; k++)
  {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ 
  // delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(fDetector);
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  
  // save Rndm status
  ////  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  if (isMaster) G4Random::showEngineStatus();
     
  // keep run condition
  if ( fPrimary ) { 
    G4ParticleDefinition* particle 
      = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(particle, energy);
  }
  
  //histograms
  //        
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // if ( analysisManager->IsActive() ) {
    
  //   analysisManager->OpenFile();
  // }
  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  // print Run summary
  //
  if (isMaster) fRun->EndOfRun();    
      
  // save histograms
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  // if ( analysisManager->IsActive() ) {    
  //   analysisManager->Write();
  //   analysisManager->CloseFile();
  // }  

  // show Rndm status
  if (isMaster) G4Random::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
