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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "EcalHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Utility function which finds a hit collection with the given Id
// and print warnings if not found
G4VHitsCollection *GetHC(const G4Event *event, G4int collId)
{
  auto hce = event->GetHCofThisEvent();
  if (!hce)
  {
    G4ExceptionDescription msg;
    msg << "No hits collection of this event found." << G4endl;
    G4Exception("EventAction::EndOfEventAction()",
                "Code001", JustWarning, msg);
    return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if (!hc)
  {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl;
    G4Exception("EventAction::EndOfEventAction()",
                "Code001", JustWarning, msg);
  }
  return hc;
}

EventAction::EventAction(DetectorConstruction* det,PrimaryGeneratorAction* prim)
:detector(det), primary(prim), fNsec(0)
{
  nbOfModules = detector->GetNbModules();	 	
  nbOfLayers  = detector->GetNbLayers();
  kLayerMax = nbOfModules*nbOfLayers + 1;
  
  EtotCalor = EvisCalor = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

void EventAction::SetProcessID(const G4int &p)
{
  if (fProcessID != -1 && fProcessID != p)
  {
    G4cout << "ovewriting process ID from " << fProcessID << " to " << p << G4endl;
  }
  fProcessID = p;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  EtotLayer.resize(kLayerMax);
  EvisLayer.resize(kLayerMax);			
  for (G4int k=0; k<kLayerMax; k++) {
    EtotLayer[k] = EvisLayer[k] = 0.0;
  }
  EtotCalor = EvisCalor = 0.;
  EvisScint.clear();

  // Find hit collections and histogram Ids by names (just once)
  // and save them in the data members of this class
  if (fEcalHCID == -1)
  {
    auto sdManager = G4SDManager::GetSDMpointer();
    // auto analysisManager = G4AnalysisManager::Instance();
    // hits collections name
    G4String EcalName = "EcalSD/EcalCalorimeterLayer";


    fEcalHCID = sdManager->GetCollectionID(EcalName);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumDeStep(G4int iModule, G4int iLayer, G4int iScint,
                            G4double deStep )
{
  if (iModule > 0) EtotCalor += deStep;
  		
  G4int kLayer = 0; G4int kScint = 0;
  if (iLayer > 0) {
	kLayer = (iModule-1)*nbOfLayers + iLayer;
	EtotLayer[kLayer] += deStep;
  }
  
  if (iScint > 0) {
	EvisLayer[kLayer] += deStep;
	EvisCalor += deStep;
    kScint = 1000*kLayer + iScint;
	EvisScint[kScint] += deStep;	  	
  }	  
}

void EventAction::SetElectronEnergy(G4int parentID, G4double E)
{
  int index = -1;
  for (int i = 0; i < fMuonIDs.size(); i++) {
    if (parentID == fMuonIDs[i]) {
      index = i; // we don't do anything if muon already in list
    }
  }
  if (index == -1) {
    fMuonIDs.push_back(parentID);
    fElectronEnergies.push_back(E);
  } else {
    fElectronEnergies[index] = E;
  }
}

void EventAction::SetVertexEnergy(G4int parentID, G4double E)
{
  int index = -1;
  for (int i = 0; i < fMuonIDs.size(); i++)
  {
    if (parentID == fMuonIDs[i])
    {
      index = i; // we don't do anything if muon already in list
    }
  }
  if (index == -1)
  {
    fMuonIDs.push_back(parentID);
    fVertexEnergies.push_back(E);
  }
  else if (index >= fVertexEnergies.size()) {
    fVertexEnergies.push_back(E);
  } else  {
    fVertexEnergies[index] = E;
  }
}

void EventAction::SetElectronLayer(G4int parentID, G4int layer)
{}
void EventAction::SetElectronRow(G4int parentID, G4int row)
{}
void EventAction::SetMuonID(G4int ID)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  //pass informations to Run
  //
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  	
  for (G4int k=0; k<kLayerMax; k++) {
     run->SumEvents_1(k,EtotLayer[k],EvisLayer[k]);   
  }
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
  analysisManager->FillH1(1,EtotCalor / MeV);
  analysisManager->FillH1(2,EvisCalor / MeV);

  G4double Ebeam = primary->GetParticleGun()->GetParticleEnergy();
  G4double Eleak = Ebeam - EtotCalor;
  run->SumEvents_2(EtotCalor,EvisCalor,Eleak);
  
  
  std::map<G4int,G4double>::iterator it;         
  for (it = EvisScint.begin(); it != EvisScint.end(); it++) {
     G4int kScint = it->first;
	 G4int iScint = kScint%1000;
     G4double Evis = it->second;
	 analysisManager->FillH1(5,iScint+0.5,Evis / MeV);
  }

  // Ecalorimeters hits
  G4int totalCalHit = 0;
  G4double totalCalEdep = 0.;

   auto hc = GetHC(event, fEcalHCID);
   if (!hc) return;

  totalCalHit = 0;
  totalCalEdep = 0.;
  for (unsigned long i = 0; i < hc->GetSize(); ++i) {
    G4double edep = 0.;

    EcalHit* hit = static_cast<EcalHit *>(hc->GetHit(i));
    edep = hit->GetEdep();

    if (edep > 0.) {
      totalCalHit++;
      totalCalEdep += edep;
    }

    fEcalEdep[i] = edep / MeV;
    fEdepElectron[i] = hit->GetElectronEdep() / MeV;
    fEdepMuon[i] = hit->GetMuonEdep() / MeV;
    fEcalHits[i] = hit->GetHits();
    fEcalLayerID[i] = hit->GetLayerID();
    fEcalBarID[i] = hit->GetBarID();
    fEcalPDG[i] = hit->GetPDG();
    fEcalCopyNo[i] = hit->GetCopyNo();
  }

  for (unsigned long i = 0; i < hc->GetSize(); ++i) {
    G4double edep = 0.;
    G4int Nhits = 0;
    G4int l = -1;
    G4int r = -1;
    auto hit = static_cast<EcalHit *>(hc->GetHit(i));
    edep = hit->GetEdep();
    Nhits = hit->GetHits();
    l = hit->GetBarID(); // layer ID
    r = hit->GetLayerID();

    if (edep > 0. && l > -1 && r > -1) {
      fEcalEdepMatrix[l][r] += edep;
      fEcalHitsMatrix[l][r] += Nhits;
    }
  }

  analysisManager->FillNtupleDColumn(1, 0, primary->GetEnergy() / MeV);
  analysisManager->FillNtupleDColumn(1, 1, totalCalEdep / MeV);
  analysisManager->FillNtupleDColumn(1, 2, primary->GetTx());
  analysisManager->FillNtupleDColumn(1, 3, primary->GetTy());
  analysisManager->FillNtupleDColumn(1, 4, primary->GetX0() / cm);
  analysisManager->FillNtupleDColumn(1, 5, primary->GetY0() / cm);
  analysisManager->FillNtupleIColumn(1, 6, fNsec);
  analysisManager->FillNtupleIColumn(1, 7, fProcessID);
  analysisManager->FillNtupleDColumn(1, 8, fPosDecay[0] / cm);
  analysisManager->FillNtupleDColumn(1, 9, fPosDecay[1] / cm);
  analysisManager->FillNtupleDColumn(1, 10, fPosDecay[2] / cm);
  analysisManager->FillNtupleDColumn(1, 11, primary->GetSizeXY() / cm);
  analysisManager->FillNtupleDColumn(1, 12, primary->GetSizeZ() / cm);
  analysisManager->AddNtupleRow(0);
  analysisManager->AddNtupleRow(1);
  analysisManager->AddNtupleRow(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
        
#include <fstream>

void EventAction::WriteScints(const G4Event* evt)
{
  // event is appended on a file
  //
  G4String name = G4AnalysisManager::Instance()->GetFileName();
  G4String fileName = name + ".scints.ascii";
  
  std::ofstream File(fileName, std::ios::app);
  std::ios::fmtflags mode = File.flags();  
  File.setf( std::ios::scientific, std::ios::floatfield );
  G4int prec = File.precision(3);
    
  //write event number  
  //
  File << evt->GetEventID() << G4endl;
  
  //gun particle informations
  //
  G4ParticleGun* gun = primary->GetParticleGun();
  G4double ekin = gun->GetParticleEnergy();
  G4ThreeVector direction = gun->GetParticleMomentumDirection();
  G4ThreeVector position  = gun->GetParticlePosition();
  File << ekin << " " << direction << " " << position << G4endl;  
  
  //write scints
  //
  File << EvisScint.size() << G4endl;
  //
  std::map<G4int,G4double>::iterator it;         
  for (it = EvisScint.begin(); it != EvisScint.end(); it++) {
     G4int kScint = it->first;
     G4double Evis = it->second;
     File << " " << std::setw(7) << kScint << " "<< std::setw(10) << Evis
            << G4endl;
  }
           
  File << G4endl;
    
  // restaure default formats
  File.setf(mode,std::ios::floatfield);
  File.precision(prec);         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

