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

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "DetectorConstruction.hh"
#include "G4SDManager.hh"
#include "Constants.hh"

#include <vector>
#include <array>
#include <map>

class PrimaryGeneratorAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// named constants
const G4int kHad = 1;

class EventAction : public G4UserEventAction
{
  public:  
    EventAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);
    
    void SumDeStep(G4int, G4int, G4int, G4double);
	
	  void WriteScints(const G4Event*);
    std::vector<G4double>& GetEcalEdep() {return fEcalEdep;}
    std::vector<G4int>& GetEcalHits() {return fEcalHits;}
    std::vector<G4int>& GetEcalLayers() {return fEcalLayerID;}
    std::vector<G4int>& GetEcalBars() {return fEcalBarID; }
    std::vector<G4int>& GetEcalPDG() {return fEcalPDG;}
    std::vector<G4int>& GetEcalCopyNo() {return fEcalCopyNo;}
    std::vector<G4String>& GetEcalParticleNames() {return fEcalParticleNames;}
    void AddNsec(const G4int& n) {fNsec += n;}
    G4int GetNsec() {return fNsec;}
    void SetProcessID(const G4int&);
    G4int GetProcessID() {return fProcessID;};
    void SetDecayPosition(const G4ThreeVector& pos) {fPosDecay = pos;}

    std::vector<G4double>& GetElectronEnergies() { return fElectronEnergies;}
    std::vector<G4double>& GetVertexEnergies() { return fVertexEnergies;}
    std::vector<G4int> &GetElectronLayers() { return fElectronLayers; }
    std::vector<G4int> &GetElectronRows() { return fElectronRows; }
    std::vector<G4int> &GetMuonIDs() { return fMuonIDs; }
    void SetElectronEnergy(G4int parentID, G4double E);
    void SetVertexEnergy(G4int parentID, G4double E);
    void SetElectronLayer(G4int parentID, G4int layer);
    void SetElectronRow(G4int parentID, G4int row);
    void SetMuonID(G4int ID);

  private:
    DetectorConstruction *detector;
    PrimaryGeneratorAction* primary;
	
	  G4int nbOfModules, nbOfLayers, kLayerMax;     
    std::vector<G4double>   EtotLayer;
    std::vector<G4double>   EvisLayer;
	
	  G4double EtotCalor;
	  G4double EvisCalor;
    G4int fNsec = 0;
	
	  std::map<G4int, G4double> EvisScint;

    // hit collections Ids
    G4int fEcalHCID = -1;
    // energy deposit in calorimeters cells
    std::vector<G4double> fEcalEdep = std::vector<G4double>(int(kNofEcalCells), 0.);
    std::vector<G4double> fEcalE = std::vector<G4double>(int(kNofEcalCells), 0.);
    std::vector<G4int> fEcalHits = std::vector<G4int>(int(kNofEcalCells), 0.); // number of hits
    std::vector<G4int> fEcalLayerID = std::vector<G4int>(int(kNofEcalCells), -1);
    std::vector<G4int> fEcalBarID = std::vector<G4int>(int(kNofEcalCells), -1);
    std::vector<G4int> fEcalPDG = std::vector<G4int>(int(kNofEcalCells), -1);
    std::vector<G4int> fEcalCopyNo = std::vector<G4int>(int(kNofEcalCells), -1);
    std::vector<G4double> fElectronEnergies;// = std::vector<G4double>(int(kNofSecParticles), -1);
    std::vector<G4double> fVertexEnergies; // = std::vector<G4double>(int(kNofSecParticles), -1);
    std::vector<G4int> fElectronLayers;// = std::vector<G4int>(int(kNofEcalCells), -1);
    std::vector<G4int> fElectronRows;// = std::vector<G4int>(int(kNofEcalCells), -1);
    std::vector<G4int> fMuonIDs;// = std::vector<G4int>(int(kNofEcalCells), -1);
    std::vector<G4String> fEcalParticleNames = std::vector<G4String>(int(kNofEcalCells), "-");
    std::vector<std::vector<G4double> > fEcalEdepMatrix = std::vector<std::vector<G4double> >(int(kNofEcalLayers), std::vector<G4double>(int(kNofEcalBars), -1));
    std::vector<std::vector<G4int> > fEcalHitsMatrix = std::vector<std::vector<G4int> >(int(kNofEcalLayers), std::vector<G4int>(int(kNofEcalBars), -1));
    G4int fProcessID = -1;
    G4ThreeVector fPosDecay;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
