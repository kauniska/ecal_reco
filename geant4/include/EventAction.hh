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
const G4int kDim = 2;

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

  private:  
    DetectorConstruction*   detector;
    PrimaryGeneratorAction* primary;
	
	  G4int nbOfModules, nbOfLayers, kLayerMax;     
    std::vector<G4double>   EtotLayer;
    std::vector<G4double>   EvisLayer;
	
	  G4double EtotCalor;
	  G4double EvisCalor;
	
	  std::map<G4int, G4double> EvisScint;

    // hit collections Ids
    std::array<G4int, kDim> fEcalHCID = {-1, -1};
    // energy deposit in calorimeters cells
    std::vector<G4double> fEcalEdep = std::vector<G4double>(int(kNofEcalCells), 0.);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
