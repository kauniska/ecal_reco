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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4Track.hh"
#include "G4Positron.hh"
#include "G4VProcess.hh"
#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, EventAction* eA)
:G4UserTrackingAction(),detector(det), fEventAction(eA)
{ }
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  
  // if decay, we're done
  if (track->GetCreatorProcess() != nullptr) {
    if (track->GetCreatorProcess()->GetProcessName() == "Decay") {
      fEventAction->SetProcessID(0);
      fEventAction->SetDecayPosition(track->GetPosition());
      if (track->GetParticleDefinition()->GetParticleName() == "e-")
      {
        fEventAction->SetElectronEnergy(track->GetParentID(), track->GetTotalEnergy() / MeV);
        fEventAction->SetVertexEnergy(track->GetParentID(), track->GetVertexKineticEnergy() / MeV);
      }
      return;
    }
  }
  // we only consider muons and neutral pions
  if (track->GetParticleDefinition()->GetParticleName() == "mu-"\
    || track->GetParticleDefinition()->GetParticleName() == "mu+"\
    || track->GetParticleDefinition()->GetParticleName() == "pi0")
    //|| track->GetParticleDefinition()->GetParticleName() == "e-"\
    //|| track->GetParticleDefinition()->GetParticleName() == "e+")
  {
    const G4VProcess *proc = track->GetCreatorProcess();
    if (proc != nullptr) // nullptr if created by gun
    {
      G4String name = proc->GetProcessName();
          // In order: electron ionisation, muon ionisation, electron Bremsstrahlung, photoelectric effect, compton effect
          G4int proc_id = -1;
      if (name == "muPairProd")
      {
        proc_id = 1;
      } else if (name == "annihil") {
        proc_id = 2;
      } else if (name == "conv") {
        proc_id = 3;
      } else if (name == "muIoni") {
        proc_id = 4;
      } else {
        G4cout << "Unidentified process : " << name << " " << proc->GetProcessSubType() << G4endl;
        proc_id = 5;
      }
      // add to event
      if (proc_id != -1) {
        G4cout << "Identified process : " << name << " " << proc->GetProcessSubType() << G4endl;
        fEventAction->SetProcessID(proc_id);
      }
    }
  }
  // G4cout << track->GetCreatorProcess()->GetProcessName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{ 
  //compute leakage
  //
  // if not at World boundaries, return
  if (track->GetVolume() != detector->GetPvolWorld()) return;
     
  //get position
  G4double x = (track->GetPosition()).x();
  G4double xlimit = 0.5*(detector->GetCalorThickness());
  G4int icase = 2;
  if (x >= xlimit) icase = 0;
  else if (x <= -xlimit) icase = 1;
  
  //get particle energy
  G4double Eleak = track->GetKineticEnergy();
  if (track->GetDefinition() == G4Positron::Positron())
     Eleak += 2*electron_mass_c2;
     
  //sum leakage
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->DetailedLeakage(icase,Eleak);         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

