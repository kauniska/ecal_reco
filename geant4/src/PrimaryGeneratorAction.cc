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

#include "PrimaryGeneratorAction.hh"

#include "PrimaryGeneratorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
:Detector(det)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  SetDefaultKinematic();
  beam = 0.*mm;
  
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::GetX0() const
{
  return x0 + (Detector->GetCalorSizeXY() - 1.) / 2.; // + Detector->GetCalorThickness() * sin(theta) * cos(phi);
}

G4double PrimaryGeneratorAction::GetY0() const
{
  return y0 + (Detector->GetCalorSizeXY() - 1.) / 2.; // + Detector->GetCalorThickness() * sin(theta) * sin(phi);
}

void PrimaryGeneratorAction::SetDefaultKinematic()
{
  // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // G4String particleName;
  // G4ParticleDefinition* particle = particleTable->FindParticle(particleName="mu-");
  // particleGun->SetParticleDefinition(particle);
  // theta = 2. * M_PI / 2. * (G4UniformRand() - 0.5); // angle to vertical
  // phi = 2. * M_PI / 2. * (G4UniformRand() - 0.5);
  // G4cout << "theta " << theta << ", phi " << phi << G4endl;
  // x0 = 2 * 24 * (G4UniformRand() - 0.5); // x at the top of the ecal
  // y0 = 2 * 24 * (G4UniformRand() - 0.5); // y at the top of the ecal
  // particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), -cos(theta)));
  // particleGun->SetParticleEnergy(2.*GeV);
  // G4double position = 0.5*(Detector->GetWorldSizeZ());
  // particleGun->SetParticlePosition(G4ThreeVector(x0 * mm, y0 * mm, position));

  G4ThreeVector position = particleGun->GetParticlePosition();
  x0 = 0.; // x at the top of the ecal
  y0 = 0.; // y at the top of the ecal
  G4double z0 = position.z();
  theta = 0.; // angle to vertical
  phi = 0.;
  particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), -cos(theta)));
  particleGun->SetParticleEnergy(106 * MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //
  //randomize beam, if requested.
  //
  if (beam > 0.) 
    {
      G4ThreeVector position = particleGun->GetParticlePosition();    
      G4double maxXY = 0.49*(Detector->GetCalorSizeXY());
      x0 = 2. * maxXY * (G4UniformRand() - 0.5) / 2.; // x at the top of the ecal
      y0 = 2. * maxXY * (G4UniformRand() - 0.5) / 2.; // y at the top of the ecal
      // x0 = position.x() + (G4UniformRand() - 0.5) * beam;
      // y0 = position.y() + (G4UniformRand() - 0.5) * beam;
      G4double z0 = position.z();
      // G4double z0 = -Detector->GetCalorThickness() / 2.;
      // G4double z0 = 100. * km;
      if (std::abs(x0) > maxXY) {
        x0 = maxXY;
      }
      if (std::abs(y0) > maxXY) {
        y0 = maxXY;
      }
      theta = 2. * atan(maxXY / z0) * (G4UniformRand() - 0.5); // angle to vertical
      phi = 2. * atan(maxXY / z0) * (G4UniformRand() - 0.5);
      theta = 0.;
      phi = 0.;
      x0 = 0.;
      y0 = 5.;
      z0 = 10. * cm;
      particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
      // particleGun->SetParticlePosition(position);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), -cos(theta)));
      // particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1));
      G4double max_energy = 100; // in GeV
      G4double min_energy = 1;   // in GeV
      energy = 1.0 * GeV * exp(log(min_energy) + G4UniformRand() * (log(max_energy / min_energy)));
      // particleGun->SetParticleEnergy(energy);
      particleGun->SetParticleEnergy(106 * MeV);
      particleGun->GeneratePrimaryVertex(anEvent);
    }
    else {
      particleGun->GeneratePrimaryVertex(anEvent);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

