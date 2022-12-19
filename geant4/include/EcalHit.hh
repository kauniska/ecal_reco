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
/// \file EcalHit.hh
/// \brief Definition of the EcalHit class

#ifndef EcahlHit_h
#define EcahlHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ParticleDefinition.hh"

class G4AttDef;
class G4AttValue;

/// Hadron Calorimeter hit
///
/// It records:
/// - the cell column ID and row ID
/// - the energy deposit
/// - the cell position and rotation

class EcalHit : public G4VHit
{
  public:
    EcalHit();
    EcalHit(G4int iCol,G4int iRow);
    EcalHit(const EcalHit &right) = default;
    ~EcalHit() override;

    EcalHit& operator=(const EcalHit &right) = default;
    G4bool operator==(const EcalHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    void Draw() override;
    const std::map<G4String,G4AttDef>* GetAttDefs() const override;
    std::vector<G4AttValue>* CreateAttValues() const override;
    void Print() override;

    void SetBarID(G4int z) { fBarID = z; } // before it was column
    G4int GetBarID() const { return fBarID; }

    void SetLayerID(G4int z) { fLayerID = z; } // before it was row
    G4int GetLayerID() const { return fLayerID; }

    void SetEdep(G4double de) { fEdep = de; }
    void AddEdep(G4double de) { fEdep += de;}
    G4double GetEdep() const { return fEdep; }

    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }

    void AddNHit(G4int nhit) {fHits += nhit;}
    G4int GetHits() const {return fHits;}

    void SetPD(const G4ParticleDefinition *pd) { fPDs.push_back(pd); }
    G4int GetPDG() const;
    void SetCopyNo(G4int cNo) {fCopyNo = cNo;}
    G4int GetCopyNo() {return fCopyNo;}

  private:
    G4int fBarID = -1; // used to be Column, along x and y
    G4int fLayerID = -1; // used to be Row, along z (vertical)
    G4double fEdep = 0.;
    G4ThreeVector fPos;
    G4RotationMatrix fRot;
    G4int fHits = 0;
    std::vector<const G4ParticleDefinition*> fPDs;
    G4int fCopyNo = -1;
};

using EcalHitsCollection = G4THitsCollection<EcalHit>;

extern G4ThreadLocal G4Allocator<EcalHit>* EcalHitAllocator;

inline void* EcalHit::operator new(size_t)
{
  if (!EcalHitAllocator) {
       EcalHitAllocator = new G4Allocator<EcalHit>;
  }
  return (void*)EcalHitAllocator->MallocSingle();
}

inline void EcalHit::operator delete(void* aHit)
{
  EcalHitAllocator->FreeSingle((EcalHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
