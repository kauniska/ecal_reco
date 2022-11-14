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

#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "EcalSD.hh"
#include "G4VisAttributes.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:scintMat(0),lvol_scint(0), absorberMat(0),lvol_layer(0),
 moduleMat(0),lvol_module(0), calorimeterMat(0),lvol_calorimeter(0),
 worldMat(0),pvol_world(0), defaultMat(0)
{
  // materials
  DefineMaterials();
  
  // default parameter values of calorimeter
  //
  scintDiameter       = 1.13*mm; 	//1.08*mm
  nbOfScints          = 24;		//490
  scintWidth          = 16.*mm;	//1.35*mm
  scintHeight         = 4.*mm;
  scintLength         = 384. * mm;
  gapSize             = 1. * mm;                // 662.175*mm
  nbOfLayers          = 1;		    //10
  nbOfModules         = 16;		    //9
  leadThickness       = 4.*mm;
  aluThickness        = 1.*mm;
  layerThickness      = scintHeight + 2*aluThickness + leadThickness; // 1.68*mm
  shieldThickness     = 10. *mm;
  moduleThickness     = layerThickness;
  calorThickness      = moduleThickness * nbOfModules;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // define Elements
  //
  G4Element* H  = new G4Element("Hydrogen","H", 1,  1.01*g/mole);
  G4Element* C  = new G4Element("Carbon",  "C", 6, 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen","N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",  "O", 8, 16.00*g/mole);

  G4int natoms, ncomponents;
  G4double density, massfraction;				     

  // Lead
  //
  G4Material* Pb =   
  new G4Material("Lead", 82., 207.20*g/mole, density= 0.98*11.20*g/cm3);

  //Aluminium
  //
  G4Material* Al =
      new G4Material("aluminium", 13.,26.98 * g / mole, density = 2.7 * g / cm3);

  // Scintillator
  //
  G4Material* Sci = 
  new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(C, natoms=8);
  Sci->AddElement(H, natoms=8);
  
  Sci->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  // Air
  //
  G4Material* Air = 
  new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, massfraction=70*perCent);
  Air->AddElement(O, massfraction=30.*perCent);

  // example of vacuum
  //
  density     = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  G4Material* Vacuum =   
  new G4Material("Galactic", 1., 1.008*g/mole, density,
                             kStateGas,temperature,pressure);

  //attribute materials
  //
  defaultMat     = Vacuum;  
  scintMat       = Sci;
  absorberMat    = Pb;
  moduleMat      = defaultMat;
  calorimeterMat = defaultMat;
  worldMat       = defaultMat;
  shieldingMat   = Al;

  // print table
  //      
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{
  // Cleanup old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // scints
  //
  G4Box* svol_scint = new G4Box("scint",			//name
                         0.5 * scintWidth, 0.5 * scintLength, 0.5 * scintHeight);
			 
  lvol_scint = new G4LogicalVolume(svol_scint,		//solid
                                   scintMat,		//material
                                   "scint");		//name

  // Lead
  //
  G4double sizeZ = leadThickness;
  G4double sizeY = scintLength + gapSize;
  G4double sizeX = scintLength + gapSize;

  G4Box* svol_lead = new G4Box("lead",                                // name
                             0.5 * sizeX, 0.5 * sizeY, 0.5 * sizeZ); // size

  lvol_lead = new G4LogicalVolume(svol_lead,  // solid
                                   absorberMat, // material
                                   "lead");    // name

  // Al
  //
  sizeZ = aluThickness;
  sizeY = scintLength + gapSize;
  sizeX = scintLength + gapSize;

  G4Box* svol_Al_layer = new G4Box("Al_layer",                                // name
                             0.5 * sizeX, 0.5 * sizeY, 0.5 * sizeZ); // size

  lvol_Al_layer = new G4LogicalVolume(svol_Al_layer,  // solid
                                   absorberMat, // material
                                   "Al_layer");    // name

  // layer
  //
  sizeZ = layerThickness;
  sizeY = scintWidth*nbOfScints;
  sizeX = scintLength;
  
  G4Box* svol_layer = new G4Box("layer",			//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size


  lvol_layer = new G4LogicalVolume(svol_layer,		//solid
                                   defaultMat,		//material
                                   "layer");		//name

  // put scints within layer
  //
  G4double Zcenter = -0.5 * layerThickness;
  G4double Xcenter = -0.5 * (sizeX - scintWidth);

  for (G4int k = 0; k < nbOfScints; k++)
  {
    G4ThreeVector pos = G4ThreeVector(Xcenter + k * scintWidth, 0., Zcenter + 0.5 * scintHeight);
    new G4PVPlacement(0,
                      pos,        // rotation+position
                      lvol_scint, // logical volume
                      "scint",    // name
                      lvol_layer, // mother
                      false,      // no boulean operat
                      k + 1);
  }

  Zcenter += scintHeight;
  // put alu 1 in layer
  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., Zcenter + 0.5 * aluThickness), // rotation+position
                    lvol_Al_layer,                      // logical volume
                    "Al_layer",                         // name
                    lvol_layer,                     // mother
                    false,                          // no boulean operat
                    0);

  // put Lead in layer
  Zcenter += aluThickness;
  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., Zcenter + 0.5 * leadThickness), // rotation+position
                    lvol_lead,                                       // logical volume
                    "lead",                                          // name
                    lvol_layer,                                      // mother
                    false,                                           // no boulean operat
                    0);

  // put Aluminium 2 in layer
  Zcenter += leadThickness;
  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., Zcenter + 0.5 * aluThickness), // rotation+position
                    lvol_Al_layer,                  // logical volume
                    "Al_layer",                     // name
                    lvol_layer,                     // mother
                    false,                          // no boulean operat
                    0);

  // modules
  //
  moduleThickness = layerThickness*nbOfLayers;       
  sizeZ = moduleThickness;
  sizeY = scintLength;
  sizeX = scintLength;
  
  G4Box* svol_module = new G4Box("module",			//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size

  lvol_module = new G4LogicalVolume(svol_module,	//solid
                                   absorberMat,		//material
                                   "module");		//name

  // put layers within module
  //
  Zcenter = -0.5*(nbOfLayers+1)*layerThickness;
  G4double Ycenter =  0.25*scintWidth;
  
  for (G4int k=0; k<nbOfLayers; k++) {
    Zcenter += layerThickness;
    Ycenter  = - Ycenter;
    new G4PVPlacement(0,		   		//no rotation
      		  G4ThreeVector(0.,Ycenter,Zcenter),    //position
                      lvol_layer,     		   	//logical volume	
                      "layer",	   			//name
                      lvol_module,        		//mother
                      false,             		//no boulean operat
                      k+1);               		//copy number

  }

  // create top and bottom shielding plates
  sizeY = scintLength + shieldThickness;
  sizeX = scintLength + shieldThickness;
  G4Box *svol_topPlate = new G4Box("topPlate", 0.5 * sizeX, 0.5 * sizeY, 0.5 * shieldThickness);
  lvol_topPlate = new G4LogicalVolume(svol_topPlate, shieldingMat, "topPlate");

  // create side shielding plates
  sizeZ = calorThickness;
  sizeX = scintLength + shieldThickness;
  G4Box* svol_sidePlateXZ = new G4Box("sidePlateXZ", 0.5 * sizeX, 0.5 * shieldThickness, 0.5 * sizeZ);
  lvol_sidePlateXZ = new G4LogicalVolume(svol_sidePlateXZ, shieldingMat, "sidePlateXZ");
  sizeY = scintLength;
  G4Box* svol_sidePlateYZ = new G4Box("sidePlateYZ", 0.5 * shieldThickness, 0.5 * sizeY, 0.5 * sizeZ);
  lvol_sidePlateYZ = new G4LogicalVolume(svol_sidePlateYZ, shieldingMat, "sidePlateYZ");

  // calorimeter
  sizeZ = calorThickness + shieldThickness;
  sizeY = scintLength + shieldThickness;
  sizeX = scintLength + shieldThickness;
  
  G4Box* svol_calorimeter = new G4Box("calorimeter",		//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size


  lvol_calorimeter = new G4LogicalVolume(svol_calorimeter,	//solid
                                   calorimeterMat,		//material
                                   "calorimeter");		//name  

  // put modules inside calorimeter
  Zcenter = -0.5*(calorThickness + moduleThickness);
  for (G4int k=0; k<nbOfModules; k++) {
    Zcenter += moduleThickness;		  
    G4RotationMatrix rotm;                    //rotation matrix to place modules    
    if ((k+1)%2 == 0) rotm.rotateZ(90*deg);
	  G4Transform3D transform(rotm, G4ThreeVector(0.,0.,Zcenter));    
    new G4PVPlacement(transform,		   		//rotation+position
                      lvol_module,	     		//logical volume	
                      "module", 	   		    //name
                      lvol_calorimeter,        	//mother
                      false,             		//no boulean operat
                      k+1);               		//copy number
  }
  // put top and bottom shielding plates
  Zcenter = 0.5 * (calorThickness + shieldThickness);
  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., Zcenter), // top
                    lvol_topPlate,
                    "topPlate",
                    lvol_calorimeter,
                    false,
                    0);
  new G4PVPlacement(0,
                    G4ThreeVector(0., 0., -Zcenter), // bottom
                    lvol_topPlate,
                    "topPlate",
                    lvol_calorimeter,
                    false,
                    1);

  // put in XZ sides plates
  Ycenter = 0.5 * (scintLength + gapSize + shieldThickness);
  new G4PVPlacement(0,
                    G4ThreeVector(0., Ycenter, 0.), // bottom
                    lvol_sidePlateXZ,
                    "sidePlateXZ",
                    lvol_calorimeter,
                    false,
                    0);
  new G4PVPlacement(0,
                    G4ThreeVector(0., -Ycenter, 0.), // bottom
                    lvol_sidePlateXZ,
                    "sidePlateXZ",
                    lvol_calorimeter,
                    false,
                    1);

  // put in YZ sides plates
  Xcenter = 0.5 * (scintLength + gapSize + shieldThickness);
  new G4PVPlacement(0,
                    G4ThreeVector(Xcenter, 0., 0.), // bottom
                    lvol_sidePlateYZ,
                    "sidePlateYZ",
                    lvol_calorimeter,
                    false,
                    0);
  new G4PVPlacement(0,
                    G4ThreeVector(-Xcenter, 0., 0.), // bottom
                    lvol_sidePlateYZ,
                    "sidePlateYZ",
                    lvol_calorimeter,
                    false,
                    1);

  // world
  sizeZ = 2.*calorThickness;
  sizeY = 2.*scintLength;
  sizeX = 2.*scintLength;
  
  worldSizeZ = sizeZ;
  
  G4Box*      
  svol_world = new G4Box("world",			//name
                  0.5*sizeX, 0.5*sizeY, 0.5*sizeZ);	//size

  lvol_world = new G4LogicalVolume(svol_world,		//solid
                                   worldMat,		//material
                                   "world");		//name 
				    
  pvol_world = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 lvol_world,		//logical volume
                                 "world",		//name
                                 0,			//mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  //put calorimeter in world
  //  
  new G4PVPlacement(0,				//no rotation
                    G4ThreeVector(),		//at (0,0,0)
                    lvol_calorimeter,		//logical volume
                    "calorimeter",		//name
                    lvol_world,			//mother  volume
                    false,			//no boolean operation
                    0);				//copy number
		    				 
  PrintCalorParameters();
  
  // Visualization attributes
  //
  lvol_scint->SetVisAttributes (G4VisAttributes::GetInvisible());  
  lvol_layer->SetVisAttributes (G4VisAttributes::GetInvisible());
  lvol_world->SetVisAttributes (G4VisAttributes::GetInvisible());
    
  //always return the physical World
  //
  return pvol_world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4UnitsTable.hh"

void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n-------------------------------------------------------------"
     << "\n ---> The calorimeter is " << nbOfModules << " Modules"
     << "\n ---> A Module is " << nbOfLayers << " Layers + 1 milled Layer";
     
  G4cout  
     << "\n ---> A Layer is " << G4BestUnit(layerThickness,"Length")  
     << " thickness of " << absorberMat->GetName();    
     
  G4cout 
     << "\n ---> A Layer includes " << nbOfScints << " scints of " 
     << scintMat->GetName();
     
  G4cout 
     << "\n      ---> diameter : " << G4BestUnit(scintDiameter,"Length")
     << "\n      ---> length   : " << G4BestUnit(scintLength,"Length")
     << "\n      ---> distance : " << G4BestUnit(scintWidth,"Length");
     
  G4cout 
   << "\n\n ---> Module thickness " << G4BestUnit(moduleThickness,"Length");
  
  G4cout 
   << "\n\n ---> Total calor thickness " << G4BestUnit(calorThickness,"Length")
   <<   "\n      Tranverse size        " << G4BestUnit(scintLength,"Length");

  G4cout << "\n-------------------------------------------------------------\n";
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

void DetectorConstruction::ConstructSDandField()
{
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg =
        new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );
        
    }
    // sensitive detectors -----------------------------------------------------
    auto sdManager = G4SDManager::GetSDMpointer();
    G4String SDname;
    auto ecal = new EcalSD(SDname = "/EcalSD");
    sdManager->AddNewDetector(ecal);
    // lvol_scint->SetSensitiveDetector(ecal);
    lvol_scint->SetSensitiveDetector(ecal);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
