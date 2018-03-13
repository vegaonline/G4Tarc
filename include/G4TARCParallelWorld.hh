/********************************************************************
 * @file     G4TARCParallelWorld.hh
 * @author   Abhijit Bhattacharyya
 * @brief    Here use parallel world
 ********************************************************************/

#ifndef G4TARC_PARALLELWORLD_H
#define G4TARC_PARALLELWORLD_H
#include "G4VUserParallelWorld.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4PVReplica.hh"
#include "G4Material.hh"
#include "G4TARCleadTargSD.hh"
#include "G4TARCVirtualSD.hh"
#include "G4SystemOfUnits.hh"


class G4LogicalVolume;
class G4VPhysicalVolume;
class G4TARCleadTargSD;
class G4TARCVirtualSD;

class G4TARCParallelWorld : public G4VUserParallelWorld {
public:
  G4TARCParallelWorld (G4String&);
  virtual ~G4TARCParallelWorld();

public:
  virtual void Construct();
  virtual void ConstructSD();

private:
  G4bool             fConstructed;
  G4double           fXblockB       = 300.0 * mm;
  G4double           fYblockB       = 300.0 * mm;
  G4double           fZblockB       = 600.0 * mm;
  G4double           fHalfXblockB   = 0.5 * fXblockB;
  G4double           fHalfYblockB   = 0.5 * fYblockB;
  G4double           fHalfZblockB   = 0.5 * fZblockB;
  G4double           fHalfXVBox     = 0.5 * 150 * mm;  // virtual box around each hole 15cms X 15cms X 300cms
  G4double           fHalfYVBox     = 0.5 * 150 * mm;
  G4double           fHalfZVBox     = 0.5 * 300 * mm;
  G4double           fZposVBox      = -1500.0 * mm + fHalfZVBox;
  G4int              fNslicesVBox   = 20; // Virtual box is divided into 20 slices
  G4double           fDiaMAxSphere  = 3300.0 * mm;
  G4double           fRadMaxSphere  = 0.5 * fDiaMAxSphere;
  G4VSolid*          LeadTargetS;
  G4VSolid*          VirtualBoxS;
  G4VSolid*          VirtualSphereS;
  G4LogicalVolume*   LeadTargetLV;
  G4LogicalVolume*   VirtualSphereLV;
  G4LogicalVolume*   VirtualBox1LV;
  G4LogicalVolume*   VirtualBox2LV;
  G4LogicalVolume*   VirtualBox3LV;
  G4LogicalVolume*   VirtualBox4LV;
  G4LogicalVolume*   VirtualBox5LV;
  G4LogicalVolume*   VirtualBox6LV;
  G4LogicalVolume*   VirtualBox7LV;
  G4LogicalVolume*   VirtualBox8LV;
  G4LogicalVolume*   VirtualBox9LV;
  G4LogicalVolume*   VirtualBox10LV;
  G4LogicalVolume*   VirtualBox11LV;
  G4LogicalVolume*   VirtualBox12LV;
  G4LogicalVolume*   ghostWorldLog;
  G4VPhysicalVolume* ghostWorld;
  G4VPhysicalVolume* VB1PV;
  G4VPhysicalVolume* VB2PV;
  G4VPhysicalVolume* VB3PV;
  G4VPhysicalVolume* VB4PV;
  G4VPhysicalVolume* VB5PV;
  G4VPhysicalVolume* VB6PV;
  G4VPhysicalVolume* VB7PV;
  G4VPhysicalVolume* VB8PV;
  G4VPhysicalVolume* VB9PV;
  G4VPhysicalVolume* VB10PV;
  G4VPhysicalVolume* VB11PV;
  G4VPhysicalVolume* VB12PV;

};

#endif
