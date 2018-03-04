/**********************************************************
 * @file      G4TARCParallelWorld.cxx
 * @author    Abhijit Bhattacharyya
 * @brief     Parallelization processes
 *********************************************************/
#include "G4TARCParallelWorld.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

G4TARCParallelWorld::G4TARCParallelWorld( G4String& tarcParallelWorld)
: G4VUserParallelWorld( tarcParallelWorld), fConstructed(false) {}

G4TARCParallelWorld::~G4TARCParallelWorld() {}

void G4TARCParallelWorld::Construct() {
  if (fConstructed) return;
  fConstructed = true;

  G4VPhysicalVolume* ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
}
