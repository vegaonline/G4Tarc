/**********************************************************
 * @file      G4TARCParallelWorld.cxx
 * @author    Abhijit Bhattacharyya
 * @brief     Parallelization processes
 *********************************************************/
#include "G4TARCParallelWorld.hh"



G4TARCParallelWorld::G4TARCParallelWorld( G4String& tarcParallelWorld)
: G4VUserParallelWorld( tarcParallelWorld), fConstructed(false) {}

G4TARCParallelWorld::~G4TARCParallelWorld() {}

void G4TARCParallelWorld::Construct() {
  if (fConstructed) return;
  // A dummy material is used to fill the volulmes of the readout geometry.
  G4Material* dummyMat = nullptr;

  // declaring the readout world
  G4VPhysicalVolume* ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
  fConstructed = true;

  
}

void G4TARCParallelWorld::ConstructSD() {

}
