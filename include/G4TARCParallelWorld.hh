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
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class G4TARCParallelWorld : public G4VUserParallelWorld {
public:
  G4TARCParallelWorld (G4String&);
  virtual ~G4TARCParallelWorld();

public:
  virtual void Construct();
  virtual void ConstructSD();

private:
  G4bool fConstructed;
};

#endif
