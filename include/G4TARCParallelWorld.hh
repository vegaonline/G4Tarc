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

class G4LogicalVolume;
class G4VPhysicalVolume;

class G4TARCParallelWorld : public G4VUserParallelWorld {
public:
  G4TARCParallelWorld (G4String& tarcParallelWorld);
  virtual ~G4TARCParallelWorld();

public:
  virtual void Construct();

private:
  G4bool fConstructed;
};

#endif
