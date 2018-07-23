/********************************************************************
 * @file     G4TARCParallelWorld.hh
 * @author   Abhijit Bhattacharyya
 * @brief    Here use parallel world
 ********************************************************************/

#ifndef G4TARC_PARALLELWORLD_H
#define G4TARC_PARALLELWORLD_H

#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <string>

#include "G4TARCDetectorConstruction.hh"
#include "G4TARCPVolumeStore.hh"
#include "G4TARCVolumeSD.hh"

// primitive scorers related
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSNofCollision.hh"
#include "G4PSPopulation.hh"
#include "G4PSTrackCounter.hh"
#include "G4PSTrackLength.hh"

//importance biasing related
#include "G4IStore.hh"

#include "G4GeometryCell.hh"

#include "G4VUserParallelWorld.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
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
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4SolidStore.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VIStore;
class G4TARCVolumeSD;

class G4TARCParallelWorld : public G4VUserParallelWorld {
public:
  G4TARCParallelWorld (G4String&);
  virtual ~G4TARCParallelWorld();

public:
  virtual void Construct();
  virtual void ConstructSD();
  virtual void ConstructSDandField(){};
  void DefineShellsBlocks();
  inline G4VPhysicalVolume& GetShellVolumeAddress() const { return *fShellPhys; }
  inline G4VPhysicalVolume& GetWorldVolumeAddress() const { return *ghostWorld; }
  inline  G4VPhysicalVolume* GetWorldVolume() {return ghostWorld;}
  G4GeometryCell GetGeometryCell(G4int);
  G4String GetCellName(G4int);
  G4VIStore* CreateImportanceStore();


private:
    G4TARCPVolumeStore                fPVolumeStore;
    G4MultiFunctionalDetector*        fTARCNeutronDet;
    G4MultiFunctionalDetector*        fTARCProtonDet;
    G4int                             fShellNumber;
    G4double                          fShellThickness;
    G4double                          fMaxOuterRadiusofShell;
    G4double                          fMinInnerRadiusofShell;
    G4double                          fInnerRadProtonShell;
    G4double                          fOuterRadProtonShell;
    std::vector<G4double>             fOuterRadiusofShell;
    std::vector<G4double>             fInnerRadiusofShell;
    std::vector<G4LogicalVolume*>     fLVvector;
    G4bool                            fConstructed;
    G4double                          fHalfXBlockB;
    G4double                          fHalfYBlockB;
    G4double                          fHalfZBlockB;
    G4double                          fHalfXVBox;
    G4double                          fHalfYVBox;
    G4double                          fHalfZVBox;
    G4double                          fNewHalfZProt;
    G4double                          fZposProt;
    G4double                          fDiaMaxSphere;
    G4double                          fRadMaxSphere;
    G4LogicalVolume*                  fVBoxLogProton;
    G4VPhysicalVolume*                fVBoxPVProton;
    G4LogicalVolume*                  fShellLog;
    G4VPhysicalVolume*                fShellPhys;
    G4LogicalVolume*                  ghostWorldLog;
    G4VPhysicalVolume*                ghostWorld;


    G4LogicalVolumeStore*             fLVS;
    G4PhysicalVolumeStore*            fPVS;
    std::vector<G4VPhysicalVolume*>::const_iterator fPVciter;
    std::vector<G4LogicalVolume*>  ::const_iterator fLVciter;
};

#endif
