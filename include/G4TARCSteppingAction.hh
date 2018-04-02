#ifndef G4TARCSteppingAction_HH
#define G4TARCSteppingAction_HH

#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "globals.hh"
#include <map>

#include "G4TARCParallelWorld.hh"

class G4TARCEventAction;
class G4TARCParallelWorld;

class G4TARCSteppingAction : public G4UserSteppingAction {
public:
  G4TARCSteppingAction(G4TARCEventAction*);
  virtual ~G4TARCSteppingAction(){};

  //virtual void UserSteppingAction(const G4Step*);
public:
  G4double GetPrimaryEnergy() { return startEnergy; }
  G4double GetPrimaryTime() { return startTime; }
  void DefineShellsBlocks();
  // virtual void UserSteppingAction(const G4Step*);
private:
  G4double startEnergy;
  G4double startTime;
  G4double fInnershellRadius;

  G4DataVector  radii;

  G4TARCEventAction*   fEventAction;
  G4bool flag;
  std::map<G4int, G4double, std::less<G4int> > parent_energy;
  std::map<G4int, G4String, std::less<G4int> > parent_particle;
  std::map<G4int, G4int, std::less<G4int> > parent_particleID;
  G4int number_generations;
  G4double                              fHalfXBlockB;
  G4double                              fHalfYBlockB;
  G4double                              fHalfZBlockB;
  G4double                              fHalfXVBox;
  G4double                              fHalfYVBox;
  G4double                              fHalfZVBox;
  G4double                              fNewHalfZProt;
  G4double                              fZposProt;
  G4double                              fDiaMAxSphere;
  G4double                              fRadMaxSphere;
  G4int                                 fShellNumber;
  G4double                              fShellThickness;
  G4double                              fMaxOuterRadiusofShell;
  G4double                              fMinInnerRadiusofShell;
  G4double                              fInnerRadProtonShell;
  G4double                              fOuterRadProtonShell;
  std::vector<G4double>                 fOuterRadiusofShell;
  std::vector<G4double>                 fInnerRadiusofShell;

};

#endif
