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
#include "G4TARCEventAction.hh"

class G4TARCEventAction;
class G4TARCParallelWorld;

class G4TARCSteppingAction : public G4UserSteppingAction {
public:
  G4TARCSteppingAction(G4TARCEventAction*);
  virtual ~G4TARCSteppingAction(){};

  virtual void UserSteppingAction(const G4Step*);
  void ProcessStepping(const G4Step*);


private:
  //G4TARCHistoManager*  fHisto;
  G4TARCEventAction*   fEventAction;
  G4bool flag;
  G4int number_generations, fNmax;

  G4double fEnergy0, fTime0;

  std::map<G4int, G4double, std::less<G4int> > fParentEnergy;
  std::map<G4int, G4String, std::less<G4int> > fParentParticle;
  std::map<G4int, G4int, std::less<G4int> > fParentParticleID;
};

#endif
