#ifndef LEADTARGSD_H
#define LEADTARGSD_H

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4VTouchable.hh"
#include "G4ios.hh"
#include "G4TARCHistoManager.hh"
#include "G4ParticleDefinition.hh"

class G4Step;
class G4TouchableHistory;
class G4HCofThisEvent;
class G4TARCHistoManager;

class G4TARCleadTargSD: public G4VSensitiveDetector {
public:
  G4TARCleadTargSD(G4String&);
  virtual ~G4TARCleadTargSD(){};

  virtual void Initialize(G4HCofThisEvent*){};
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*){};

private:
  G4TARCHistoManager* fHisto;
};


#endif
