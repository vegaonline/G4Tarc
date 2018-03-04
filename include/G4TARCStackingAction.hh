#ifndef G4TARCStackingAction_H
#define G4TARCStackingAction_H

#include "G4UserStackingAction.hh"
#include "globals.hh"
#include "G4TARCHistoManager.hh"
#include "G4TARCStackingMessenger.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

class G4TARCHistoManager;
class G4TARCStackingMessenger;
class G4Track;
class G4ParticleDefinition;

class G4TARCStackingAction : public G4UserStackingAction {
public:
  G4TARCStackingAction();
  virtual ~G4TARCStackingAction();

  void SetKillStatus(G4bool);
  void SetKill(const G4String&);
  virtual G4ClassificationOfNewTrack ClassfyNewTrack(const G4Track*);

private:
  G4TARCHistoManager*         fHistoManager;
  G4TARCStackingMessenger*    fStackMessenger;
  G4bool                      fKillSecondary;

  const G4ParticleDefinition* fParticle;
};


#endif
