#include "G4TARCStackingAction.hh"

G4TARCStackingAction::G4TARCStackingAction()
: G4UserStackingAction(), fHistoManager(0), fStackMessenger(0), fParticle(0) {
  fStackMessenger = new G4TARCStackingMessenger(this);
  fHistoManager   = G4TARCHistoManager::GetPointer();
  fKillSecondary  = false;
  fParticle       = 0;
}

G4TARCStackingAction::~G4TARCStackingAction() {
  delete fStackMessenger;
}

G4ClassificationOfNewTrack G4TARCStackingAction::ClassfyNewTrack(const G4Track* myTrack) {
  G4ClassificationOfNewTrack status = fUrgent;

  if (myTrack->GetTrackStatus() == fAlive)
    fHistoManager->ScoreNewTrack(myTrack);

  const G4ParticleDefinition* part = myTrack->GetDefinition();

  if (fHistoManager->GetVerbose() > 1) {
    G4cout << " Track #"                      << myTrack->GetTrackID()
           << " of "                          << part->GetParticleName()
           << " E (MeV)= "                    << myTrack->GetKineticEnergy() / MeV
           << " produced by Track ID= "       << myTrack->GetParentID()
           << G4endl;
  }
  if (myTrack->GetTrackID() > 1) {
    if (fKillSecondary || part == fParticle)
      status = fKill;
  }
  return status;
}

void G4TARCStackingAction::SetKillStatus( G4bool value ) {
  fKillSecondary = value;
}

void G4TARCStackingAction::SetKill( const G4String& name ) {
  fParticle = G4ParticleTable::GetParticleTable()->FindParticle(name);
}
