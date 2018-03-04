#include "G4TARCCheckVolumeSD.hh"

G4TARCCheckVolumeSD::G4TARCCheckVolumeSD(const G4String& name)
: G4VSensitiveDetector(name), fHisto(0) {
  fHisto = G4TARCHistoManager::GetPointer();
}

G4bool G4TARCCheckVolumeSD::ProcessHits(G4Step* myStep, G4TouchableHistory*){
  const G4Track* myTrack = myStep->GetTrack();

  fHisto->AddLeakingParticle(myTrack);

  if (myTrack->GetTrackID() <= 1){
    fHisto->AddNzero(myTrack, myStep);
  } else {
    fHisto->NeutFinalState(myTrack, myStep);
    fHisto->TargetProfile(myTrack, myStep);
  }
  return true;
}
