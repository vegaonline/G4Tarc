#include "G4TARCVirtualSD.hh"

G4TARCVirtualSD::G4TARCVirtualSD(G4String& name)
: G4VSensitiveDetector(name), fHisto(0) {
  fHisto=G4TARCHistoManager::GetPointer();
}


G4bool G4TARCVirtualSD::ProcessHits(G4Step* myStep, G4TouchableHistory*){
  const G4Track* myTrack = myStep->GetTrack();
  fHisto->AddLeakingParticle(myTrack);
  fHisto->AddTargetStep(myStep);
  fHisto->ScoreNewTrack(myTrack);
  if (   (myTrack->GetParticleDefinition()->GetParticleName()=="neutron")
      && (myTrack->GetTrackID() > 1)
    )
    fHisto->AddEnergyTimeHole(myTrack, myStep);
  return true;
}
