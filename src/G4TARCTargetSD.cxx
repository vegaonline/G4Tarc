#include "G4TARCTargetSD.hh"

G4TARCTargetSD::G4TARCTargetSD(const G4String& name)
: G4VSensitiveDetector(name), fHisto(0) {
  fHisto = G4TARCHistoManager::GetPointer();
}
G4bool G4TARCTargetSD::ProcessHits(G4Step* myStep, G4TouchableHistory*){
  const G4Track* myTrack = myStep->GetTrack();

  G4StepPoint* point1 = myStep->GetPreStepPoint();
  G4TouchableHandle touch1 = point1->GetTouchableHandle();
  G4StepPoint* point2 = myStep->GetPostStepPoint();
  G4TouchableHandle touch2 = point2->GetTouchableHandle();

  G4cout << "PreStep "  << touch1->GetVolume()->GetName()
         << " ****** "
         << "PostStep " << touch2->GetVolume()->GetName() << G4endl;


  fHisto->AddTargetStep(myStep);
  fHisto->ScoreNewTrack(myTrack);
  if (myTrack->GetTrackID() > 1) {
    fHisto->AddEnergyTime(myTrack, myStep);
  } else {
    fHisto->GunParticleDistribution(myTrack, myStep);
  }
  return true;
}
