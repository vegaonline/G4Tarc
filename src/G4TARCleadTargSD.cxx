
#include "G4TARCleadTargSD.hh"

G4TARCleadTargSD::G4TARCleadTargSD(G4String& name)
: G4VSensitiveDetector(name), fHisto(0) {
  fHisto = G4TARCHistoManager::GetPointer();
}

G4bool G4TARCleadTargSD::ProcessHits(G4Step* myStep, G4TouchableHistory*){
  G4double edep = myStep->GetTotalEnergyDeposit();
  if (edep == 0.0) return false;
  const G4Track* myTrack = myStep->GetTrack();
  G4TouchableHandle touch1 = myStep->GetPreStepPoint()->GetTouchableHandle();
  G4int copyNum = touch1->GetReplicaNumber(1);

  if ( (myTrack->GetTrackID()==1 ) && (myTrack->GetDefinition()->GetParticleName()== "proton") )
    fHisto->TotalProtonIn();
  return true;
}
