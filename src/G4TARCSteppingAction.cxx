#include "G4TARCSteppingAction.hh"

G4TARCSteppingAction::G4TARCSteppingAction(G4TARCEventAction* anEvent)
: G4UserSteppingAction()//, fEventAction(anEvent), fHisto(0)
{
  fEventAction = anEvent;
  fHisto = G4TARCHistoManager::GetPointer();
}

void G4TARCSteppingAction::UserSteppingAction(const G4Step* myStep){
  fHisto->ProcessStepping(myStep);
}
