#include "G4TARCEventAction.hh"

G4TARCEventAction::G4TARCEventAction(): fNeutronStack(0)  {
    fDebugStarted = false;
    fEventMessenger = new G4TARCEventActionMessenger (this);
    fUITARC = G4UImanager::GetUIpointer();
    //fHisto = G4TARCHistoManager::GetPointer();
    auto fAnalysisManager = G4AnalysisManager::Instance();
    fSelected = 0;
    SetPrintModulo(1);
}

G4TARCEventAction::~G4TARCEventAction() {
  delete fEventMessenger;
}

void G4TARCEventAction::BeginOfEventAction( const G4Event* evt ){
  G4cout << "BOEventAction" << G4endl;
  fEventID = evt->GetEventID();
  fNeutronStack = 0;

  if (fEventID % fPrintModulo == 0) G4cout << " Begin of Event:  " << fEventID << G4endl;

/*
  if ( fSelected > 0 ){
    for (G4int i = 0; i < fSelected; ++i ) {
      if ( nEvt == fSelectedEvents[i] ){
        fUITARC->ApplyCommand("/random/saveThisEvent");
        fUITARC->ApplyCommand("/tracking/verbose 2");
        fDebugStarted = true;
        break;
      }
    }
  }
  if ( G4int( nEvt / fPrintModulo ) * fPrintModulo == nEvt ){
    G4cout << "EventAction: Event # " << nEvt << " started "  << G4endl;
  }
  fHisto->BeginOfEvent(nEvt);
*/
}


void G4TARCEventAction::EndOfEventAction( const G4Event* ) {
  G4cout << " Ending of event " << fEventID << G4endl;
  //G4TARCRun* thisRun = static_cast<G4TARCRun*> (G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  auto  fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->FillH1(6, fNeutronStack);
  if (fEventID % fPrintModulo == 0) G4cout << " End of event: " << fEventID << G4endl;
  /*
  G4int nEvt = evt->GetEventID();
  if ( fDebugStarted ){
    fUITARC->ApplyCommand("/tracking/verbose 0");
    fDebugStarted = false;
    // G4cout << " EventAction: Event # " << evt->GetEventID() << " ended" << G4endl;
  }
  fHisto->EndOfEvent();
  if (fHisto->GetVerbose() > 1) G4cout << "   EventAction: Event # " << nEvt << " ended." <<G4endl;
  */
}
