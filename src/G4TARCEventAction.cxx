#include "G4TARCEventAction.hh"

G4TARCEventAction::G4TARCEventAction()
: G4UserEventAction(),
  fEventMessenger(0),
  //fTARCHisto(0),
  fUITARC(0),
  fSelectedEvents(0),
  fPrintModulo(1),
  fSelected(0),
  fDebugStarted(false) {
    fEventMessenger = new G4TARCEventActionMessenger (this);
    //fTARCHisto = G4TARCHistoManager::GetPointer();
    fUITARC = G4UImanager::GetUIpointer();
}

G4TARCEventAction::~G4TARCEventAction() {
  delete fEventMessenger;
  //delete fTARCHisto;
}

void G4TARCEventAction::BeginOfEventAction( const G4Event* evt ){
  G4int nEvt = evt->GetEventID();

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
    G4cout << "EventAction: Event # " << nEvt << " started " << G4endl;
  }
  //fTARCHisto->SetVerbose(0);
  G4TARCHistoManager::GetPointer()->BeginOfEvent();
}


void G4TARCEventAction::EndOfEventAction( const G4Event* evt ) {
  if ( fDebugStarted ){
    fUITARC->ApplyCommand("/tracking/verbose 0");
    fDebugStarted = false;
    G4cout << " EventAction: Event # " << evt->GetEventID() << " ended" << G4endl;
  }

  G4TARCHistoManager* fTARCHisto = G4TARCHistoManager::GetPointer();
  fTARCHisto->EndOfEvent();
  if ( fTARCHisto->GetVerbose() > 1 ) {
    G4cout << " EventAction: Event " << evt->GetEventID() << " ended. " << G4endl;
  }
}
