#include "G4TARCEventAction.hh"

G4TARCEventAction::G4TARCEventAction()
: G4UserEventAction()
//  ,fEventMessenger(0),
//  fUITARC(0),
//  fHisto(0),
//  fSelectedEvents(0),
//  fPrintModulo(1),
//  fSelected(0),
//  fDebugStarted(false)
  {
    fDebugStarted = false;
    fEventMessenger = new G4TARCEventActionMessenger (this);
    fUITARC = G4UImanager::GetUIpointer();
    fHisto = G4TARCHistoManager::GetPointer();
    fSelected = 0;
    SetPrintModulo(1);
}

G4TARCEventAction::~G4TARCEventAction() {
  delete fEventMessenger;
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
    G4cout << "EventAction: Event # " << nEvt << " started " ;
  }
  fHisto->BeginOfEvent(nEvt);
}


void G4TARCEventAction::EndOfEventAction( const G4Event* evt ) {
  if ( fDebugStarted ){
    fUITARC->ApplyCommand("/tracking/verbose 0");
    fDebugStarted = false;
    G4cout << " EventAction: Event # " << evt->GetEventID() << " ended" << G4endl;
  }
  fHisto->EndOfEvent();
  if (fHisto->GetVerbose() > 0) G4cout << "   EventAction: Event # " << evt->GetEventID() << " ended." <<G4endl;
}
