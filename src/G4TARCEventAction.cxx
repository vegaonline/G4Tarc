#include "G4TARCEventAction.hh"

G4TARCEventAction::G4TARCEventAction(): fNeutronStack(0)  {
    fDebugStarted = false;
    fEventMessenger = new G4TARCEventActionMessenger (this);
    fUITARC = G4UImanager::GetUIpointer();
    //fHisto = G4TARCHistoManager::GetPointer();
    fSelected = 0;
    SetPrintModulo(1);
}

G4TARCEventAction::~G4TARCEventAction() {
  delete fEventMessenger;
}

void G4TARCEventAction::BeginOfEventAction( const G4Event* evt ){
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


void G4TARCEventAction::EndOfEventAction( const G4Event* evt ) {
  G4TARCRun* thisRun = static_cast<G4TARCRun*> (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  //thisRun->FillPerEvent(999.);
  //G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
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



void G4TARCEventAction::NeutronEnergyTime(G4double thisE, G4double thisT, G4double E0){
  G4double tempT = thisT / microsecond;
  G4double tempE = thisE / eV;
  G4double tempE0 = E0 / eV;

  if (tempT > 0.0 && tempE > 0.0) fAnalysisManager->FillH2(1, log10(tempT), log10(tempE), 1.0);
  fAnalysisManager->FillNtupleDColumn(1, 0, tempE);
  fAnalysisManager->FillNtupleDColumn(1, 1, tempT);
  fAnalysisManager->FillNtupleDColumn(1, 2, tempE0);
  fAnalysisManager->AddNtupleRow(1);
}

void G4TARCEventAction::otherEnergyTime(G4double thisE, G4double thisT, G4double E0){
  //G4cout << fNtuple_full << "  OTHERS-------> thisE : " << thisE << "  thisT : " << thisT << G4endl;
  G4double tempT = thisT / microsecond;
  G4double tempE = thisE / eV;
  G4double tempE0 = E0 / eV;

  if (tempT > 0.0 && tempE > 0.0) fAnalysisManager->FillH2(2, log10(tempT), log10(tempE), 1.0);
    fAnalysisManager->FillNtupleDColumn(4, 0, tempE);
    fAnalysisManager->FillNtupleDColumn(4, 1, tempT);
    fAnalysisManager->FillNtupleDColumn(4, 2, tempE0);
    fAnalysisManager->AddNtupleRow(4);
}


void G4TARCEventAction::exitingTally(G4bool exiting_flag, G4double energyL){
  if(exiting_flag) {
    CalcExitingFlux(energyL);
    fAnalysisManager->FillNtupleDColumn(2, 0, energyL / eV);
    fAnalysisManager->AddNtupleRow(2);
  }
}

void G4TARCEventAction::analysePS(G4double fParticleEnergy, G4String fParticleName, G4double fParticleMomentum
){

  if(fParticleName == "gamma") {
    fAnalysisManager->FillH1(1, fParticleEnergy/eV);
  } else if(fParticleName == "neutron") {
    //++fNeutCap;
    //G4cout << fNeutCap << "     " << fNeutronStack << G4endl;
    //fAnalysisManager->FillH2(3, fNeutCap * 1e9, fParticleTime / microsecond, 1.0);
    fAnalysisManager->FillH1(2, fParticleEnergy / eV, 1.0 / fParticleMomentum);
    if(fParticleEnergy / MeV < 2.0) {
      fNeutflux[0] += 1.0;
      fENflux[0] += fParticleEnergy;
    } else if(fParticleEnergy / MeV >= 2.0 && fParticleEnergy / MeV < 20.0) {
      fNeutflux[1] += 1.0;
      fENflux[1] += fParticleEnergy / MeV;
    } else if(fParticleEnergy / MeV >= 20.0) {
      fNeutflux[2] += 1.0;
      fENflux[2] += fParticleEnergy / MeV;
    }
    if(fParticleEnergy / MeV >= 1000.0) {
      fNeutflux[3] += 1.0;
      fENflux[3] += fParticleEnergy / MeV;
    }
  } else if(fParticleName == "e-") {
    fAnalysisManager->FillH1(3, fParticleEnergy / eV);
  } else if(fParticleName == "e+") {
    fAnalysisManager->FillH1(4, fParticleEnergy / eV);
  } else {   //(fParticleName == "other") {
    fAnalysisManager->FillH1(5, fParticleEnergy / eV);
  }
}
