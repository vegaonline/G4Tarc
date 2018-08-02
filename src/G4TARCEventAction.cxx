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

void G4TARCEventAction::NeutronEnergyTime(G4double thisE, G4double thisT, G4double E0){
  G4cout << " Entering Neutron ET. " << G4endl;
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
  G4double tempT = thisT / microsecond;
  G4double tempE = thisE / eV;
  G4double tempE0 = E0 / eV;

  if (fAnalysisManager->IsActive()){
    if (tempT > 0.0 && tempE > 0.0) fAnalysisManager->FillH2(1, log10(tempT), log10(tempE), 1.0);
    fAnalysisManager->FillNtupleDColumn(1, 0, tempE);
    fAnalysisManager->FillNtupleDColumn(1, 1, tempT);
    fAnalysisManager->FillNtupleDColumn(1, 2, tempE0);
    fAnalysisManager->AddNtupleRow(1);
  }

  G4cout << " Exiting Neutron ET. " << G4endl;
}

void G4TARCEventAction::otherEnergyTime(G4double thisE, G4double thisT, G4double E0){
  G4cout << " Entering Other ET. " << G4endl;
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
  //G4cout << fNtuple_full << "  OTHERS-------> thisE : " << thisE << "  thisT : " << thisT << G4endl;
  G4double tempT = thisT / microsecond;
  G4double tempE = thisE / eV;
  G4double tempE0 = E0 / eV;

  if (fAnalysisManager->IsActive()){
  if (tempT > 0.0 && tempE > 0.0) fAnalysisManager->FillH2(2, log10(tempT), log10(tempE), 1.0);
    fAnalysisManager->FillNtupleDColumn(4, 0, tempE);
    fAnalysisManager->FillNtupleDColumn(4, 1, tempT);
    fAnalysisManager->FillNtupleDColumn(4, 2, tempE0);
    fAnalysisManager->AddNtupleRow(4);
    G4cout << " Exiting Neutron ET. " << G4endl;
  }
}


void G4TARCEventAction::exitingTally(G4bool exiting_flag, G4double energyL){
  G4cout << " Entering exiting_tally " << G4endl;
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
  G4cout << "exiting tally Analysis Manager activated" << G4endl;
  if(exiting_flag) {
    G4TARCRun* thisRun = static_cast<G4TARCRun*> (G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    G4cout << "exiting tally Run activated" << G4endl;
    thisRun->CalcExitingFlux(energyL);
    G4cout << "exiting tally calced" << G4endl;
    if (fAnalysisManager->IsActive()) {
    fAnalysisManager->FillNtupleDColumn(2, 0, energyL / eV);
    fAnalysisManager->AddNtupleRow(2);
  }
    G4cout << "Exiting from exiting_tally." << G4endl;
  }
}

void G4TARCEventAction::exitingTallyCheck(G4bool exiting_flag_check){
  G4cout << " Entering exiting_tally_check" << G4endl;
  if (exiting_flag_check){
    G4TARCRun* thisRun = static_cast<G4TARCRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    thisRun->exitingTallyCheck(exiting_flag_check);
  }
  G4cout << " Exiting exiting_tally_check" << G4endl;
}

void G4TARCEventAction::analysePS(G4double fParticleEnergy, G4String fParticleName, G4double fParticleMomentum
){
  // G4cout << " Inside analysePS" << G4endl;
  auto fAnalysisManager = G4AnalysisManager::Instance();
  if(fParticleName == "gamma") {
    // fAnalysisManager->FillH1(1, fParticleEnergy/eV);
  } else if(fParticleName == "neutron") {
    //++fNeutCap;
    //G4cout << fNeutCap << "     " << fNeutronStack << G4endl;
    //fAnalysisManager->FillH2(3, fNeutCap * 1e9, fParticleTime / microsecond, 1.0);
    //                 fAnalysisManager->FillH1(2, fParticleEnergy / eV, 1.0 / fParticleMomentum);
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
    //                fAnalysisManager->FillH1(3, fParticleEnergy / eV);
  } else if(fParticleName == "e+") {
    //         fAnalysisManager->FillH1(4, fParticleEnergy / eV);
  } else {   //(fParticleName == "other") {
    //                fAnalysisManager->FillH1(5, fParticleEnergy / eV);
  }
  //  G4cout << " Exiting  analysePS" << G4endl;
}



void G4TARCEventAction::analyseNeutronRadialFluence(G4double fParticleEnergyL, //G4double fParticleTimeL,
  G4double StepLengthL, G4int ishellL){
    G4cout << "Entering Neutron_Radial_Fluence" << G4endl;
    G4TARCRun* thisRun = static_cast<G4TARCRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    thisRun->analyseNeutronRadialFluence(fParticleEnergyL, StepLengthL, ishellL);
    G4cout << "Exiting Neutron_Radial_Fluence" << G4endl;
}


void G4TARCEventAction::analyseNeutronShellFluence(G4double energyL, G4double StepLengthL){
  G4cout << "Entering Neutron_Shell_Fluence" << G4endl;
  G4TARCRun* thisRun = static_cast<G4TARCRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  thisRun->analyseNeutronShellFluence(energyL, StepLengthL);
  G4cout << "Exiting Neutron_Shell_Fluence" << G4endl;
  }



  //void G4TARCHistoManager::analyseNeutronFluence(G4double energyL, G4String& nameL, G4double timeL, G4double momentumL,
  //  G4int thisTrackIDL, G4double radiusL, G4double thisStepL,  G4int ParentIDL, G4double primaryEnergyL,
  //  G4double parentEnergyL, G4String& parentParticleL, G4bool reduced_fluxL,  G4int number_generationsL){


  void G4TARCEventAction::analyseNeutronFluence(G4double energyL, G4double timeL, G4int thisTrackIDL,
    G4double radiusL, G4double thisStepL,  G4int ParentIDL, G4double parentEnergyL, G4String& parentParticleL){
      G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
      G4TARCRun* thisRun = static_cast<G4TARCRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
      thisRun->analyseNeutronFluence( energyL, thisStepL);   //,  timeL,  thisTrackIDL,  radiusL, thisStepL,  ParentIDL, parentEnergyL,parentParticleL);


      G4int iParent = 0;
      G4double fTempT    = timeL / microsecond;
      G4double fTempE0   = fEnergy0 / eV;
      G4double fTempE    = energyL / eV;
      if (parentParticleL == "gamma") iParent = 1;
      if (parentParticleL == "neutron") iParent = 2;
      if (parentParticleL == "e-") iParent = 3;
      if (parentParticleL == "pi-") iParent = 4;
      if (parentParticleL == "pi+") iParent = 5;
      if (parentParticleL == "pi0") iParent = 6;
      if (parentParticleL == "e+") iParent = 7;
      if (parentParticleL == "proton") iParent = 8;
/*
      //if (fTempE >= fFlux_Energy[fFlux_Energy.size()-1]){
      if (fTempE >= fFlux_Energy[0]){
        for (G4int ii = 0; ii < fMaxTestFluxData; ii++){
          if (fTempE > fFlux_Energy[ii] && fTempE < fFlux_Energy[ii + 1]) fFluence_step[ii] += thisStepL;
          // G4cout << fTempE << "   " << fFlux_Energy[0] << "   " << thisStepL / mm << "       " << fFluence_step[ii] << G4endl;
        }
      }
*/
      //fAnalysisManager->FillNtupleDColumn(7, 0, fTempE);
      fAnalysisManager->FillNtupleDColumn(3, 0, fTempE);
      fAnalysisManager->FillNtupleDColumn(3, 1, fTempT);
      fAnalysisManager->FillNtupleDColumn(3, 2, fTempE0);
      fAnalysisManager->FillNtupleIColumn(3, 3, thisTrackIDL);
      fAnalysisManager->FillNtupleIColumn(3, 4, ParentIDL);
      fAnalysisManager->FillNtupleDColumn(3, 5, 0.0);
      fAnalysisManager->FillNtupleDColumn(3, 6, 0.0);
      fAnalysisManager->FillNtupleDColumn(3, 7, 0.0);   // zMomentum
      fAnalysisManager->FillNtupleDColumn(3, 8, fTime0 / microsecond);
      fAnalysisManager->FillNtupleDColumn(3, 9, radiusL / mm);
      fAnalysisManager->FillNtupleDColumn(3, 10, parentEnergyL / eV);
      fAnalysisManager->FillNtupleIColumn(3, 11, iParent);
      fAnalysisManager->FillNtupleDColumn(3, 12, thisStepL);
      fAnalysisManager->FillNtupleIColumn(3, 13, 0);
      fAnalysisManager->AddNtupleRow(3);
  }


  void G4TARCEventAction::analyseNeutronFlux(G4double n_EnergyL, G4int thisTrackIDL, G4double radiusL, G4double cosAngleL, G4String fParticleNameL)
    //G4double zPosL,G4double cosAngleL, G4String fParticleNameL)
    {
      G4cout << "Entering Neutron_Flux" << G4endl;
      G4TARCRun* thisRun = static_cast<G4TARCRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
      thisRun->analyseNeutronFlux(n_EnergyL,  thisTrackIDL, radiusL, cosAngleL, fParticleNameL);
      G4cout << "Entering Neutron_Flux" << G4endl;
  }
