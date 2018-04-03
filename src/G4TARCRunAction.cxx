/*********************************************************
 * @file        G4TARCRunAction.cxx
 * @author      Abhijit Bhattacharyya
 * @brief       This is for the run Action
 ********************************************************/
#include "G4TARCRunAction.hh"

//G4TARCRunAction::G4TARCRunAction(G4TARCDetectorConstruction* det, G4TARCPrimaryGeneratorAction* prim)
//: G4UserRunAction(), fDetector(det), fPrimary(prim), fHisto(0){
G4TARCRunAction::G4TARCRunAction(): G4UserRunAction(){

}

G4TARCRunAction::~G4TARCRunAction() {

}

void G4TARCRunAction::BeginOfRunAction( const G4Run* aRun ) {
  auto id = aRun->GetRunID();
  G4cout << "Run # " << id << " starts." << G4endl;
  G4NuclearLevelData::GetInstance();
  (G4TARCHistoManager::GetPointer())->BeginOfRun();

#ifdef G4VIS_USE
  auto UI = G4UImanager::GetUIpointer();
  auto pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager) {
    UI->ApplyCommand("/vis/scene/notifyHandlers");    // this crashes the code .....
  }
#endif
}


void G4TARCRunAction::EndOfRunAction( const G4Run* aRun ){
  //auto analysisManager = G4AnalysisManager::Instance();
  G4cout << " RunAction: End of run actions for # " << aRun->GetRunID() << " is started" << G4endl;

  (G4TARCHistoManager::GetPointer())->EndOfRun();
  #ifdef G4VIS_USE
    if (G4VVisManager::GetConcreteInstance())
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  #endif
}
