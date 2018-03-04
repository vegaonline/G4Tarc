/****************************************************
 *  @file    G4TARCActionInitialization.hh
 *  @author  Abhijit Bhattacharyya
 *  @brief   Action Initialization
 ***************************************************/

#include "G4TARCActionInitialization.hh"

//G4TARCActionInitialization::G4TARCActionInitialization(G4TARCDetectorConstruction* det)
//: G4VUserActionInitialization(), fGeomConstruct(det){}
G4TARCActionInitialization::G4TARCActionInitialization(): G4VUserActionInitialization(){}

G4TARCActionInitialization::~G4TARCActionInitialization() {}


void G4TARCActionInitialization::BuildForMaster() const {
  SetUserAction(new G4TARCRunAction() );
}

void G4TARCActionInitialization::Build() const {
  SetUserAction(new G4TARCPrimaryGeneratorAction());
  SetUserAction(new G4TARCRunAction());
  SetUserAction(new G4TARCEventAction());
  SetUserAction(new G4TARCStackingAction());

}
