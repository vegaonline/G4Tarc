/*******************************************************
 * @file      G4TARCPrimaryGeneratorAction.cxx
 * @author    Abhijit Bhattacharyya
 * @brief     Primary Generator using GPS
 ******************************************************/

#include "G4TARCPrimaryGeneratorAction.hh"

G4TARCPrimaryGeneratorAction::G4TARCPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(), fGPS(0), fEnergy(0.0*MeV)//,  fHisto(0)
 ,fGunMess(0), fBeamFlag(false), fCurrent(0)
{
    fGPS = new G4GeneralParticleSource();

    //fHisto = G4TARCHistoManager::GetPointer();
    fGunMess = new G4TARCPrimaryMessenger(this);
    //fGunMess = new G4TARCPrimaryMessenger();
}

G4TARCPrimaryGeneratorAction::~G4TARCPrimaryGeneratorAction() {
  delete fGPS;
}

void G4TARCPrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent ){
  if (DefaultBeamPosition()){
    G4cout << " Ok in Gun Mess\n";
  }
  fGPS->GeneratePrimaryVertex(anEvent);
  G4double tmp = fGPS->GetParticleEnergy();

  if (!GetEnergy()){
     SetEnergy(tmp);
     //(G4TARCHistoManager::GetPointer())->SetGPSEnergyIN(tmp);
   }

}
