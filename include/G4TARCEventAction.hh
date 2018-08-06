/*****************************************************************
 * @file    G4TARCEventAction.hh
 * @author  Abhijit Bhattacharyya
 * @brief   data members holds energy deposit and track length
 ****************************************************************/
#ifndef G4TARC_EVENTACTION_H
#define G4TARC_EVENTACTION_H

#include <vector>
#include <fstream>
#include <iomanip>

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4THitsMap.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4Run.hh"

#include "G4ios.hh"
#include "globals.hh"

#include "G4TARCRun.hh"
#include "G4TARCEventActionMessenger.hh"
#include "G4TARCAnalysis.hh"
//#include "G4TARCHistoManager.hh"


class G4Event;
class G4UImanager;
class G4TARCEventActionMessenger;
//class G4TARCHistoManager;

class G4TARCEventAction : public G4UserEventAction {
public:
  G4TARCEventAction();
  virtual ~G4TARCEventAction();

  virtual void BeginOfEventAction( const G4Event* );
  virtual void EndOfEventAction( const G4Event* );

  void SetPrintModulo ( G4int val) {fPrintModulo = val;}
  void AddEventToDebug ( G4int val){ G4cout << " start AddEvent2Debug" << G4endl;  fSelectedEvents.push_back(val); ++fSelected;  G4cout << " start AddEvent2Debug" << G4endl; }

  inline void analyseSecondaries(G4double, G4String, G4double, G4double, G4int, G4double, G4double, G4String, G4bool, G4int);
  inline void NeutronEnergyTime(G4double, G4double, G4double);
  inline void otherEnergyTime(G4double, G4double, G4double);
  inline void exitingTally(G4bool, G4double);
  inline void exitingTallyCheck(G4bool exiting_flag_check);
  void Adding2NeutronStack() {++fNeutronStack;}

  inline void analyseNeutronFlux(G4double, G4int, G4double, G4double,  G4String);
  inline void analyseNeutronShellFluence(G4double, G4double);
  inline void analyseNeutronRadialFluence(G4double, G4double, G4int); //G4double, G4int);
  inline void analysePS(G4double, G4String, G4double); // , G4double, G4double);
  inline void analyseNeutronFluence(G4double energyL, G4double timeL, G4int thisTrackIDL,
    G4double radiusL, G4double thisStepL,  G4int ParentIDL, G4double parentEnergyL, G4String& parentParticleL);


private:

  //    G4TARCEventAction& operator=(const G4TARCEventAction& right);
  //      G4TARCEventAction ( const G4TARCEventAction& );

  //G4TARCHistoManager*                        fHisto;
  //G4AnalysisManager*                         fAnalysisManager;

  G4TARCEventActionMessenger*                  fEventMessenger;
  G4UImanager*                                 fUITARC;
  std::vector<G4int>                           fSelectedEvents;
  G4bool                                                  fDebugStarted;
  G4int                                                     fPrintModulo;
  G4int                                                     fSelected;
  G4int                                                     fEventID;
  G4int                                                     fNeutronStack;
  G4int                                                      fGamma_flux;
  G4int                                                      fNeutron_flux;
  G4int                                                      fNeutron_check;
  G4int                                                      fElectron_flux;
  G4int                                                      fPiMinus_flux;
  G4int                                                      fPiPlus_flux;
  G4int                                                      fPiZero_flux;
  G4int                                                      fPositron_flux;
  G4int                                                      fProton_flux;
  G4int                                                      fMuon_flux;
  G4int                                                      fNeutron_fluence;
  G4int                                                      fOther_flux;
  G4int                                        fExiting_Flux;

  G4int                                  fMaxFluxData;
  G4int                                  fMaxFluenceData;
  G4int                                  fMaxTestFluxData;
  G4int                                  fIFluxCountRef;

  G4int                                  fMaxRadCount;
  G4int                                  fRefShellNumber;

  G4double                               fEnergy0, fTime0;
  G4double                               fExiting_Energy;
  std::vector<G4double>                  fEflux_Data;
  //std::vector<G4double>                fFine_Energy;
  std::vector<G4double>                  fFlux_Energy;
  std::vector<G4double>                  fFluence_step;

  std::vector<G4double>                  fENflux;
  std::vector<G4double>                  fNeutflux;
  std::vector<G4double>                  fLithium_Radial_Energy_Lower;
  std::vector<G4double>                  fLithium_Radial_Energy_Upper;
  std::vector<std::vector<G4double> >    fRadialFluenceStep;

};


inline void G4TARCEventAction::analyseSecondaries(G4double energyL, G4String nameL, G4double timeL, G4double momentumL,
  G4int ParentIDL, G4double primaryEnergyL, G4double parentEnergyL, G4String parentParticleL, G4bool reduced_fluxL,
  G4int number_generationsL){
    //G4cout << " Analyse Secondary started" << G4endl;

    G4AnalysisManager*  fAnalysisManager = G4AnalysisManager::Instance();
    G4TARCRun* thisRun = static_cast<G4TARCRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4int Iparticle=-9;
  G4double temp_time = timeL / microsecond;
  G4double temp_energy = energyL / eV;
  G4double temp_momentum = momentumL;

  if(nameL == "gamma") {
    thisRun->AddFlux(nameL);
    //fGamma_flux++;
    Iparticle = 1;
  } else if(nameL == "neutron") {
    if(!reduced_fluxL) thisRun->AddFlux("neutron_check"); // fNeutron_check++;
    //fNeutron_flux++;
    Iparticle = 2;
    if(reduced_fluxL) Iparticle = -2;
    if(energyL > 0.1*eV && energyL < 10.0*keV) thisRun->AddFlux("neutron_fluence"); // fNeutron_fluence++;
  } else if(nameL == "e-") {
    thisRun->AddFlux(nameL);
    //fElectron_flux++;
    Iparticle = 3;
    //    return;
  } else if(nameL == "pi-") {
    thisRun->AddFlux(nameL);
    //fPiMinus_flux++;
    Iparticle = 4;
  } else if(nameL == "pi+") {
    thisRun->AddFlux(nameL);
    // fPiPlus_flux++;
    Iparticle = 5;
  } else if(nameL == "pi0") {
    thisRun->AddFlux(nameL);
    // fPiZero_flux++;
    Iparticle = 6;
  } else if(nameL == "e+") {
    thisRun->AddFlux(nameL);
    // fPositron_flux++;
    Iparticle = 7;
  } else if(nameL == "proton") {
    thisRun->AddFlux(nameL);
    // fProton_flux++;
    Iparticle = 8;
  } else if(nameL == "mu-") {
    thisRun->AddFlux(nameL);
    // fMuon_flux++;
    Iparticle = 9;
  } else if(nameL == "mu+") {
    thisRun->AddFlux(nameL);
    // fMuon_flux++;
    Iparticle = 10;
  } else {
    thisRun->AddFlux("other");
    fOther_flux++;
    Iparticle = 99;
    return;
  }
  G4int iParent = 0;
  if (parentParticleL == "gamma")        iParent = 1;
  else if (parentParticleL == "neutron") iParent = 2;
  else if(reduced_fluxL)                 iParent = -2;
  else if (parentParticleL == "e-")      iParent = 3;
  else if (parentParticleL == "pi-")     iParent = 4;
  else if (parentParticleL == "pi+")     iParent = 5;
  else if (parentParticleL == "pi0")     iParent = 6;
  else if (parentParticleL == "e+")      iParent = 7;
  else if (parentParticleL == "proton")  iParent = 8;
  else if (parentParticleL == "proton")  iParent = 8;
  else if (parentParticleL == "mu-")     iParent = 9;
  else if (parentParticleL == "mu+")     iParent = 10;

  if (fAnalysisManager->IsActive()) {
    fAnalysisManager->FillNtupleDColumn(0,0, temp_energy);
    fAnalysisManager->FillNtupleDColumn(0,1, temp_time);
    fAnalysisManager->FillNtupleIColumn(0,2, Iparticle);
    fAnalysisManager->FillNtupleDColumn(0,3, temp_momentum);
    fAnalysisManager->FillNtupleIColumn(0,4, ParentIDL);
    fAnalysisManager->FillNtupleDColumn(0,5, primaryEnergyL);
    fAnalysisManager->FillNtupleIColumn(0,6, iParent);
    fAnalysisManager->FillNtupleDColumn(0,7, parentEnergyL);
    fAnalysisManager->FillNtupleIColumn(0,8, number_generationsL);
    fAnalysisManager->FillNtupleIColumn(0,9, fEventID);
    fAnalysisManager->AddNtupleRow(0);
  }

    //G4cout << " Analyse Secondary exiting" << G4endl;

}


inline void G4TARCEventAction::NeutronEnergyTime(G4double thisE, G4double thisT, G4double E0){
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

inline void G4TARCEventAction::otherEnergyTime(G4double thisE, G4double thisT, G4double E0){
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


inline void G4TARCEventAction::exitingTally(G4bool exiting_flag, G4double energyL){
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

inline void G4TARCEventAction::exitingTallyCheck(G4bool exiting_flag_check){
  G4cout << " Entering exiting_tally_check" << G4endl;
  if (exiting_flag_check){
    G4TARCRun* thisRun = static_cast<G4TARCRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    thisRun->exitingTallyCheck(exiting_flag_check);
  }
  G4cout << " Exiting exiting_tally_check" << G4endl;
}

inline void G4TARCEventAction::analysePS(G4double fParticleEnergy, G4String fParticleName, G4double fParticleMomentum
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



inline void G4TARCEventAction::analyseNeutronRadialFluence(G4double fParticleEnergyL, //G4double fParticleTimeL,
  G4double StepLengthL, G4int ishellL){
    G4cout << "Entering Neutron_Radial_Fluence" << G4endl;
    G4TARCRun* thisRun = static_cast<G4TARCRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    thisRun->analyseNeutronRadialFluence(fParticleEnergyL, StepLengthL, ishellL);
    G4cout << "Exiting Neutron_Radial_Fluence" << G4endl;
}


inline void G4TARCEventAction::analyseNeutronShellFluence(G4double energyL, G4double StepLengthL){
  G4cout << "Entering Neutron_Shell_Fluence" << G4endl;
  G4TARCRun* thisRun = static_cast<G4TARCRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  thisRun->analyseNeutronShellFluence(energyL, StepLengthL);
  G4cout << "Exiting Neutron_Shell_Fluence" << G4endl;
  }



  //void G4TARCHistoManager::analyseNeutronFluence(G4double energyL, G4String& nameL, G4double timeL, G4double momentumL,
  //  G4int thisTrackIDL, G4double radiusL, G4double thisStepL,  G4int ParentIDL, G4double primaryEnergyL,
  //  G4double parentEnergyL, G4String& parentParticleL, G4bool reduced_fluxL,  G4int number_generationsL){


  inline void G4TARCEventAction::analyseNeutronFluence(G4double energyL, G4double timeL, G4int thisTrackIDL,
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


  inline void G4TARCEventAction::analyseNeutronFlux(G4double n_EnergyL, G4int thisTrackIDL, G4double radiusL, G4double cosAngleL, G4String fParticleNameL)
    //G4double zPosL,G4double cosAngleL, G4String fParticleNameL)
    {
      G4cout << "Entering Neutron_Flux" << G4endl;
      G4TARCRun* thisRun = static_cast<G4TARCRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
      thisRun->analyseNeutronFlux(n_EnergyL,  thisTrackIDL, radiusL, cosAngleL, fParticleNameL);
      G4cout << "Entering Neutron_Flux" << G4endl;
  }


#endif
