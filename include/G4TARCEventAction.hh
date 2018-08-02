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

  void analyseSecondaries(G4double, G4String, G4double, G4double, G4int, G4double, G4double, G4String, G4bool, G4int);
  void NeutronEnergyTime(G4double, G4double, G4double);
  void otherEnergyTime(G4double, G4double, G4double);
  void exitingTally(G4bool, G4double);
  void exitingTallyCheck(G4bool exiting_flag_check);
  void Adding2NeutronStack() {++fNeutronStack;}

  void analyseNeutronFlux(G4double, G4int, G4double, G4double,  G4String);
  void analyseNeutronShellFluence(G4double, G4double);
  void analyseNeutronRadialFluence(G4double, G4double, G4int); //G4double, G4int);
  void analysePS(G4double, G4String, G4double); // , G4double, G4double);
  void analyseNeutronFluence(G4double energyL, G4double timeL, G4int thisTrackIDL,
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


#endif
