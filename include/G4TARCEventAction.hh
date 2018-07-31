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


//class G4Event;
//class G4UImanager;
class G4TARCEventActionMessenger;
//class G4TARCHistoManager;

class G4TARCEventAction : public G4UserEventAction {
public:
  G4TARCEventAction();
  virtual ~G4TARCEventAction();

  virtual void BeginOfEventAction( const G4Event* );
  virtual void EndOfEventAction( const G4Event* );

  inline void SetPrintModulo ( G4int val) {fPrintModulo = val;}
  inline void AddEventToDebug ( G4int val){fSelectedEvents.push_back(val); ++fSelected;}

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
  G4TARCEventAction& operator=(const G4TARCEventAction& right);
  G4TARCEventAction ( const G4TARCEventAction& );
  //G4TARCHistoManager*                      fHisto;
  G4AnalysisManager*                           fAnalysisManager;
  G4TARCEventActionMessenger*        fEventMessenger;
  G4UImanager*                                      fUITARC;
  std::vector<G4int>                                fSelectedEvents;
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

  G4int                                                      fMaxRadCount;
  G4int                      fRefShellNumber;

  G4double fEnergy0, fTime0;
  G4double                         fExiting_Energy;
  std::vector<G4double>                  fEflux_Data;
  //std::vector<G4double>                  fFine_Energy;
  std::vector<G4double>                  fFlux_Energy;
  std::vector<G4double>                          fFluence_step;

  std::vector<G4double>                  fENflux;
  std::vector<G4double>                  fNeutflux;
  std::vector<G4double>                          fLithium_Radial_Energy_Lower;
  std::vector<G4double>                          fLithium_Radial_Energy_Upper;
  std::vector<std::vector<G4double> >    fRadialFluenceStep;

};

#endif
