#ifndef G4TARCRUN_HH
#define G4TARCRUN_HH


#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4SystemOfUnits.hh"
#include "G4THitsMap.hh"
#include "globals.hh"
#include "G4TARCRun.hh"
#include <vector>

G4TARCRun::G4TARCRun():G4Run() {

public:
  G4TARCRun();
  virtual ~G4TARCRun();

public:
  void DefineShellBlocks();
  void ReadExperimentalDataFromFile(G4String&);

  void StartProcessing();
  virtual void RecordEvent(const G4Event*);
  G4int GetNumberOfHitsMap()  const (return fRunMap.size());
  /*
  G4THitsMap<G4double>* GetHitsMap(G4int idx) { return fRunMap[idx]};
  G4THitsMap<G4double>* GetHitsMap(const G4String& detName, const G4String& colName);
  G4THitsMap<G4double>* GetHitsMap(const G4String& fullName);

  void DumpAllScorer();
  */

  virtual void Merge(const G4Run*);

  inline void CalcExitingFlux(G4double exitingE)                   { fExiting_Flux++;   fExiting_Energy += exitingE;}
  inline void exitingTallyCheck(G4bool exiting_flag_check)   {if(exiting_flag_check) fExiting_check_Flux++; }

  void AddFlux(G4String);

  void analyseNeutronFlux(G4double, G4int, G4double, G4double,  G4String);
  void analyseNeutronShellFluence(G4double, G4double);
  void analyseNeutronRadialFluence(G4double, G4double, G4int); //G4double, G4int);

  G4int GetExitingFlux() const              {return fExiting_Flux;}
  G4int GetExitingEnergy() const          {return fExiting_Energy;}
  G4int GetExitingCheckFlux() const    {return fExiting_check_Flux;}
  G4int GetIntegralFlux_46cm() const   {return fIntegral_flux_46cm;}
  G4int GetIntegralEFlux_46cm() const {return fTARC_Integral_Eflux_46cm;}

private:
  std::vector<G4String fCollName;
  std::vector<G4int> fCollID;
  std::vector<G4THitsMap<G4double> *> fRunMap;

public:
  unsigned                                                 fMeanEnergyTable = 40;

  G4int                                                      fMaxRadCount;
  G4int                                                      fMaxTestFluxData;
  G4int                                                      fMaxFluxData;
  G4int                                                      fMaxFluenceData;

  G4int                                                      fTotal_flux;
  G4int                                                      fNMax;
  G4int                                                      fExiting_Flux;
  G4int                                                      fExiting_check_Flux;
  G4int                                                      fIntegral_flux_5cm;
  G4int                                                      fIntegral_flux_10cm;
  G4int                                                      fIntegral_flux_46cm;
  G4int                                                      fIntegral_flux_70cm;
  G4int                                                      fIntegral_flux_100cm;
  G4int                                                      fIntegral_flux_120cm;
  G4int                                                      fTARC_Integral_Eflux_46cm
  G4int                                                      fDuplicate_neutrons;
  G4int                                                      fOldTrackID;
  G4int                                                      fGamma_flux;
  G4int                                                      fNeutron_flux;
  G4int                                                      fNeutron_check;
  G4int                                                      fElectron_flux;
  G4int                                                      fPiminus_flux;
  G4int                                                      fPiplus_flux;
  G4int                                                      fPizero_flux;
  G4int                                                      fPositron_flux;
  G4int                                                      fProton_flux;
  G4int                                                      fMuon_flux;
  G4int                                                      fOther_flux;
  G4int                                                      fNeutron_fluence;
  G4int                                                      fNeutron_fluence_46cm;

  G4double                                               fTARC_Integral;
  G4double                                               fTARC_Integral_E;
  G4double                                               fTARC_lithium;
  G4double                                               fTARC_lithium_IntegralData;
  G4double                                               fTARC_lithium_E;
  G4double                                               fTARC_helium;
  G4double                                               fTARC_helium_E;
  G4double                                               fEflux_Integral;
  G4double                                               fTARC_Integral_Eflux_46cm;
  G4double                                               fExiting_Energy;
  G4double                                               fFracBinWidth;

  std::vector<G4int>                                 fFluxTableList {36, 38};     // , 40}; the energy supplied is E_low

  std::vector<G4double>                          fMeanEnergyT40List;
  std::vector<G4double>                          fFlux_Energy;
  std::vector<G4double>                          fFluence_Spectrum;
  std::vector<G4double>                          fRadList;
  std::vector< G4double>                         fExptEnergyBin;
  std::vector<G4double>                          fFlux_Data;
  std::vector<G4double>                          fEflux_Data;
  std::vector<G4double>                          fFlux_Syst_Err;
  std::vector<G4double>                          fFlux;
  std::vector<G4double>                          fCos_Flux;
  std::vector<G4double>                          fFluence_Data;
  std::vector<G4double>                          fFluence_step;
  std::vector<G4double>                          fFluence_Cyl;
  std::vector<G4double>                          fFluence_Step_Shell;
  std::vector<G4double>                          fFluence_Step_cyll;
  std::vector<G4double>                          fEFlux;
  std::vector<G4double>                          fLocal_Energy_Integral;
  std::vector<G4double>                          fENflux;
  std::vector<G4double>                          fNeutflux;
  std::vector<G4double>                          fLow_Flux;
  std::vector<G4double>                          fLow_Flux_Energy;
  std::vector<G4double>                          fLow_Flux_Data;
  std::vector<G4double>                          fLow_Flux_Syst_Err;
  std::vector<G4double>                          fLow_Cos_Flux;
  std::vector<G4double>                          fLow_Fluence;
  std::vector<G4double>                          fLow_Fluence_step;
  std::vector<G4double>                          fLow_Fluence_cyl;
  std::vector<G4double>                          fLow_Fluence_Step_Shell;
  std::vector<G4double>                          fLithium_Radial_Energy_Lower;
  std::vector<G4double>                          fLithium_Radial_Energy_Upper;
  std::vector<G4double>                          fLithium_Radial_Mean;
  std::vector<G4double>                          fLithium_Radial_True_Mean;
  std::vector<G4double>                          fLithium_Fluence_Step_Shell;
  std::vector<G4double>                          fLithium_Fluence_Step;
  std::vector<G4double>                          fLithium_Flux_Energy;
  std::vector<G4double>                          fLithium_Flux_Data;
  std::vector<G4double>                          fLithium_Flux_Syst_Err;

  std::vector< std::vector<G4double> >   fFlux_Radius;

  std::vector< std::vector<G4double> >   fExptRadiiTables;
  std::vector< std::vector<G4double> >   fExptFluenceTables;
  std::vector< std::vector<G4double> >   fExptErrTables;
  std::vector< std::vector<G4double> >   fExptEnergyTables;
  std::vector< std::vector<G4double> >   fExptFluxTables;
  std::vector< std::vector<G4double> >   fExptFluxErrTables;
  std::vector< std::vector<G4double> >   fFlux_Radius;
  std::vector<std::vector<G4double> >    fRadialFluenceStep;


};


#endif
