#ifndef G4TARCRUN_HH
#define G4TARCRUN_HH


#include "G4Run.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4SystemOfUnits.hh"
#include "G4THitsMap.hh"
#include "globals.hh"
#include "G4TARCRun.hh"
#include <vector>
#include <regex>

class G4TARCRun : public G4Run {

public:
  G4TARCRun();
  virtual ~G4TARCRun();

public:
  void DefineShellBlocks();
  void ReadExperimentalDataFromFile(G4String&);

  void StartProcessing();
  virtual void RecordEvent(const G4Event*);
  G4int GetNumberOfHitsMap()  const {return fRunMap.size();}

  G4THitsMap<G4double>* GetHitsMap(G4int idx) { return fRunMap[idx];}
  G4THitsMap<G4double>* GetHitsMap(const G4String& detName, const G4String& colName);
  G4THitsMap<G4double>* GetHitsMap(const G4String& fullName);

  void DumpAllScorer();


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
  G4int GetRefShellNumber() const       {return fRefShellNumber;}

private:
  std::vector<G4String> fCollName;
  std::vector<G4int> fCollID;
  std::vector<G4THitsMap<G4double> *> fRunMap;

public:
  G4bool                      fStartHisto;
  G4bool                      fNtuple_full;
  G4bool                      fReadData;
  G4bool                      fInitialized;

  unsigned                                                 fMeanEnergyTable = 40;

  G4int                                                      fMaxRadCount;
  G4int                                                      fMaxTestFluxData;
  G4int                                                      fMaxFluxData;
  G4int                                                      fMaxFluenceData;

  G4int                                                      fTotal_flux;
  G4int                                                      fNmax;
  G4int                                                      fExiting_Flux;
  G4int                                                      fExiting_check_Flux;
  G4int                                                      fIntegral_flux_5cm;
  G4int                                                      fIntegral_flux_10cm;
  G4int                                                      fIntegral_flux_46cm;
  G4int                                                      fIntegral_flux_70cm;
  G4int                                                      fIntegral_flux_100cm;
  G4int                                                      fIntegral_flux_120cm;
  G4int                                                      fTARC_Integral_Eflux_46cm;
  G4int                                                      fDuplicate_neutrons;
  G4int                                                      fOldTrackID;
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
  G4double                                               fExiting_Energy;
  G4double                                               fFracBinWidth;

G4double floatDummy=0;

    G4double fTotVolVBox  = (150.0 * mm) * (150.0 * mm) * (300.0 * mm);  // volume of virtual box around holes
    G4double fMaxLVal     = 3000.0 * mm;
    G4double fMaxEVal     = 3000.0 * CLHEP::MeV;
    G4double fEVal0       = 3000.0 * CLHEP::MeV;
    //G4int    fNumMax      = 200;  // for fE/Msecond etc.
    G4int    fMaxBin      = 500;
    G4int    fNbin        = fMaxBin;
    // G4int    fMaxEBin     = 200;
    //G4int    fStepE       = (fMaxEVal / fMaxBin);
    G4int    fMaxSlices   = 3 * fMaxBin;
    //G4int    fNHisto      = 25;
    //G4int    fMaxNdx      = 500;
    G4int    fMaxFluenceTable = 0;
    G4double fMyTol       = 1.0e-9*mm;
    G4double fMyRadTol    = 1.0e-6*mm;


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

  G4double            fShellThickness;
  G4double            fRefShellThickness;
  G4double            fRefShellOuterRad;
  G4double         fRefShellInnerRad;
  G4double fRefShellVol;
  G4double fMinInnerRadiusofShell;
  G4double fMaxOuterRadiusofShell;
  G4double fInnerRadProtonShell;
  G4double fOuterRadProtonShell;
  G4double fRadHole;
  G4double fRadCyl;
  G4double fTestSphereRadius;
  G4double fTestSphereVolume;
  G4double fTestSphereSurfaceArea;
  G4double fTestShellVol;
  G4double fEnergy0;
  G4double                               fLenCyl;

 G4int                                  fIFluxCountRef;
 G4int                                  fRefShellNumber = fRadiusReference.size();
 G4int                                  fShellNumber;

 std::vector<G4double>                  fOuterRadiusofShell;
 std::vector<G4double>                  fInnerRadiusofShell;
 std::vector<G4double>                  fRadiusReference {200.0 * cm, 190.0 * cm, 185.0 * cm, 175.0 * cm, 165.0 * cm, 150.0 * cm,
   140.0 * cm, 130.0 * cm, 120.0 * cm, 110.0 * cm, 100.0 * cm, 90.0 * cm, 80.0 * cm, 70.0 * cm, 60.0 * cm, 50.0 * cm, 45.7 * cm,
   40.0 * cm, 30.0 * cm, 25.0 * cm, 20.0 * cm, 15.0 * cm, 10.0 * cm, 8.0 * cm, 5.0 * cm, 3.0 * cm};

  std::vector< std::vector<G4double> >   fExptRadiiTables;
  std::vector< std::vector<G4double> >   fExptFluenceTables;
  std::vector< std::vector<G4double> >   fExptErrTables;
  std::vector< std::vector<G4double> >   fExptEnergyTables;
  std::vector< std::vector<G4double> >   fExptFluxTables;
  std::vector< std::vector<G4double> >   fExptFluxErrTables;
  std::vector< std::vector<G4double> >   fFlux_Radius;
  std::vector<std::vector<G4double> >    fRadialFluenceStep;

  //std::vector< G4double>                 fExptEnergyBin;
  //  std::vector<G4double>                  fFluxRadTables;
  //std::vector<G4double>                  fRadList;
  //std::vector<G4double>                  fFlux;
  // std::vector<G4double>                  fFlux_Energy;
  //std::vector<G4double>                  fFlux_Data;
  //std::vector<G4double>                  fFlux_Syst_Err;
  std::vector<G4double>                  fFlux_Energy_in;
  std::vector<G4double>                  fFlux_Data_in;
  std::vector<G4double>                  fFlux_Syst_Err_in;

  //std::vector<G4double>                  fFlux_Low;
  //std::vector<G4double>                  fFlux_Low_Radius;
  std::vector<G4double>                  fFlux_Low_Energy;
  std::vector<G4double>                  fFlux_Low_Energy_in;
  std::vector<G4double>                  fFlux_Low_Data;
  std::vector<G4double>                  fFlux_Low_Data_in;
  std::vector<G4double>                  fFlux_Low_Syst_Err;
  std::vector<G4double>                  fFlux_Low_Syst_Err_in;

  //  std::vector<G4double>                  fFlux_Lithium;
  //std::vector<G4double>                  fFlux_Lithium_Radius;
  std::vector<G4double>                  fFlux_Lithium_Energy;
  std::vector<G4double>                  fFlux_Lithium_Energy_in;
  std::vector<G4double>                  fFlux_Lithium_Data;
  std::vector<G4double>                  fFlux_Lithium_Data_in;
  std::vector<G4double>                  fFlux_Lithium_Syst_Err;
  std::vector<G4double>                  fFlux_Lithium_Syst_Err_in;

  // std::vector<G4double>                  fFlux_He3;
  //std::vector<G4double>                  fFlux_He3_Energy;
  //std::vector<G4double>                  fFlux_He3_Energy_in;
  //std::vector<G4double>                  fFlux_He3_Data;
  //std::vector<G4double>                  fFlux_He3_Syst_Err;

  //std::vector<G4double>                  fFluence1D;
  //std::vector<G4double>                  fFluence_Radius;
  //std::vector<G4double>                  fFluence_Energy;
  //std::vector<G4double>                  fFluence_Data;
  //std::vector<G4double>                  fFluence_Syst_Err;

  //std::vector<G4double>                  fEflux_Data;
  //std::vector<G4double>                  fFine_Energy;

  //std::vector<G4double>                  fENflux;
  //std::vector<G4double>                  fNeutflux;

  // std::vector<G4double>                  fFluence_Spectrum;
  //std::vector<G4double>                  fLithium_Radial_Energy_Lower;
  //std::vector<G4double>                  fLithium_Radial_Energy_Upper;
  //std::vector<G4double>                  fLithium_Radial_Mean;
  //std::vector<G4double>                  fLithium_Radial_True_Mean;
  //std::vector<G4double>                  fLithium_Fluence_Step_Shell;
  //std::vector<G4double>                  fLithium_Fluence_Step;
  //std::vector<G4double>                  fLow_Fluence_Step_Shell;
  //std::vector<G4double>                  fFluence_Step_Shell;
  //std::vector<G4double>                  fFluence_step;
  std::vector<G4double>                  fLithium_Flux;
  //std::vector<G4double>                  fHe3_Flux;
  // std::vector<G4double>                  fCos_He3_Flux;
  std::vector<G4double>                  fCos_Lithium_Flux;
  //std::vector<G4double>                  fLow_Flux;
  std::vector<G4double>                  fCos_Low_Flux;
  //std::vector<G4double>                  fCos_Flux;
  //std::vector<G4double>                  fEFlux;
  //std::vector<G4double>                  fFluence_Cyl;
  //std::vector<G4double>                  fLow_Fluence_step;



};


#endif
