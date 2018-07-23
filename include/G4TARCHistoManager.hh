/***************************************************
 * @file           G4TARCHistoManager.hh
 * @author         Abhijit Bhattacharyya
 * @brief          This is for the histogram
 ***************************************************/
#ifndef G4TARC_HISTOMANAGER_H
#define G4TARC_HISTOMANAGER_H

#include <cstdlib>
#include <string>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sstream>
#include <regex>
#include <fstream>


#include "G4TARCDetectorConstruction.hh"
// #include "G4TARCHisto.hh"
#include "G4GeneralParticleSource.hh"
#include "G4TARCAnalysis.hh"

#include "globals.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4UnitsTable.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZeroLong.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Version.hh"
#include "G4TrackStatus.hh"
#include "G4HadronicProcessType.hh"
#include "G4VProcess.hh"
#include "G4NuclearLevelData.hh"


class G4TARCHisto;
class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4TARCEventAction;

//class G4TARCDetectorConstruction;
// class G4AnalysisManager;

class G4TARCHistoManager {
public:
  static G4TARCHistoManager* GetPointer ();
  ~G4TARCHistoManager ();

private:
  G4TARCHistoManager ();

  void save();
  void FillHisto(G4int, G4double, G4double);
  void Normalize(G4int, G4double);
  void FillNTuple(G4double, G4double, G4double, G4double);

public:
  void BeginOfRun ();
  void EndOfRun ();
  void BeginOfEvent (G4int);
  void EndOfEvent ();

  void AddTargetStep(const G4Step*);
  void ScoreNewTrack(const G4Track*);
  void AddLeakingParticle(const G4Track*);
  void NeutFinalState(const G4Track*);   // , const G4Step*);

  void TargetProfile(const G4Track*);   // , const G4Step*);
  void AddEnergyTime(const G4Track*);   //, const G4Step*);
  void AddEnergyTimeHole(const G4Track*);    //      , const G4Step*);
  void AddNzero( const G4Track*);
  void GunParticleDistribution ( const G4Step* );
  void WriteEventRange(G4ThreeVector, G4double, G4double);
  G4double GetPrimaryEnergy() { return fEnergy0; }
  G4double GetPrimaryTime() { return fTime0; }

  void TrackRun(G4double);
  void NeutronRun(G4double);
  void GunParticleRun(G4double);
  //void Fill(G4int, G4double, G4double);

  void InitVectors();
  void ReadExperimentalDataFromFile(G4String&);
  void FillRadialExperimentalData();
  void BookHistogram();
  void CreateTuples();
  void NeutronFluxHistogram();
  void RadialFluxHistogram();
  void StartProcessing();
  void ProcessStepping(const G4Step*);
  void analysePS(G4double, G4String, G4double); // , G4double, G4double);
  void analyseNeutronFluence(G4double, G4double, G4int, G4double, G4double, G4int, G4double,  G4String&);

  void analyseNeutronRadialFluence(G4double, G4double, G4int); //G4double, G4int);
  void analyseNeutronShellFluence(G4double, G4double);
  void analyseNeutronFlux(G4double, G4int, G4double, G4double,  G4String);
  void analyseSecondaries(G4double, G4String, G4double, G4double, G4int, G4double, G4double, G4String, G4bool, G4int);
  void NeutronEnergyTime(G4double, G4double, G4double);
  void otherEnergyTime(G4double, G4double, G4double);
  void exitingTally(G4bool, G4double);

  void DefineShellBlocks();

  template <typename T>  void Check10s(T, T&, G4String&);

  inline void exitingTallyCheck(G4bool exiting_flag_check){if(exiting_flag_check) fExiting_check_Flux++; }
  inline void CalcExitingFlux(G4double exitingE)           { fExiting_Flux++;   fExiting_Energy += exitingE;}
  inline G4int GetVerbose()               const           { return fVerbose; }
  inline G4double GetLength()             const           { return fLength; }
  inline G4double GetGPSEnergy()          const           { return fPrimaryKineticEnergy; }
  //inline G4double GetBeamEnergy()         const           { return G4GeneralParticleSource::GetParticleEnergy(); }

  inline void SetEventID(G4int nVal)                      { fNevent_id = nVal;}
  inline void SetNumberOfBinsE(G4int val)                 { fNBinsE = val;}
  inline void SetMaxEnergyDeposit(G4double val)           { fEdepMax = val;}
  inline void SetVerbose ( G4int val)                     { fVerbose = val;}
  inline void SetGPSEnergyIN (const G4double value)       { fPrimaryKineticEnergy = value; }
  inline void SetGPSMomentumIN (const G4double value)     { fPrimaryMomentum = value; }
  inline void TotalProtonIn ()                            { fProtonIN++; }
  inline void TotalNCount()                               { fNCountTotal++; }
  inline void SetHistoBooked( G4bool val)                 { fHistoBooked = val;}

private:
  G4bool                      fStartHisto;
  G4bool                      fNtuple_full;
  G4bool                      fReadData;
  G4bool                      fInitialized;

  G4String                    fRootFileName;
  static G4TARCHistoManager*  fHistoManager;
  G4TARCHisto*                fHisto;

  const G4ParticleDefinition* fPrimaryDef;
  const G4ParticleDefinition* fNeutron;
  const G4ParticleDefinition* fProton;

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

  G4int fOldTrackID, fTrackID, fDuplicate_neutrons = 0;

  G4double fNstepEnergy;
  G4double fEdepMax;
  G4double fEdepEvt;
  G4double fEdepEM;
  G4double fEdepProton;
  G4double fEdepPI;
  G4double fEdepSum;
  G4double fEdepSum2;
  G4double fTcut;
  G4double fTkin;
  G4double fTmax;
  G4double fTmin;
  G4double fLength;
  G4double fHLength;
  G4double fRange;
  G4double fRho;
  G4double fAbsX0;
  G4double fAbsY0;
  G4double fAbsZ0;
  G4double fPrimaryKineticEnergy;
  G4double fPrimaryMomentum;
  G4double fVirtualDia;
  G4double fVirtVol;
  G4double fExiting_Energy;

  G4int fExiting_check_Flux;
  G4int fExiting_Flux;
  G4int fNevent_id;
  G4int fVerbose;
  G4int fNBinsE;
  G4int fNSlices;
  G4int fNinelastic;
  G4int fNelastic;
  G4int fNsecondary;
  G4int fNzero;
  G4int fNevt;
  G4int fNelec;
  G4int fNposit;
  G4int fNgam;
  G4int fNmuons;
  G4int fNions;
  G4int fNdeut;
  G4int fNalpha;
  G4int fNneutron;
  G4int fNproton;
  G4int fNprot_leak;
  G4int fNPionleak;
  G4int fNaproton;
  G4int fNneu_forw;
  G4int fNneu_leak;
  G4int fNneu_back;
  G4int fNeutronStack;
  G4int fNeutCap;
  G4int fNpions;
  G4int fNpi0;
  G4int fNkaons;
  G4int fNstep;
  G4int fLMax;
  G4int fLBin;
  G4int fProtonIN;
  G4int fNCountTotal;
  G4int fNeutron_check, fGamma_flux, fNeutron_flux, fElectron_flux, fPiminus_flux, fPiPlus_flux, fPizero_flux, fPositron_flux;
  G4int fProton_flux, fMuon_flux, fOther_flux, fNEUTRON_fluence;
  G4int fNumber_newTrack;

  G4bool fHistoBooked;

  G4double fEbin;
  G4double fNeutronInit,fNeutronSum, fNeutronBreed, fTimeMin, fTimeMax;
  std::vector<std::vector<G4double> > fET;
  std::vector<std::vector<G4double> > fNSpectra;
  // std::vector<std::vector<G4double> > fEdNdE;
  //std::vector<std::vector<G4double> > fFluence;
  std::vector<std::vector<G4double> > fETVirtual;


  G4DataVector       fGunParticleX;
  G4DataVector       fGunParticleY;
  G4DataVector       fGunParticleZ;
  G4DataVector       fGunParticleTLX;
  G4DataVector       fGunParticleTLY;
  G4DataVector       fGunParticleTLZ;
  G4DataVector       fGunParticleDep;
  G4DataVector       fGunParticleRho;
  G4DataVector       fRangeVector;
  G4DataVector       fRhoVector;
  G4DataVector       fDeltaVector;
  G4DataVector       fEsecond;
  G4DataVector       fMsecond;
  G4DataVector       fNETsum;
  //G4DataVector       fNEfluxBin;

  G4DataVector       fNSecondSum1;
  G4DataVector       fNSecondSum2;
  G4DataVector       fNSecondSum3;
  G4DataVector       fNSecondLow;
  G4DataVector                           fLocal_Energy_Integral;

  G4ThreeVector      fRangeSum;
  G4double           fStepSum;
  G4double           fDeltaSum;

  G4PhysicsLogVector fNEsecond;
  G4PhysicsLogVector fNTsecond;
//--------------------------------------------------------
  G4AnalysisManager*                    fAnalysisManager;
  G4String                                        fExptlDataFileName = "Data/TARC_EXPT_DATA/TARC_EXPTL_DATA.txt";
  G4String                                        fAnalysisFileName = "G4TARC_output";
  G4TARCEventAction*                 fEventAction;

  G4double                               fTestSphereRadius;
  G4double                               fTestSphereVolume;
  G4double                               fTestSphereSurfaceArea;
  G4double                               fTestShellVol;
  G4double                               fHalfXBlockB;
  G4double                               fHalfYBlockB;
  G4double                               fHalfZBlockB;
  G4double                               fHalfXVBox;
  G4double                               fHalfYVBox;
  G4double                               fHalfZVBox;
  G4double                               fNewHalfZProt;
  G4double                               fZposProt;
  G4int                                  fShellNumber;
  G4double                               fShellThickness;
  G4int                                  fRefShellNumber;
  G4double                               fRefShellThickness;
  G4double                               fRefShellOuterRad;
  G4double                               fRefShellInnerRad;
  G4double                               fRefShellVol;
  G4double                               fMaxOuterRadiusofShell;
  G4double                               fMinInnerRadiusofShell;
  G4double                               fInnerRadProtonShell;
  G4double                               fOuterRadProtonShell;
  G4double                               fRadHole;
  G4double                               fRadCyl;
  G4double                               fLenCyl;

  G4int                                  fMaxFluxData;
  G4int                                  fMaxFluenceData;
  G4int                                  fMaxTestFluxData;
  G4int                                  fIFluxCountRef;
  G4int                                  fMaxRadCount;
  G4int                                  fIntegral_flux_5cm, fIntegral_flux_10cm, fIntegral_flux_46cm,
                                  fIntegral_flux_70cm, fIntegral_flux_100cm, fIntegral_flux_120cm;
                                  //fLithium_flux_5cm;

  G4double                               fTARC_Integral, fTARC_Integral_E, fTARC_lithium, fTARC_lithium_IntegralData, fTARC_lithium_E;
  G4double                               fTARC_helium, fTARC_helium_E, fEflux_Integral, fTARC_Integral_Eflux_46cm;
  G4double                               fTotal_flux;
  G4double fAbs_Int_Scint_E;

  G4bool flag;
  std::map<G4int, G4double, std::less<G4int> > fParentEnergy;
  std::map<G4int, G4String, std::less<G4int> > fParentParticle;
  std::map<G4int, G4int, std::less<G4int> > fParentParticleID;
  G4int number_generations, fNmax;
  G4double fEnergy0, fTime0, fFracBinWidth;

  std::vector<G4double>                  fOuterRadiusofShell;
  std::vector<G4double>                  fInnerRadiusofShell;

  std::vector<G4double>                  fRadiusReference {200.0 * cm, 190.0 * cm, 185.0 * cm, 175.0 * cm, 165.0 * cm, 150.0 * cm,
    140.0 * cm, 130.0 * cm, 120.0 * cm, 110.0 * cm, 100.0 * cm, 90.0 * cm, 80.0 * cm, 70.0 * cm, 60.0 * cm, 50.0 * cm, 45.7 * cm,
    40.0 * cm, 30.0 * cm, 25.0 * cm, 20.0 * cm, 15.0 * cm, 10.0 * cm, 8.0 * cm, 5.0 * cm, 3.0 * cm};

  G4double                                               fAbsolute_TotalFlux, fAbsolute_Flux;

  G4double floatDummy=0.0;
  G4int      intDummy = 0;

  unsigned                                                 fMeanEnergyTable = 40;
  std::vector<G4double>                          fMeanEnergyT40List;
  std::vector<G4int>                                 fFluxTableList {36, 38};     // , 40}; the energy supplied is E_low

  std::vector< std::vector<G4double> >   fExptRadiiTables;
  std::vector< std::vector<G4double> >   fExptFluenceTables;
  std::vector< std::vector<G4double> >   fExptErrTables;
  std::vector< std::vector<G4double> >   fExptEnergyTables;
  std::vector< std::vector<G4double> >   fExptFluxTables;
  std::vector< std::vector<G4double> >   fExptFluxErrTables;
  std::vector< std::vector<G4double> >   fFlux_Radius;
  std::vector<std::vector<G4double> >    fRadialFluenceStep;


  std::vector< G4double>                 fExptEnergyBin;
  //  std::vector<G4double>                  fFluxRadTables;
  std::vector<G4double>                  fRadList;
  std::vector<G4double>                  fFlux;
  std::vector<G4double>                  fFlux_Energy;
  std::vector<G4double>                  fFlux_Data;
  std::vector<G4double>                  fFlux_Syst_Err;
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

  std::vector<G4double>                  fEflux_Data;
  //std::vector<G4double>                  fFine_Energy;

  std::vector<G4double>                  fENflux;
  std::vector<G4double>                  fNeutflux;

  //std::vector<G4double>                  fFluence_Spectrum;
  std::vector<G4double>                  fLithium_Radial_Energy_Lower;
  std::vector<G4double>                  fLithium_Radial_Energy_Upper;
  std::vector<G4double>                  fLithium_Radial_Mean;
  std::vector<G4double>                  fLithium_Radial_True_Mean;
  std::vector<G4double>                  fLithium_Fluence_Step_Shell;
  std::vector<G4double>                  fLithium_Fluence_Step;
  std::vector<G4double>                  fLow_Fluence_Step_Shell;
  std::vector<G4double>                  fFluence_Step_Shell;
  std::vector<G4double>                  fFluence_step;
  std::vector<G4double>                  fLithium_Flux;
  //std::vector<G4double>                  fHe3_Flux;
  // std::vector<G4double>                  fCos_He3_Flux;
  std::vector<G4double>                  fCos_Lithium_Flux;
  std::vector<G4double>                  fLow_Flux;
  std::vector<G4double>                  fCos_Low_Flux;
  std::vector<G4double>                  fCos_Flux;
  std::vector<G4double>                  fEFlux;
  //std::vector<G4double>                  fFluence_Cyl;
  std::vector<G4double>                  fLow_Fluence_step;

//  G4double testMax1 = -99999999.99e8, testMax2 = testMax1, testMin1 = -testMax1, testMin2 = testMin1;    // for histo check

};

#endif
