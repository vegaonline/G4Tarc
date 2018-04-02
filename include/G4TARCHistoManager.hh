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


#include "G4TARCDetectorConstruction.hh"
#include "G4TARCHisto.hh"
#include "G4GeneralParticleSource.hh"

// #include "G4TARCAnalysis.hh"
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
#include "G4GeneralParticleSource.hh"

#include "g4root.hh"

//#include "TH1D.h"
//#include "TFile.h"
//#include "TTFree.h"

class G4TARCHisto;
class G4Track;
class G4Step;
class G4ParticleDefinition;
class G4GeneralParticleSource;
class G4TARCEventAction;

//class G4TARCDetectorConstruction;
// class G4AnalysisManager;

class G4TARCHistoManager {
public:
  static G4TARCHistoManager* GetPointer ();
  ~G4TARCHistoManager ();

private:
  G4TARCHistoManager ();

  void BookHisto();
  void save();
  void FillHisto(G4int, G4double, G4double);
  void Normalize(G4int, G4double);
  void FillNTuple(G4double, G4double, G4double, G4double);

public:
  void BeginOfRun ();
  void EndOfRun ();
  void BeginOfEvent ();
  void EndOfEvent ();

  void AddTargetStep(const G4Step*);
  void ScoreNewTrack(const G4Track*);
  void AddLeakingParticle(const G4Track*);
  void NeutFinalState(const G4Track*, const G4Step*);

  void TargetProfile(const G4Track*, const G4Step*);
  void AddEnergyTime(const G4Track*, const G4Step*);
  void AddEnergyTimeHole(const G4Track*, const G4Step*);
  void AddNzero( const G4Track*, const G4Step*);
  void GunParticleDistribution ( const G4Track*, const G4Step* );
  void WriteEventRange(G4ThreeVector, G4double, G4double);

  void TrackRun(G4double);
  void NeutronRun(G4double);
  void GunParticleRun(G4double);
  void Fill(G4int, G4double, G4double);

  inline G4int GetVerbose()               const           { return fVerbose; }
  inline G4double GetLength()             const           { return fLength; }
  inline G4double GetGPSEnergy()          const           { return fPrimaryKineticEnergy; }
  //inline G4double GetBeamEnergy()         const           { return G4GeneralParticleSource::GetParticleEnergy(); }

  inline void SetNumberOfBinsE(G4int val)                 { fNBinsE = val;}
  inline void SetMaxEnergyDeposit(G4double val)           { fEdepMax = val;}
  inline void SetVerbose ( G4int val)                     { fVerbose = val;}
  inline void SetGPSEnergyIN (const G4double value)       { fPrimaryKineticEnergy = value; }
  inline void SetGPSMomentumIN (const G4double value)     { fPrimaryMomentum = value; }
  inline void TotalProtonIn ()                            { fProtonIN++; }
  inline void TotalNCount()                               { fNCountTotal++; }

private:
  G4String                    fRootFileName;
  static G4TARCHistoManager*  fHistoManager;
  G4TARCHisto*                fHisto;
  //G4TARCDetectorConstruction* fDetector;

  const G4ParticleDefinition* fPrimaryDef;
  const G4ParticleDefinition* fNeutron;
  const G4ParticleDefinition* fProton;

  G4double fTotVolVBox  = (150.0 * mm) * (150.0 * mm) * (300.0 * mm);  // volume of virtual box around holes
  G4double fMaxLVal     = 5000.0 * mm;
  G4double fMaxEVal     = 8000.0 * CLHEP::MeV;
  G4double fEVal0       = 4000.0 * CLHEP::MeV;
  G4int    fNumMax      = 1000;  // for fE/Msecond etc.
  G4int    fMaxBin      = 100;
  G4int    fMaxEBin     = 10000;
  G4int    fStepE       = (fMaxEVal / fMaxBin);
  G4int    fMaxSlices   = 3 * fMaxBin;
  G4int    fNHisto      = 25;
  G4int    fMaxNdx      = 10000;


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
  G4int fNcpions;
  G4int fNpi0;
  G4int fNkaons;
  G4int fNstep;
  G4int fLMax;
  G4int fLBin;
  G4int fProtonIN;
  G4int fNCountTotal;

  G4bool fHistoBooked;

  G4double fEbin;
  G4double fNeutronInit,fNeutronSum, fNeutronBreed, fTimeMin, fTimeMax;
  std::vector<std::vector<G4double> > fET;
  std::vector<std::vector<G4double> > fNSpectra;
  std::vector<std::vector<G4double> > fEdNdE;
  std::vector<std::vector<G4double> > fFluence;
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
  G4DataVector       fnETsum;
  G4DataVector       fNEfluxBin;




  G4DataVector       fNSecondSum1;
  G4DataVector       fNSecondSum2;
  G4DataVector       fNSecondSum3;
  G4DataVector       fNSecondLow;

  G4ThreeVector      fRangeSum;
  G4double           fStepSum;
  G4double           fDeltaSum;

  G4PhysicsLogVector fnEsecond;
  G4PhysicsLogVector fnTsecond;

  size_t             fNbin;
};

#endif
