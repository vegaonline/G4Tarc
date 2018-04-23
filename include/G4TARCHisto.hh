/***************************************************
 * @file           G4TARCHisto.hh
 * @author         Abhijit Bhattacharyya
 * @brief          This is for the histogram
 ***************************************************/

#ifndef G4TARC_HISTO_H
#define G4TARC_HISTO_H

#include "globals.hh"
#include <vector>
#include <string>
#include "G4TARCHistoMessenger.hh"
#include "G4TARCAnalysis.hh"
#include "G4TARCParallelWorld.hh"

#include "G4DataVector.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4SystemOfUnits.hh"


class G4TARCParallelWorld;
class G4TARCHistoMessenger;


class G4TARCHisto {
public:
  G4TARCHisto(){};
  ~G4TARCHisto();

public:
  //G4bool CreateTuples(G4AnalysisManager*);
  void SetGeomParam(G4int&, G4double&, std::vector<G4double>&, std::vector<G4double>&);
  void SetExptDataParam(G4int&, G4int&, G4int&, std::vector<G4double>&, std::vector<std::vector<G4double> >&,
      std::vector<std::vector<G4double> >&, std::vector<std::vector<G4double> >& );
  void InitializeParams();
  //void FillRadialExperimentalData(G4AnalysisManager*);
  void StartProcessing(G4AnalysisManager*);
  G4THitsMap<G4double>* GetHitsMap(G4int ndx) { return fRunMap[ndx]; } // get hits map by: (a) sequential Number
  G4THitsMap<G4double>* GetHitsMap(const G4String&, const G4String&);  // (b) multifun det name, (c) collection name
  G4THitsMap<G4double>* GetHitsMap(const G4String&);                   // (d) collection name with full path
  //virtual void Merge(const G4Run*);

  inline void SetVerbose(G4int val) { fVerbose = val; };


public:
  G4int                                 fVerbose;
  G4int                                 fShellNumber;
  G4double                              fShellThickness;
  G4double                              fMinInnerRadiusofShell;
  G4int                                 fTotal_flux;
  G4int                                 fN_max;
  G4int                                 fMaxRadIndex;  // number of shells consulted by fShellNumber
  G4int                                 fMaxDimEnInt      = 4; //  time and space (r) ?
  G4int                                 fMaxFluxData;     // = 100; // 21; was used in Alex's code for 21 input data
  G4int                                 fMaxlowFluxData;  // = 100;
  G4int                                 fMaxFluenceData;  // = 1000;
  G4int                                 fMaxTestFluxData;  // 21
  G4int                                 fMaxFluenceSpectrumData; // 1000
  std::vector<G4String>                 fCollName;
  std::vector<G4double>                 fInnerRadiusofShell;
  std::vector<G4int>                    fCollID;
  std::vector<G4double>                 fOuterRadiusofShell;
  std::vector<std::vector<G4double> >   fRadial_fluence_step;    // [26][10]
  std::vector<std::vector<G4double> >   fFlux_radius;
  std::vector<G4double>                 fExptEnergyBin;
  std::vector< std::vector<G4double> >  fExptRadiiTables;
  std::vector< std::vector<G4double> >  fExptFluenceTables;
  std::vector< std::vector<G4double> >  fExptErrTables;
  std::vector< std::vector<G4double> >  fExptEnergyTables;
  std::vector< std::vector<G4double> >  fExptFluxTables;
  std::vector< std::vector<G4double> >  fExptFluxErrTables;
  std::vector<G4THitsMap<G4double>* >   fRunMap;
  std::vector<G4int>                    fFluxTableList {36, 38, 40};
  std::vector<G4double> fRadiiBin {16.8*cm, 40.4*cm, 45.6*cm, 69.1*cm, 81.1*cm,
                           98.6*cm, 105.3*cm, 113.5*cm, 124.8*cm, 153.9*cm };

  G4int fGamma_flux, fNeutron_flux, fNeutron_check, fElectron_flux, fPositron_flux,
         fPiminus_flux, fPiplus_flux, fPizero_flux,  fProton_flux, fMuon_flux,
         fOther_flux, fExiting_flux, fExiting_check_flux, fNeutron_fluence, fNeutron_fluence_cyl;

  G4double                              fExiting_energy;
  G4int                                 fIntegral_scintillation;
  G4double                              fIntegral_scintillation_E;
  G4int                                 fIntegral_lithium;
  G4double                              fIntegral_lithium_E;
  G4int                                 fIntegral_helium;
  G4double                              fIntegral_helium_E;
  G4int                                 fDuplicate_neutron;
  G4int                                 fOldTrackID;
  G4double                              fFractional_bin_width;
  G4double                              fEnergy, fTime;
  G4double                              fName;
  G4double                              fNeutron_energy, fNeutron_time;
  G4double                              fIntegral_Eflux_p;
  G4double                              fIntegral_flux_p;
  G4double                              fIntegral_Eflux;
  G4double                              fLithium_integral;
  G4DataVector                          fFlux_energy;            // 22 ?
  G4DataVector                          fRadii;                  // 10
  G4DataVector                          fRadial_energies;        // 10
  G4DataVector                          fFlux_stat_error;        // 21
  G4DataVector                          fFlux_systematic_error;  // 21
  G4DataVector                          fFine_energy;
  G4DataVector                          fFlux;
  G4DataVector                          fFlux_data;
  G4DataVector                          fCos_flux;
  G4DataVector                          fFluence;
  G4DataVector    fFluence_spectrum;        // 1000
  G4DataVector    fFluence_step;
  G4DataVector    fFluence_front_step;
  G4DataVector    fFluence_cyl;
  G4DataVector    fFluence_step_cyl;
  G4DataVector    fFluence_step_shell;
  G4DataVector    fEFlux;
  G4DataVector    fFine_eFlux;
  G4DataVector    fEnergy_integral;        // 4
  G4DataVector    fEnFlux;                 // 4
  G4DataVector  fNeutFlux;               // 4
  G4DataVector  fLowEnergy;              // 101
  G4DataVector  fLowFluxData;            // 100
  G4DataVector  fLowFlux;
  G4DataVector  fCos_low_flux;
  G4DataVector  fLow_fluence;
  G4DataVector  fLow_fluence_step;
  G4DataVector  fLow_fluence_front_step;
  G4DataVector  fLow_fluence_cyl;
  G4DataVector  fLow_fluence_step_cyl;
  G4DataVector  fLow_fluence_step_shell;
  G4DataVector  fLow_stat;
  G4DataVector  fLow_syst;
  G4DataVector  fLithium_energy;
  G4DataVector  fLithium_radial_energy_lower;
  G4DataVector  fLithium_radial_energy_upper;
  G4DataVector  fLithium_radial_mean;
  G4DataVector  fLithium_radial_fluence_step;
  G4DataVector  fLithium_flux;
  G4DataVector  fLithium_flux_data;
  G4DataVector  fLithium_stat;
  G4DataVector  fCos_lithium_flux;
  G4DataVector  fLithium_fluence;
  G4DataVector  fLithium_fluence_Step;
  G4DataVector  fLithium_fluence_front_Step;
  G4DataVector  fLithium_fluence_cyl;
  G4DataVector  fLithium_Zflux;

};

#endif
