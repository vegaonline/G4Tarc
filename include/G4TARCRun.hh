/*************************************************
* @file      G4TARCRun.hh
* @author    Abhijit Bhattacharyya
* @brief     This is for run i.e. to
*                 calculate energy deposition
************************************************/
#ifndef G4TARC_RUN_HH
#define G4TARC_RUN_HH 1

#include "G4Run.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4THitsMap.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4SystemOfUnits.hh"

#include <vector>
#include <regex>

class G4TARCRun : public G4Run {

public:
  G4TARCRun();
  virtual ~G4TARCRun(){};

public:
  ReadExperimentalDataFromFile(G4String& );

private:
  G4String                         fExptlDataFileName = "Data/TARC_EXPT_DATA/TARC_EXPTL_DATA.txt";
s
  std::vector< std::vector<G4double> >   fExptRadiiTables;
  std::vector< std::vector<G4double> >   fExptFluenceTables;
  std::vector< std::vector<G4double> >   fExptErrTables;
  std::vector< std::vector<G4double> >   fExptEnergyTables;
  std::vector< std::vector<G4double> >   fExptFluxTables;
  std::vector< std::vector<G4double> >   fExptFluxErrTables;
  std::vector< std::vector<G4double> >   fFlux_Radius;
  std::vector< std::vector<G4double> >    fRadialFluenceStep;


};
#endif
