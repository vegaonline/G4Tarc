/*************************************************
 * @file      G4TARCRunAction.hh
 * @author    Abhijit Bhattacharyya
 * @brief     This is for run Action i.e. to
 *                 calculate energy deposition
 ************************************************/
#ifndef G4TARC_RUNACTION_H
#define G4TARC_RUNACTION_H

#include "G4UserRunAction.hh"
#include "G4NuclearLevelData.hh"
#include "G4Run.hh"
#include "G4ios.hh"
#include "G4UserRunAction.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4SystemOfUnits.hh"

#include "globals.hh"

#include "G4TARCRun.hh"
#include "G4TARCEventAction.hh"
#include "G4TARCAnalysis.hh"

#include "Randomize.hh"
#include <fstream>
#include <iomanip>

class G4Run;
class G4TARCRun;
class G4TARCEventAction;
//class G4TARCDetectorConstruction;
class G4TARCPrimaryGeneratorAction;
class G4UserRunAction;


class G4TARCRunAction : public G4UserRunAction {
public:
  // G4TARCRunAction(G4TARCEventAction* eventAction);
  //G4TARCRunAction(G4TARCDetectorConstruction*, G4TARCPrimaryGeneratorAction*);
  G4TARCRunAction();
  virtual ~G4TARCRunAction();

public:
  virtual G4Run* GenerateRun();
  virtual void BeginOfRunAction( const G4Run* ); // booking of histogram
  virtual void EndOfRunAction( const G4Run* );  // method fills histogram

  void DefineShellBlocks();
  void ReadExperimentalDataFromFile(G4String&);
  void BookHistogram();
  void CreateTuples();
  void FillRadialExperimentalData();
  void NeutronFluxHistogram(G4int, const G4TARCRun*);
  void RadialFluxHistogram(G4int, const G4TARCRun*);


private:
  G4bool    fHistoBooked;
  G4bool   fReadData;

  G4String                                             fExptlDataFileName = "Data/TARC_EXPT_DATA/TARC_EXPTL_DATA.txt";
  G4String                                             fAnalysisFileName = "G4TARC_output";
  G4AnalysisManager*                         fAnalysisManager;
  G4TARCEventAction*                      fEventAction;
  //G4TARCDetectorConstruction*      fDetector;
  G4TARCPrimaryGeneratorAction*  fPrimary;
  G4TARCRunAction*                        fRun;
  //G4TARCHistoManager*                   fHistoM;

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

  G4int                                  fMaxFluxData;
  G4int                                  fMaxFluenceData;
  G4int                                  fMaxTestFluxData;
  G4int                                  fIFluxCountRef;
  G4int                                  fMaxRadCount;
  G4int                                  fRefShellNumber = fRadiusReference.size();
  G4int                                  fShellNumber;

  std::vector<G4double>                  fOuterRadiusofShell;
  std::vector<G4double>                  fInnerRadiusofShell;
  std::vector<G4double>                  fRadiusReference {200.0 * cm, 190.0 * cm, 185.0 * cm, 175.0 * cm, 165.0 * cm, 150.0 * cm,
    140.0 * cm, 130.0 * cm, 120.0 * cm, 110.0 * cm, 100.0 * cm, 90.0 * cm, 80.0 * cm, 70.0 * cm, 60.0 * cm, 50.0 * cm, 45.7 * cm,
    40.0 * cm, 30.0 * cm, 25.0 * cm, 20.0 * cm, 15.0 * cm, 10.0 * cm, 8.0 * cm, 5.0 * cm, 3.0 * cm};

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

};

#endif
