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
#include "globals.hh"
#include "G4Run.hh"
#include "G4UserRunAction.hh"
#include "G4TARCHistoManager.hh"
#include "G4TARCEventAction.hh"
// #include "G4TARCAnalysis.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <iomanip>

class G4Run;
class G4TARCEventAction;
//class G4TARCDetectorConstruction;
class G4TARCPrimaryGeneratorAction;
class G4UserRunAction;
class G4TARCHistoManager;

class G4TARCRunAction : public G4UserRunAction {
public:
  // G4TARCRunAction(G4TARCEventAction* eventAction);
  //G4TARCRunAction(G4TARCDetectorConstruction*, G4TARCPrimaryGeneratorAction*);
  G4TARCRunAction();
  virtual ~G4TARCRunAction();

public:
  //virtual G4Run* GenerateRun();
  virtual void BeginOfRunAction( const G4Run* ); // booking of histogram
  virtual void EndOfRunAction( const G4Run* );  // method fills histogram
  void CreateTuples();
  void BookHistogram();

  G4bool                        fHistoBooked;

private:
  G4String                      fAnalysisFileName = "G4TARC_output";
  G4AnalysisManager*            fAnalysisManager;
  G4TARCEventAction*            fEventAction;
  //G4TARCDetectorConstruction*   fDetector;
  G4TARCPrimaryGeneratorAction* fPrimary;
  G4TARCRunAction*              fRun;
  G4TARCHistoManager*           fHistoM;
  G4bool                        fInMaster;

};

#endif
