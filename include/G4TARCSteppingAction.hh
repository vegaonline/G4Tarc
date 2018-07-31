#ifndef G4TARCSteppingAction_HH
#define G4TARCSteppingAction_HH

#include "G4UserSteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "globals.hh"
#include <map>

#include "G4TARCParallelWorld.hh"
#include "G4TARCEventAction.hh"

class G4TARCEventAction;
class G4TARCParallelWorld;

class G4TARCSteppingAction : public G4UserSteppingAction {
public:
  G4TARCSteppingAction(G4TARCEventAction*);
  virtual ~G4TARCSteppingAction(){};

  virtual void UserSteppingAction(const G4Step*);
  void ProcessStepping(const G4Step*);


private:
  //G4TARCHistoManager*  fHisto;
  G4TARCEventAction*   fEventAction;
  G4bool flag;

  G4double            fRefShellOuterRad;
  G4double         fRefShellInnerRad;
  G4double     fMaxOuterRadiusofShell;
  G4double fMinInnerRadiusofShell;
  G4double fEnergy0, fTime0;
  G4double fMyTol       = 1.0e-9*mm;
  G4double fMyRadTol    = 1.0e-6*mm;
  G4double fRefShellThickness;
  G4double fShellThickness;

  G4int number_generations, fNmax;
   G4int                                  fIFluxCountRef;
   G4int                                  fRefShellNumber = fRadiusReference.size();
   G4int                                  fShellNumber;
   G4int                                                      fMaxRadCount;

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


  std::map<G4int, G4double, std::less<G4int> > fParentEnergy;
  std::map<G4int, G4String, std::less<G4int> > fParentParticle;
  std::map<G4int, G4int, std::less<G4int> > fParentParticleID;
};

#endif
