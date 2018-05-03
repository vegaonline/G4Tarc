#ifndef G4TARC_HISTO_MESSENGER_H
#define G4TARC_HISTO_MESSENGER_H

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"
#include <sstream>
#include "G4TARCHistoManager.hh"
//#include "G4TARCHisto.hh"

//class G4TARCHisto;
class G4TARCHistoManager;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithABool;


class G4TARCHistoMessenger : public G4UImessenger {
public:
  G4TARCHistoMessenger(G4TARCHisto*);
  virtual ~G4TARCHistoMessenger();
  virtual void SetNewValue(G4UIcommand*, G4String);

private:
  G4TARCHisto*                fHisto;
  G4TARCHistoManager*         fHistoM;
  G4UIdirectory*              fHistoDir;
  G4UIcmdWithAString*         fFactoryCmd;
  G4UIcmdWithAString*         fFileCmd;
  G4UIcommand*                fHistoCmd;
  G4UIcmdWithADoubleAndUnit*  fEdepCmd;
  G4UIcmdWithAnInteger*       fBinCmd;
  G4UIcmdWithAnInteger*       fVerbCmd;

};

#endif
