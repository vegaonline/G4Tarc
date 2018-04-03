/*****************************************************************
 * @file    G4TARCEventAction.hh
 * @author  Abhijit Bhattacharyya
 * @brief   data members holds energy deposit and track length
 ****************************************************************/
#ifndef G4TARC_EVENTACTION_H
#define G4TARC_EVENTACTION_H

#include <vector>

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "G4TARCEventActionMessenger.hh"
#include "G4TARCHistoManager.hh"


class G4Event;
class G4UImanager;
class G4TARCEventActionMessenger;
class G4TARCHistoManager;

class G4TARCEventAction : public G4UserEventAction {
public:
  G4TARCEventAction();
  virtual ~G4TARCEventAction();

  virtual void BeginOfEventAction( const G4Event* );
  virtual void EndOfEventAction( const G4Event* );

  inline void SetPrintModulo ( G4int val) {fPrintModulo = val;}
  inline void AddEventToDebug ( G4int val){fSelectedEvents.push_back(val); ++fSelected;}

private:
  G4TARCEventAction& operator=(const G4TARCEventAction& right);
  G4TARCEventAction ( const G4TARCEventAction& );
  G4TARCHistoManager*                fHisto;
  G4TARCEventActionMessenger*        fEventMessenger;
  G4UImanager*                       fUITARC;
  std::vector<G4int>                 fSelectedEvents;
  G4bool                             fDebugStarted;
  G4int                              fPrintModulo;
  G4int                              fSelected;
};

#endif
