/******************************************************************************
 * @file      G4TARCCheckVolumeSD.hh
 * @author    Abhijit Bhattacharyya
 * @brief     This file is sensitive detector for the whole volume
 *****************************************************************************/
#ifndef G4TARCCHECKVOLUMESD_H
#define G4TARCCHECKVOLUMESD_H

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

#include "G4TARCHistoManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"

class G4Step;
class G4TouchableHistory;
class G4HCofThisEvent;
class G4TARCHistoManager;

class G4TARCCheckVolumeSD: public G4VSensitiveDetector {
public:
  G4TARCCheckVolumeSD(const G4String&);
  virtual ~G4TARCCheckVolumeSD(){};

  virtual void Initialize(G4HCofThisEvent*){};
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  virtual void EndofEvent(G4HCofThisEvent*){};
  virtual void clear() {};
  virtual void PrintAll() {};

private:
  G4TARCHistoManager* fHisto;
};

#endif
