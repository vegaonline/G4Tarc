/***************************************************
 * @file           G4TARCHisto.hh
 * @author         Abhijit Bhattacharyya
 * @brief          This is for the histogram
 ***************************************************/

#ifndef G4TARC_HISTO_H
#define G4TARC_HISTO_H

#include "globals.hh"
#include "G4DataVector.hh"
#include <vector>
#include <string>
#include "G4TARCHistoMessenger.hh"
#include "G4RootAnalysisManager.hh"

class G4RootAnalysisManager;
class G4TARCHistoMessenger;


class G4TARCHisto {
public:
  G4TARCHisto();
  ~G4TARCHisto();

  void Book();   // Book predefined Histograms
  void Save();   // Save histograms to file

  void Add1D(const G4String&, const G4String&, G4int, G4double, G4double, G4double);  // 1D Histograms
  void SetHisto1D(G4int, G4int, G4double, G4double, G4double);                        // change bin/boundary
  void Activate(G4int, G4bool);                                                      //(de-)activate
  void Fill(G4int, G4double, G4double);                                               // Filled Histo
  void ScaleH1(G4int, G4double);                                                      // Histo scaled
  void AddTuple(const G4String&);                                                     // nTuple is booked
  void AddTupleI(const G4String&);                                                    // nTuple is booked
  void AddTupleF(const G4String&);
  void AddTupleD(const G4String&);
  void FillTupleI(G4int, G4int);                                                      // Fill nTuple parameter
  void FillTupleF(G4int, G4float);
  void FillTupleD(G4int, G4double);
  void AddRow();                                                                     // Save tuple event
  void SetFileName(const G4String&);                                                 // Set output file
  void SetFileType(const G4String&);
  void SetFullFileNameWithPath(const G4String inputStr)   { fFullHistoFileName = inputStr;}

  inline void SetVerbose(G4int val) { fVerbose = val; };
  inline G4bool IsActive() const { return fHistoActive; };

private:
  G4RootAnalysisManager*        fRMan;
  G4TARCHistoMessenger*         fMess;

  G4String                      fHistName;
  G4String                      fFullHistoFileName;
  G4String                      fHistType;
  G4String                      fTupleName;
  G4String                      fTupleTitle;
  G4int                         fNHisto;
  G4int                         fVerbose;
  G4bool                        fDefaultAct;
  G4bool                        fHistoActive;
  G4bool                        fNtupleActive;

  std::vector<G4int>            fHisto;
  std::vector<G4int>            fTupleI;
  std::vector<G4int>            fTupleF;
  std::vector<G4int>            fTupleD;
  std::vector<G4int>            fBins;
  std::vector<G4bool>           fActive;
  std::vector<G4double>         fXmin;
  std::vector<G4double>         fXmax;
  std::vector<G4double>         fUnit;
  std::vector<G4String>         fIds;
  std::vector<G4String>         fTitles;

  std::vector<G4String>         fNtupleI;
  std::vector<G4String>         fNtupleF;
  std::vector<G4String>         fNtupleD;
};

#endif
