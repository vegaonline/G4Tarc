#include "G4TARCHisto.hh"

G4TARCHisto::G4TARCHisto()
: fRMan(0), fMess(0), fHistName("test"), fHistType("root"),
  fTupleName("tuple"), fTupleTitle("test"), fNHisto(0), fVerbose(0),
  fDefaultAct(true), fHistoActive(false), fNtupleActive(false) {
    fMess = new G4TARCHistoMessenger(this);
}

G4TARCHisto::~G4TARCHisto() {
  delete fMess;
  delete fRMan;
}

void G4TARCHisto::Book() {
  if (!(fHistoActive || fNtupleActive)) return;

  fRMan = G4RootAnalysisManager::Instance();  // Create

  G4String fPath = std::string(std::getenv("dateStr"));
  G4String nameHistoFile = fPath + "/" + fHistName + "." + fHistType;   // creating a new hbook file
  SetFullFileNameWithPath(nameHistoFile);

  if (!fRMan->OpenFile(nameHistoFile)){       // open histogram file
      G4cout << "Histo::Book Error to open file < " << nameHistoFile <<  " >" << G4endl;
      fHistoActive = false;
      fNtupleActive = false;
      return;
  }
  G4cout << "### Histo::Save: Opened file < " << nameHistoFile << " > for " << fNHisto << " histograms " << G4endl;

  for (G4int i = 0; i < fNHisto; ++i) {   // creating 1D Histogram in root dir of the tree
    if (fActive[i]){
      G4String ss = "h" + fIds[i];
      fHisto[i] = fRMan->CreateH1(ss, fTitles[i], fBins[i], fXmin[i], fXmax[i]);
      if (fVerbose > 0)
        G4cout << "Created histogram #" << i << "  id= " << fHisto[i]
               << "  "  << ss << "  " << fTitles[i] << G4endl;
    }
  }

  if (fNtupleActive){   // creting tuple factory: tuples handled by tree
      fRMan->CreateNtuple(fTupleName, fTupleTitle);
      G4int i;
      G4int n;

      n = fNtupleI.size();
      for (i = 0; i < n; ++i)
        if (fTupleI[i] == -1)
          fTupleI[i] = fRMan->CreateNtupleIColumn(fNtupleI[i]);

      n = fNtupleF.size();
      for (i = 0; i < n; ++i)
        if (fTupleF[i] == -1)
          fTupleF[i] = fRMan->CreateNtupleIColumn(fNtupleF[i]);

      n = fNtupleD.size();
      for (i = 0; i < n; ++i)
        if (fTupleD[i] == -1)
          fTupleD[i] = fRMan->CreateNtupleIColumn(fNtupleD[i]);
  }

}

void G4TARCHisto::Save() {
  if (!(fHistoActive || fNtupleActive)) return;

  G4String nameHistoFile = fFullHistoFileName; //fHistName + "." + fHistType;  // Create a tree mapped to a new HBOOK file
  if (!fRMan->Write()){      // write histo file
    G4cout << "Histo::Save: FATAL ERROR writing ROOT file" << G4endl;
    exit(1);
  }
  if(fVerbose > 0)
    G4cout << "### Histo::Save: Histograms and Ntuples are saved" << G4endl;
  if(fRMan->CloseFile() && fVerbose > 0)
    G4cout << "                 File is closed" << G4endl;
  delete G4RootAnalysisManager::Instance();
  fRMan = 0;
}

void G4TARCHisto::Add1D(const G4String& id, const G4String& name,
  G4int nb, G4double x1, G4double x2, G4double u){
  if(fVerbose > 0)
    G4cout << "Histo::Add1D: New histogram will be booked: #"
           << id << "  <" << name
           << "  " << nb << "  " << x1 << "  " << x2 << "  " << u
           << G4endl;
  ++fNHisto;
  x1 /= u;
  x2 /= u;
  fActive.push_back(fDefaultAct);
  fBins.push_back(nb);
  fXmin.push_back(x1);
  fXmax.push_back(x2);
  fUnit.push_back(u);
  fIds.push_back(id);
  fTitles.push_back(name);
  fHisto.push_back(-1);
}

void G4TARCHisto::SetHisto1D(G4int i, G4int nb, G4double x1, G4double x2, G4double u){
  if(i>=0 && i<fNHisto){
    if(fVerbose > 0)
      G4cout << "Histo::SetHisto1D: #" << i
             << "  " << nb << "  " << x1 << "  " << x2 << "  " << u
             << G4endl;
    fBins[i] = nb;
    fXmin[i] = x1;
    fXmax[i] = x2;
    fUnit[i] = u;
    fActive[i] = true;
    fHistoActive = true;
  } else {
    G4cout << "Histo::SetHisto1D: WARNING! wrong histogram index "
           << i << G4endl;
  }
}

void G4TARCHisto::Activate(G4int i, G4bool val){
  if(fVerbose > 1) {
    G4cout << "Histo::Activate: Histogram: #" << i << "   "
           << val << G4endl;
  }
  if(i>=0 && i<fNHisto) {
    fActive[i] = val;
    if(val) { fHistoActive = true; }
  }
}

void G4TARCHisto::Fill(G4int i, G4double x, G4double w){
  if(!fHistoActive) return;
  if(fVerbose > 1) {
    G4cout << "Histo::Fill: Histogram: #" << i << " at x= " << x
           << "  weight= " << w
           << G4endl;
  }
  if(i>=0 && i<fNHisto) {
    if(fActive[i]) { fRMan->FillH1(fHisto[i], x/fUnit[i], w); }
  } else {
    G4cout << "Histo::Fill: WARNING! wrong histogram index " << i << G4endl;
  }
}

void G4TARCHisto::ScaleH1(G4int i, G4double x){
  if(!fHistoActive) { return; }
  if(fVerbose > 0) {
    G4cout << "Histo::Scale: Histogram: #" << i
           << " by factor " << x << G4endl;
  }
  if(i>=0 && i<fNHisto) {
    if(fActive[i]) { fRMan->GetH1(fHisto[i])->scale(x); }
  } else {
    G4cout << "Histo::Scale: WARNING! wrong histogram index " << i << G4endl;
  }
}


void G4TARCHisto::AddTuple(const G4String& w1){
  fTupleTitle = w1;
}

void G4TARCHisto::AddTupleI(const G4String& w1){
  fNtupleActive = true;
  fNtupleI.push_back(w1);
  fTupleI.push_back(-1);
}


void G4TARCHisto::AddTupleF(const G4String& w1){
  fNtupleActive = true;
  fNtupleF.push_back(w1);
  fTupleF.push_back(-1);
}


void G4TARCHisto::AddTupleD(const G4String& w1){
  fNtupleActive = true;
  fNtupleD.push_back(w1);
  fTupleD.push_back(-1);
}

void G4TARCHisto::FillTupleI(G4int i, G4int x){
  if(!fNtupleActive) { return; }
  G4int n = fNtupleI.size();
  if(i >= 0 && i < n) {
    if(fVerbose > 1) {
      G4cout << "Histo::FillTupleI: i= " << i << "  id= " << fTupleI[i]
             << "   <" << fNtupleI[i] << "> = " << x << G4endl;
    }
    fRMan->FillNtupleIColumn(fTupleI[i], x);
  } else {
    G4cout << "Histo::FillTupleI: WARNING! wrong ntuple index "
           << i << G4endl;
  }
}

void G4TARCHisto::FillTupleF(G4int i, G4float x){
  if(!fNtupleActive) { return; }
  G4int n = fNtupleF.size();
  if(i >= 0 && i < n) {
    if(fVerbose > 1) {
      G4cout << "Histo::FillTupleF: i= " << i << "  id= " << fTupleF[i]
             << "   <" << fNtupleF[i] << "> = " << x << G4endl;
    }
    fRMan->FillNtupleFColumn(fTupleF[i], x);
  } else {
    G4cout << "Histo::FillTupleF: WARNING! wrong ntuple index "
           << i << G4endl;
  }
}

void G4TARCHisto::FillTupleD(G4int i, G4double x){
  if(!fNtupleActive) { return; }
  G4int n = fNtupleD.size();
  if(i >= 0 && i < n) {
    if(fVerbose > 1) {
      G4cout << "Histo::FillTupleD: i= " << i << "  id= " << fTupleD[i]
             << "   <" << fNtupleD[i] << "> = " << x << G4endl;
    }
    fRMan->FillNtupleDColumn(fTupleD[i], x);
  } else {
    G4cout << "Histo::FillTupleD: WARNING! wrong ntuple index "
           << i << G4endl;
  }
}

void G4TARCHisto::AddRow(){
  if(!fNtupleActive) { return; }
  fRMan->AddNtupleRow();
}


void G4TARCHisto::SetFileName(const G4String& nam){
  fHistName = nam;
  fHistoActive = true;
}

void G4TARCHisto::SetFileType(const G4String& nam){
  if(nam == "root" || nam == "ROOT" )   { fHistType = "root"; }
  else if(nam == "xml" || nam == "XML") { fHistType = "xml"; }
  else if(nam == "ascii" || nam == "ASCII" ||
          nam == "Csv" || nam == "csv" || nam == "CSV")
    { fHistType = "ascii"; }
}
