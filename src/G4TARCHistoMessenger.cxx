#include "G4TARCHistoMessenger.hh"

G4TARCHistoMessenger::G4TARCHistoMessenger(G4TARCHisto* histoG)
: G4UImessenger(), fHisto(histoG), fHistoDir(0),
  fFactoryCmd(0), fFileCmd(0), fHistoCmd(0) {

  fHistoDir = new G4UIdirectory("/tarc/histo/");
  fHistoDir->SetGuidance("Histograms control");

  fBinCmd = new G4UIcmdWithAnInteger("/tarc/histo/NumberOfBinsE", this);
  fBinCmd->SetGuidance("Set number of bins for energy");
  fBinCmd->SetParameterName("NEbins", false);
  fBinCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEdepCmd = new G4UIcmdWithADoubleAndUnit("/tarc/histo/MaxEdep", this);
  fEdepCmd->SetGuidance("set max energy in histogram");
  fEdepCmd->SetParameterName("edep", false);
  fEdepCmd->SetUnitCategory("Energy");
  fEdepCmd->SetRange("edep>0");
  fEdepCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fVerbCmd = new G4UIcmdWithAnInteger("/tarc/Verbose",this);
  fVerbCmd->SetGuidance("Set verbose for ");
  fVerbCmd->SetParameterName("verb",false);
  fVerbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFactoryCmd = new G4UIcmdWithAString("/tarc/histo/fileName",this);
  fFactoryCmd->SetGuidance("Set name for the Histograms file");

  fFileCmd = new G4UIcmdWithAString("/tarc/histo/setFileType", this);
  fFileCmd->SetGuidance("Set type (hbook, XML) for the histogram file");

  fHistoCmd = new G4UIcommand("/tarc/histo/setHisto", this);
  fHistoCmd->SetGuidance("Set binning of the Histo number ih: ");
  fHistoCmd->SetGuidance("  nBins; valMin; valMax; unit (of vmin and vmax)");

  G4UIparameter* ih = new G4UIparameter("ih", 'i', false);
  ih->SetGuidance("Histo number ; from 0 to MaxHisto - 1");
  fHistoCmd->SetParameter(ih);
  fHistoCmd->SetRange("ih>0");

  G4UIparameter* nbBins = new G4UIparameter("nbBins", 'i', false);
  nbBins->SetGuidance("number of bins");
  nbBins->SetParameterRange("nbBins > 0");
  fHistoCmd->SetParameter(nbBins);
  fHistoCmd->SetRange("nbBins>0");

  G4UIparameter* valMin = new G4UIparameter("valMin", 'd', false);
  valMin->SetGuidance("valMin, expressed in unit");
  fHistoCmd->SetParameter(valMin);

  G4UIparameter* valMax = new G4UIparameter("valMax", 'd', false);
  valMin->SetGuidance("valMax, expressed in unit");
  fHistoCmd->SetParameter(valMax);

  G4UIparameter* unit = new G4UIparameter("unit", 's', true);
  unit->SetGuidance("if ommitted valMin and valMax are assumed dimensionless");
  unit->SetDefaultValue("none");
  fHistoCmd->SetParameter(unit);
}


G4TARCHistoMessenger::~G4TARCHistoMessenger(){
  delete fFileCmd;
  delete fHistoCmd;
  delete fFactoryCmd;
  delete fHistoDir;
  delete fEdepCmd;
  delete fBinCmd;
  delete fVerbCmd;
}


void G4TARCHistoMessenger::SetNewValue(G4UIcommand* comm, G4String newVal) {
  //G4TARCHistoManager* histo = G4TARCHistoManager::GetPointer();
  if (comm == fFactoryCmd) {
    fHisto->SetFileName(newVal);
  } else if (comm == fFileCmd) {
    fHisto->SetFileType(newVal);
  } else if (comm == fBinCmd) {
    fHistoM->SetNumberOfBinsE(fBinCmd->GetNewIntValue(newVal));
  } else if (comm == fEdepCmd) {
    fHistoM->SetMaxEnergyDeposit(fEdepCmd->GetNewDoubleValue(newVal));
  } else if (comm == fVerbCmd) {
    fHistoM->SetVerbose(fVerbCmd->GetNewIntValue(newVal));
  } else if (comm == fHistoCmd) {
    G4int ih, nbBins;
    G4double vmin, vmax;

    std::istringstream iss(newVal);
    G4String fHistoUnits;
    iss >> ih >> nbBins >> vmin >> vmax >> fHistoUnits;

    G4String unit = fHistoUnits;
    G4double vUnit = 1.0 ;
    if(unit != "none") { vUnit = G4UIcommand::ValueOf(unit); }
    vUnit = (vUnit <= 0.0) ? 1.0 : vUnit;
    fHisto->SetHisto1D(ih,nbBins,vmin,vmax,vUnit);
  }
}
