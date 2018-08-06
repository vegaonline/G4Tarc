/*********************************************************
 * @file        G4TARCRunAction.cxx
 * @author      Abhijit Bhattacharyya
 * @brief       This is for the run Action
 ********************************************************/
#include "G4TARCRunAction.hh"

G4TARCRunAction::G4TARCRunAction() : G4UserRunAction() {
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();

  DefineShellBlocks();

  if (!fHistoBooked){
    BookHistogram();
    CreateTuples();
  }
}

G4TARCRunAction::~G4TARCRunAction() {

}

G4Run* G4TARCRunAction::GenerateRun() {
  return new G4TARCRun;
}

void G4TARCRunAction::BeginOfRunAction(const G4Run* thisRun){
  G4cout << " Run # " << thisRun->GetRunID() << " starts. " << G4endl;
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();

  fAnalysisManager->OpenFile(fAnalysisFileName);
}

void G4TARCRunAction::EndOfRunAction(const G4Run* thisRun) {
  //   ReadExperimentalDataFromFile(fExptlDataFileName);   use in G4TARCRun and bring as tarcRun->param[ijk] in FillRadial
  FillRadialExperimentalData();

  G4cout << "Number of events: " << thisRun->GetNumberOfEvent() << G4endl;

  auto fAnalysisManager = G4AbalysisManager::Instance();
  fAnalysisManager->Write();
  fAnalysisManager->CloseFile();
}


void G4TARCRunAction::DefineShellBlocks() {
  for (G4int ii = 0; ii < fRefShellNumber; ii++) {
    G4double radThis = fRadiusReference[ii] / mm;
    fInnerRadiusofShell.push_back(radThis - fRefShellThickness);
    fOuterRadiusofShell.push_back(radThis);
  }
  fRefShellThickness = 2.0 * mm;
}


void G4TARCRunAction::BookHistogram() {
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
  fHistoBooked = true;
  fAnalysisManager->SetFirstHistoId(1);
  //1
  fAnalysisManager->CreateH1("Gamma","Gamma Edep (eV)", 2e3, 1.0e5, 5e8);  // 0:
  //2
  fAnalysisManager->CreateH1("NeutronEnergy","Neutron energy (eV) vs. 1/mom /eV", 5e3, 1.0e5, 1.0e9); // , 0.0, 50.0);  // 100000, 0., 1000000.);
  //3
  fAnalysisManager->CreateH1("ElectronEdep","Electron Edep (eV)", 2e3, 1.0e4, 1.0e7);
  //4
  fAnalysisManager->CreateH1("PositronEdep","Positron Edep (keV)", 2e2, 10.0, 6.0e7);
  //5
  fAnalysisManager->CreateH1("OtherEdep","Other Edep (eV)", 1e3, 1.0e6, 5.0e9);
  //6
  fAnalysisManager->CreateH1("ParticleStack","Particle Stack", 50, 0.0, 1e9);
  //7
  fAnalysisManager->CreateH1("NeutronPerEvent","Neutrons/event", 30, 0.0, 1e5);
  //8
  fAnalysisManager->CreateH1("ProtonPerEvent","Protons/event", 50, 0.0, 30.0);
  //9
  fAnalysisManager->CreateH2("NeutronET","log(Neutron Energy <eV>) vs. log(Time <us> )", 50, -2, 4, 50, -3, 7); // H2:1
  //10
  fAnalysisManager->CreateH2("OtherPartET","log(OTHER particle Energy <eV>) vs. log(Time <us>)", 50, -3.0, 4.0, 50, -4.0, 6.0); // H2:2
  //11
  //fAnalysisManager->CreateH2("NeutronCapture", "Neutron Capture / 10^9 p", 700, 0, 10000, 400, 0, 1e9);  //H2:3
  G4cout << "BookHisto in RunAction done." << G4endl;
}


void G4TARCRunAction::CreateTuples(){
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->CreateNtuple("h1_Secondary", "Secondary Particle Info");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("time");
  fAnalysisManager->CreateNtupleIColumn("particle");
  fAnalysisManager->CreateNtupleDColumn("momentum");
  fAnalysisManager->CreateNtupleIColumn("parentid");
  fAnalysisManager->CreateNtupleDColumn("e_prim");
  fAnalysisManager->CreateNtupleIColumn("parent");
  fAnalysisManager->CreateNtupleDColumn("e_parent");
  fAnalysisManager->CreateNtupleIColumn("numgen");
  fAnalysisManager->CreateNtupleIColumn("event");
  fAnalysisManager->FinishNtuple(); // ntupleID: 0 - filled     : 0

  fAnalysisManager->CreateNtuple("h2_N_ET", "Neutron Time");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("time");
  fAnalysisManager->CreateNtupleDColumn("primary");
  fAnalysisManager->FinishNtuple(); // ntupleID: 1   : 1

  fAnalysisManager->CreateNtuple("h3_N_Exiting", "Neutrons Exiting");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->FinishNtuple(); // ntupleID: 2 - filled     : 2

  fAnalysisManager->CreateNtuple("h4_Flux_4002", "Neutrons G4TARC flux");
  fAnalysisManager->CreateNtupleDColumn("energy");     // 0
  fAnalysisManager->CreateNtupleDColumn("tarcflux");   // 1
  fAnalysisManager->CreateNtupleDColumn("errstat");    // 2
  fAnalysisManager->CreateNtupleDColumn("g4flux");     // 3
  fAnalysisManager->CreateNtupleDColumn("g4perp");     // 4
  fAnalysisManager->CreateNtupleDColumn("g4fluence");  // 5
  fAnalysisManager->CreateNtupleDColumn("g4err");      // 6
  fAnalysisManager->CreateNtupleDColumn("rawflux");    // 7
  fAnalysisManager->CreateNtupleDColumn("eflux");   // 8
  fAnalysisManager->CreateNtupleDColumn("g4eflux");    // 9
  fAnalysisManager->CreateNtupleDColumn("gstep");      // 10
  fAnalysisManager->CreateNtupleDColumn("gfl_cyl");    // 11
  fAnalysisManager->CreateNtupleDColumn("g4front");    // 12
  fAnalysisManager->CreateNtupleDColumn("g4shell");   // 13
  fAnalysisManager->CreateNtupleDColumn("g4shellerr");   // 14
  fAnalysisManager->FinishNtuple(); // ntupleID: 3

  fAnalysisManager->CreateNtuple("h5_Flux_4004", "Neutrons G4TARC flux");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("tarcflux");
  fAnalysisManager->CreateNtupleDColumn("errstat");
  fAnalysisManager->CreateNtupleDColumn("g4flux");
  fAnalysisManager->CreateNtupleDColumn("g4perp");
  fAnalysisManager->CreateNtupleDColumn("g4fluence");
  fAnalysisManager->CreateNtupleDColumn("g4err");
  fAnalysisManager->CreateNtupleDColumn("rawflux");
  fAnalysisManager->CreateNtupleDColumn("gstep");
  fAnalysisManager->CreateNtupleDColumn("gfl_cyl");
  fAnalysisManager->CreateNtupleDColumn("g4front");
  fAnalysisManager->CreateNtupleDColumn("g4shell");
  fAnalysisManager->CreateNtupleDColumn("g4shellerr");
  fAnalysisManager->FinishNtuple(); // ntupleID: 4


  fAnalysisManager->CreateNtuple("h6_Flux_4005", "Neutrons G4TARC flux");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("tarcflux");
  fAnalysisManager->CreateNtupleDColumn("errstat");
  fAnalysisManager->CreateNtupleDColumn("g4flux");
  fAnalysisManager->CreateNtupleDColumn("g4perp");
  fAnalysisManager->CreateNtupleDColumn("g4fluence");
  //fAnalysisManager->CreateNtupleDColumn("g4zflux");
  fAnalysisManager->CreateNtupleDColumn("g4err");
  //fAnalysisManager->CreateNtupleDColumn("flux5cm");
  //fAnalysisManager->CreateNtupleDColumn("err5cm");
  fAnalysisManager->CreateNtupleDColumn("rawflux");
  fAnalysisManager->CreateNtupleDColumn("gstep");
  fAnalysisManager->CreateNtupleDColumn("gfl_cyl");
  fAnalysisManager->CreateNtupleDColumn("g4front");
  fAnalysisManager->CreateNtupleDColumn("g4shell");
  //fAnalysisManager->CreateNtupleDColumn("g4_shell_err");
  fAnalysisManager->FinishNtuple(); // ntupleID: 5

  fAnalysisManager->CreateNtuple("h7_Created_N_A", "Created Neutrons");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("time");
  fAnalysisManager->CreateNtupleIColumn("particle");
  fAnalysisManager->CreateNtupleDColumn("momentum");
  fAnalysisManager->CreateNtupleDColumn("zmom");
  fAnalysisManager->FinishNtuple(); // ntupleID: 6


  fAnalysisManager->CreateNtuple("h8_Created_N", "Created Neutrons");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("time");
  fAnalysisManager->CreateNtupleDColumn("starte");
  fAnalysisManager->CreateNtupleIColumn("trackid");
  fAnalysisManager->CreateNtupleIColumn("parentid");
  fAnalysisManager->CreateNtupleDColumn("fluxe");
  fAnalysisManager->CreateNtupleDColumn("fluxidx");
  fAnalysisManager->CreateNtupleDColumn("zmom");
  fAnalysisManager->CreateNtupleDColumn("startt");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("e_parent");
  fAnalysisManager->CreateNtupleIColumn("parent");
  fAnalysisManager->CreateNtupleDColumn("step");
  fAnalysisManager->CreateNtupleIColumn("dupli");
  fAnalysisManager->FinishNtuple(); // ntupleID: 7


  fAnalysisManager->CreateNtuple("h9_Rad_Shell_Fluence", "Radial Shell Fluence");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("fluence");
  fAnalysisManager->CreateNtupleDColumn("true_e");
  fAnalysisManager->CreateNtupleDColumn("true_f");
  fAnalysisManager->FinishNtuple(); // ntupleID: 8


  fAnalysisManager->CreateNtuple("h10_Rad_Fluence_Expt_Li_Data", "Lithium Radial Fluence Exptl Data");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("data");
  fAnalysisManager->CreateNtupleDColumn("error");
  fAnalysisManager->FinishNtuple(); // ntupleID: 9

  fAnalysisManager->CreateNtuple("h11_Rad_Fluence_Expt_He3_Data", "He3 Radial Fluence Exptl Data");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("data");
  fAnalysisManager->CreateNtupleDColumn("error");
  fAnalysisManager->FinishNtuple(); // ntupleID: 10


  fAnalysisManager->CreateNtuple("h12_3GeV5_He3_Expt_Data", "Radial Fluence He3 Expt Data");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("data");
  fAnalysisManager->CreateNtupleDColumn("stat_err");
  fAnalysisManager->CreateNtupleDColumn("syst_err");
  fAnalysisManager->FinishNtuple(); // ntupleID: 11

  fAnalysisManager->CreateNtuple("h13_Rad_Fluence_Li", "Radial Fluence Li");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("data");
  fAnalysisManager->CreateNtupleDColumn("stat_err");
  fAnalysisManager->CreateNtupleDColumn("syst_err");
  fAnalysisManager->FinishNtuple(); // ntupleID: 12

  fAnalysisManager->CreateNtuple("h14_Rad_Fluence", "Radial Fluence");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("fluence");
  fAnalysisManager->CreateNtupleDColumn("he_data");
  fAnalysisManager->FinishNtuple(); // ntupleID: 13

  fAnalysisManager->CreateNtuple("h15_Other_ET", "OTHER Time");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("time");
  fAnalysisManager->CreateNtupleDColumn("primary");
  fAnalysisManager->FinishNtuple(); // ntupleID: 14              : 4

/*
  fAnalysisManager->CreateNtuple("h16", "log10(En)");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->FinishNtuple(); // ntupleID: 15 for neutron

  */

  G4cout << "Ntuples created." << G4endl;
}



void G4TARCRunAction::FillRadialExperimentalData(){
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();

  for (G4int ij1 = 0; ij1 < 8; ij1++) {  //  fExptRadiiTables.size(); ij1++){  0~ 41 to 7 ~ 48
    for (std::size_t ij2 = 0; ij2 < fExptRadiiTables[ij1].size(); ij2++){    //   fExptRadiiTables[ij1].size(); ij2++){
      fAnalysisManager->FillNtupleDColumn(9, 0, fExptRadiiTables[ij1][ij2] );  //  converted to mm
      fAnalysisManager->FillNtupleDColumn(9, 1, fExptEnergyBin[ij1]);
      fAnalysisManager->FillNtupleDColumn(9, 2, fExptFluenceTables[ij1][ij2] * 100.0);   // transferring to unit n/cm^2/eV/10^9p
      fAnalysisManager->FillNtupleDColumn(9, 3, fExptErrTables[ij1][ij2] * 100.0);
      fAnalysisManager->AddNtupleRow(9);
    }
  }

  for (std::size_t ij1 = 8; ij1 <= 16; ij1++){ // 8 ~ 49 to 16 ~ 57
    G4int ijE = ij1 - 7;
    for (std::size_t ij2 = 0; ij2 < fExptRadiiTables[ij1].size(); ij2++){    //   fExptRadiiTables[ij1].size(); ij2++){
      fAnalysisManager->FillNtupleDColumn(10, 0, fExptRadiiTables[ij1][ij2] );   // converted to mm
      fAnalysisManager->FillNtupleDColumn(10, 1, fExptEnergyBin[ijE]);
      fAnalysisManager->FillNtupleDColumn(10, 2, fExptFluenceTables[ij1][ij2] * 100.0);  // transferring to unit n/cm^2/eV/10^9p
      fAnalysisManager->FillNtupleDColumn(10, 3, fExptErrTables[ij1][ij2] * 100.0);
      fAnalysisManager->AddNtupleRow(10);
    }
  }

  G4cout << "Experimental data filling complete." << G4endl;

  // This is testing of erasing unused vectors
  std::vector<std::vector<G4double> >().swap(fExptFluenceTables);
  std::vector<std::vector<G4double> >().swap(fExptErrTables);
  // This is the end of the test
}
