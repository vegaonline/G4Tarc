/*********************************************************
 * @file        G4TARCRunAction.cxx
 * @author      Abhijit Bhattacharyya
 * @brief       This is for the run Action
 ********************************************************/
#include "G4TARCRunAction.hh"

//G4TARCRunAction::G4TARCRunAction(G4TARCDetectorConstruction* det, G4TARCPrimaryGeneratorAction* prim)
//: G4UserRunAction(), fDetector(det), fPrimary(prim), fHisto(0){
G4TARCRunAction::G4TARCRunAction(): G4UserRunAction(){
  fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->OpenFile(fAnalysisFileName);
  fHistoBooked = false;

  if (!fHistoBooked) {
    BookHistogram();
    CreateTuples();
  }
}

G4TARCRunAction::~G4TARCRunAction() {

}

void G4TARCRunAction::BookHistogram() {
  fHistoBooked = true;
  fAnalysisManager->SetFirstHistoId(1);
  fAnalysisManager->CreateH1("TH1","Gamma Edep /keV", 10000, 0.0, 1000*keV);
  fAnalysisManager->CreateH1("TH2","Neutron ener vs. 1/mom /eV", 100000, -1.0, 100.0);   // 100000, 0., 1000000.);
  fAnalysisManager->CreateH1("TH3","Electron Edep /keV", 100000, 0.0, 1000*keV);
  fAnalysisManager->CreateH1("TH4","Positron Edep /keV", 10000, 0.0, 1000*keV);
  fAnalysisManager->CreateH1("TH5","Other Edep /keV", 1000, 0.0, 1000*keV);
  fAnalysisManager->CreateH1("TH6","Particle Stack", 1000, 0.5, 12.5);
  fAnalysisManager->CreateH1("TH7","Neutrons/event", 30, 0.0, 30.0);
  fAnalysisManager->CreateH1("TH8","Protons/event", 30, 0.0, 30.0);
  fAnalysisManager->CreateH2("Neutron_Energy_Time","Neutron Energy vs. Time", 100, -4.0, 4.0, 100, -4.0, 8.0);
  //fAnalysisManager->CreateH2("Neutron_Energy_Mom","Neutron Energy vs. Mometum", 100, 0.0, 1500.0, 100, 0.0, 1500.0);
  fAnalysisManager->CreateH2("Other_Particle_Energy_Time","OTHER particle Energy vs. Time", 100, -4.0, 4.0, 100, -4.0, 6.0);



  /*
  fAnalysisManager->CreateH1("Protons_per_Event", "Protons/event", 50, 0.0, 100.0);
  fAnalysisManager->CreateH1("Neutrons_per_Event", "Neutrons/event", 50, 0.0, 100.0);
  fAnalysisManager->CreateH1("Neutron_Energy", "Neutron Energy vs 1/mom /eV", 100000, 0.0, 1000000.0);
  fAnalysisManager->CreateH1("proton_Energy", "Proton Energy Deposition/keV", 1000, 0.0, 1000.0 * keV);
  fAnalysisManager->CreateH1("Particle_Stack", "Particle Stack", 12, 0.5, 12.5);
  fAnalysisManager->CreateH1("log10En", "Log10 Energy (eV) of neutron", 100, -9.0, 9.0);
  fAnalysisManager->CreateH2("Neutron_Energy_Time", "Neutron Energy vs Time", 100, -1.0, 4.0, 100, -2.0, 7.0);
  fAnalysisManager->CreateH2("Other_particle_Energy_Time", "Other particle Energy vs Time", 100, -1.0, 4.0, 100, -2.0, 7.0);
  */
}


void G4TARCRunAction::CreateTuples(){
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
  fAnalysisManager->FinishNtuple(); // ntupleID: 0 - filled

  fAnalysisManager->CreateNtuple("h2_N_ET", "Neutron Time");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("time");
  fAnalysisManager->CreateNtupleDColumn("primary");
  fAnalysisManager->FinishNtuple(); // ntupleID: 1

  fAnalysisManager->CreateNtuple("h3_N_Exiting", "Neutrons Exiting");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->FinishNtuple(); // ntupleID: 2 - filled

  fAnalysisManager->CreateNtuple("h4_Flux_4002", "Neutrons G4TARC flux");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("tarcflux");
  fAnalysisManager->CreateNtupleDColumn("errstat");
  fAnalysisManager->CreateNtupleDColumn("g4flux");
  fAnalysisManager->CreateNtupleDColumn("g4perp");
  fAnalysisManager->CreateNtupleDColumn("gfluence");
  fAnalysisManager->CreateNtupleDColumn("g4err");
  fAnalysisManager->CreateNtupleDColumn("rawflux");
  fAnalysisManager->CreateNtupleDColumn("trceflux");
  fAnalysisManager->CreateNtupleDColumn("g4eflux");
  fAnalysisManager->CreateNtupleDColumn("gstep");
  fAnalysisManager->CreateNtupleDColumn("gfl_cyl");
  fAnalysisManager->CreateNtupleDColumn("g4front");
  // fAnalysisManager->CreateNtupleDColumn("g4_shell");
  fAnalysisManager->FinishNtuple(); // ntupleID: 3

  fAnalysisManager->CreateNtuple("h5_Flux_4004", "Neutrons G4TARC flux");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("tarcflux");
  fAnalysisManager->CreateNtupleDColumn("errstat");
  fAnalysisManager->CreateNtupleDColumn("g4flux");
  fAnalysisManager->CreateNtupleDColumn("g4perp");
  fAnalysisManager->CreateNtupleDColumn("gfluence");
  fAnalysisManager->CreateNtupleDColumn("g4err");
  fAnalysisManager->CreateNtupleDColumn("rawflux");
  //fAnalysisManager->CreateNtupleDColumn("gstep");
  //fAnalysisManager->CreateNtupleDColumn("gfl_cyl");
  //fAnalysisManager->CreateNtupleDColumn("g4front");
  // fAnalysisManager->CreateNtupleDColumn("g4_shell");
  fAnalysisManager->FinishNtuple(); // ntupleID: 4


  fAnalysisManager->CreateNtuple("h6 Flux 4005", "Neutrons Test15 flux");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("tarcflux");
  //fAnalysisManager->CreateNtupleDColumn("errstat");
  fAnalysisManager->CreateNtupleDColumn("errsyst");
  fAnalysisManager->CreateNtupleDColumn("g4flux");
  fAnalysisManager->CreateNtupleDColumn("g4perp");
  fAnalysisManager->CreateNtupleDColumn("gfluence");
  //fAnalysisManager->CreateNtupleDColumn("g4zflux");
  fAnalysisManager->CreateNtupleDColumn("g4err");
  //fAnalysisManager->CreateNtupleDColumn("flux5cm");
  //fAnalysisManager->CreateNtupleDColumn("err5cm");
  fAnalysisManager->CreateNtupleDColumn("rawflux");
  fAnalysisManager->CreateNtupleDColumn("gstep");
  //fAnalysisManager->CreateNtupleDColumn("gfl_cyl");
  //fAnalysisManager->CreateNtupleDColumn("g4front");
  //fAnalysisManager->CreateNtupleDColumn("g4_shell");
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
  fAnalysisManager->CreateNtupleIColumn("step");
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


  fAnalysisManager->CreateNtuple("h12_3.5GeV_He3_Expt_Data", "Radial Fluence He3 Expt Data");
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
  fAnalysisManager->FinishNtuple(); // ntupleID: 14

  fAnalysisManager->CreateNtuple("h16", "log10(En)");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->FinishNtuple(); // ntupleID: 15 for neutron

  G4cout << "Ntuples created." << G4endl;
}


void G4TARCRunAction::BeginOfRunAction( const G4Run* aRun ) {
  auto id = aRun->GetRunID();
  G4cout << "Run # " << id << " starts." << G4endl;
  G4NuclearLevelData::GetInstance();
  (G4TARCHistoManager::GetPointer())->BeginOfRun();

#ifdef G4VIS_USE
  auto UI = G4UImanager::GetUIpointer();
  auto pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager) {
    UI->ApplyCommand("/vis/scene/notifyHandlers");    // this crashes the code .....
  }
#endif
}


void G4TARCRunAction::EndOfRunAction( const G4Run* aRun ){
  //auto analysisManager = G4AnalysisManager::Instance();
  G4cout << " RunAction: End of run actions for # " << aRun->GetRunID() << " is started" << G4endl;

  (G4TARCHistoManager::GetPointer())->EndOfRun();
  #ifdef G4VIS_USE
    if (G4VVisManager::GetConcreteInstance())
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  #endif
}
