/*********************************************************
 * @file        G4TARCRunAction.cxx
 * @author      Abhijit Bhattacharyya
 * @brief       This is for the run Action
 ********************************************************/
#include "G4TARCRunAction.hh"

G4TARCRunAction::G4TARCRunAction() : G4UserRunAction() {
  //    auto fAnalysisManager = G4AnalysisManager::Instance();
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
  const G4TARCRun*  tarcRun = static_cast<const G4TARCRun*>(thisRun);
  FillRadialExperimentalData(tarcRun);
  G4int eventN = tarcRun->GetNumberOfEvents();

  G4cout << "Total Number of events: " << eventN << G4endl;
  G4double ExitFlux = tarcRun->GetExitingFlux();
  G4double ExitEnergy = tarcRun->GetExitingEnergy();
  G4double ExitCheckFlux = tarcRun->GetExitingCheckFlux();

  if (IsMaster()) {
    NeutronFluxHistogram(eventN, tarcRun);
    RadialFluxHistogram(eventN, tarcRun);
    G4cout << "Exit flux: " << ExitFlux << ",   Exit Energy: " << ExitEnergy << ",   Exit CheckFlux: " << ExitCheckFlux << G4endl;
    G4cout << "Integral Neutron Flux at 45.6 cms:  " << tarcRun->fIntegral_flux_46cm << " Integral EFLUX at 45.6 cms:  " << tarcRun->fTARC_Integral_Eflux_46cm << G4endl;
  }

  auto fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->Write();
  fAnalysisManager->CloseFile();
}


void G4TARCRunAction::DefineShellBlocks() {
  fRefShellNumber = fRadiusReference.size();
  fRefShellThickness = 2.0 * mm;
  for (G4int ii = 0; ii < fRefShellNumber; ii++) {
    G4double radThis = fRadiusReference[ii] ;
    fOuterRadiusofShell.push_back(radThis);
    radThis = fRadiusReference[ii]  - fRefShellThickness;
    fInnerRadiusofShell.push_back(radThis);
  }

  fLocal_Energy_Integral.resize(4, 0.0);

  fRefShellOuterRad = 456.0 * mm;
  fRefShellInnerRad = fRefShellOuterRad - fRefShellThickness;
  fTestSphereRadius = 45.6 * cm;
  fRadHole = 32.0 * cm;
  fLenCyl = 150.0 * mm;
  fTestSphereVolume = (4.0 / 3.0) * CLHEP::pi * (fRadHole * fRadHole * fRadHole) ;  // This is in mm3
  fTestSphereSurfaceArea = 4.0 * CLHEP::pi * (fTestSphereRadius * fTestSphereRadius) ; //  / cm2;
  fTestShellVol          = (4.0 / 3.0) * CLHEP::pi * (std::pow(fRefShellOuterRad, 3.0) - std::pow(fRefShellInnerRad, 3.0));

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

void G4TARCRunAction::FillRadialExperimentalData(const G4TARCRun* tarcRun){
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
  for (G4int ij1 = 0; ij1 < 8; ij1++) {  //  fExptRadiiTables.size(); ij1++){  0~ 41 to 7 ~ 48
    for (std::size_t ij2 = 0; ij2 < tarcRun->fExptRadiiTables[ij1].size(); ij2++){    //   fExptRadiiTables[ij1].size(); ij2++){
      fAnalysisManager->FillNtupleDColumn(9, 0, tarcRun->fExptRadiiTables[ij1][ij2] );  //  converted to mm
      fAnalysisManager->FillNtupleDColumn(9, 1, tarcRun->fExptEnergyBin[ij1]);
      fAnalysisManager->FillNtupleDColumn(9, 2, tarcRun->fExptFluenceTables[ij1][ij2] * 100.0);   // transferring to unit n/cm^2/eV/10^9p
      fAnalysisManager->FillNtupleDColumn(9, 3, tarcRun->fExptErrTables[ij1][ij2] * 100.0);
      fAnalysisManager->AddNtupleRow(9);
    }
  }

  for (std::size_t ij1 = 8; ij1 <= 16; ij1++){ // 8 ~ 49 to 16 ~ 57
    G4int ijE = ij1 - 7;
    for (std::size_t ij2 = 0; ij2 < tarcRun->fExptRadiiTables[ij1].size(); ij2++){    //   fExptRadiiTables[ij1].size(); ij2++){
      fAnalysisManager->FillNtupleDColumn(10, 0, tarcRun->fExptRadiiTables[ij1][ij2] );   // converted to mm
      fAnalysisManager->FillNtupleDColumn(10, 1, tarcRun->fExptEnergyBin[ijE]);
      fAnalysisManager->FillNtupleDColumn(10, 2, tarcRun->fExptFluenceTables[ij1][ij2] * 100.0);  // transferring to unit n/cm^2/eV/10^9p
      fAnalysisManager->FillNtupleDColumn(10, 3, tarcRun->fExptErrTables[ij1][ij2] * 100.0);
      fAnalysisManager->AddNtupleRow(10);
    }
  }
  G4cout << "Experimental data filling complete." << G4endl;
}



void G4TARCRunAction::NeutronFluxHistogram(G4int fNevents, const G4TARCRun* tarcRun){
  auto fAnalysisManager = G4AnalysisManager::Instance();
  G4double fAbsolute_TotalFlux = (tarcRun->fTotalFlux *  1.0e9 / (G4double)fNevents) / (fTestSphereSurfaceArea); // per cm^2
  for (G4int ij1 = 0; ij1 < tarcRun->fMaxTestFluxData; ij1++){
    G4double fMeanEnergy   = 0.5 * (tarcRun->fFlux_Energy[ij1 + 1] + tarcRun->fFlux_Energy[ij1]);   // eV
    G4double fAbsFlux      = (tarcRun->fFlux[ij1] *  (1.0e9 / (G4double)fNevents)) / (fTestSphereSurfaceArea);     // per cm^2
    G4double fBinWidth     = std::abs(tarcRun->fFlux_Energy[ij1 + 1] - tarcRun->fFlux_Energy[ij1]);   //eV
    G4double fAbsFluxPerp  = fMeanEnergy * (((tarcRun->fCos_Flux[ij1] * (1.0e9 / (G4double)fNevents)) / (fTestSphereSurfaceArea)) / fBinWidth);   // per cm^2
    G4double fAbsEFlux     = (tarcRun->fEFlux[ij1] * ( 1.0e9 / (G4double)fNevents)) / (fTestSphereSurfaceArea);                                   // per cm^2
    G4double fAbsFluence   = fMeanEnergy * ( (1.0e9 / (G4double)fNevents) * ((tarcRun->fFluence_Step[ij1] / 10.0) / fTestSphereVolume) / fBinWidth);       // (mm/10)/cm^3-> per cm^2
    G4double fAbsFluenceShell = fMeanEnergy * ( (1.0e9 / (G4double)fNevents) * ((tarcRun->fFluence_Step_Shell[ij1] / 10.0) / fTestShellVol) / fBinWidth);  // (mm/10)/cm^3-> per cm^2
    G4double fAbsFluenceShellErr = (tarcRun->fFluence_Step_Shell[ij1] > 0.0)
                                             ? (std::pow(tarcRun->fFluence_Step_Shell[ij1], 0.5) / tarcRun->fFluence_Step_Shell[ij1]) * fAbsFluenceShell
                                             : 0.0;
    G4double fAbsErr = (tarcRun->fFlux[ij1] != 0.0)
                                              ? (std::pow(tarcRun->fFlux[ij1], 0.5) / tarcRun->fFlux[ij1]) * fAbsFlux
                                              : 0.0;

    fTARC_helium   += tarcRun->fFlux_Data[ij1] ;
    fTARC_helium_E += tarcRun->fFlux_Data[ij1] * fMeanEnergy;

    fLocal_Energy_Integral[2] += fAbsFlux * fMeanEnergy;
    fLocal_Energy_Integral[3] += tarcRun->fFlux_Data[ij1] * fMeanEnergy;

    fAnalysisManager->FillNtupleDColumn(3, 0, fMeanEnergy);
    fAnalysisManager->FillNtupleDColumn(3, 1, tarcRun->fFlux_Data[ij1]  * 100.0);
    fAnalysisManager->FillNtupleDColumn(3, 2, tarcRun->fFlux_Syst_Err[ij1] * 100.0);
    fAnalysisManager->FillNtupleDColumn(3, 3, fAbsFlux );
    fAnalysisManager->FillNtupleDColumn(3, 4, fAbsFluxPerp);
    fAnalysisManager->FillNtupleDColumn(3, 5, fAbsFluence );
    fAnalysisManager->FillNtupleDColumn(3, 6, fAbsErr);
    fAnalysisManager->FillNtupleDColumn(3, 7, tarcRun->fFlux[ij1] );
    fAnalysisManager->FillNtupleDColumn(3, 8, tarcRun->fEflux_Data[ij1]);         // eflux_data
    fAnalysisManager->FillNtupleDColumn(3, 9, fAbsEFlux);
    fAnalysisManager->FillNtupleDColumn(3, 10, tarcRun->fFluence_Step[ij1]);
    fAnalysisManager->FillNtupleDColumn(3, 11, 0.0);                      // abs fluence cyl
    fAnalysisManager->FillNtupleDColumn(3, 12, 0.0);                     // abs fluence front
    fAnalysisManager->FillNtupleDColumn(3, 13, fAbsFluenceShell);           // abs fluence Shell
    fAnalysisManager->FillNtupleDColumn(3, 14, fAbsFluenceShellErr);           // abs fluence Shell
    fAnalysisManager->AddNtupleRow(3);

  }

  for (G4int ij1 = 0; ij1 < tarcRun->fMaxFluenceData; ij1++){
    G4double fMeanLowEnergy   = std::exp(0.5 * (std::log(tarcRun->fFlux_Low_Energy[ij1 + 1]) + std::log(tarcRun->fFlux_Low_Energy[ij1])));     // eV
    G4double fAbsLowFlux      = (tarcRun->fFlux_Low_Data[ij1] * ( 1.0e9 / (G4double)fNevents)) / fTestSphereSurfaceArea;
    G4double fBinWidth      = tarcRun->fFlux_Low_Energy[ij1 + 1] - tarcRun->fFlux_Low_Energy[ij1];
    G4double fAbsLowFluxPerp = fMeanLowEnergy * (((tarcRun->fCos_Low_Flux[ij1] * (1e9 / (G4double)fNevents)) / fTestSphereSurfaceArea) / fBinWidth);
    G4double fAbsLowFluence   =  ( 1.0e9 / (G4double)fNevents) * ((tarcRun->fLow_Fluence_Step[ij1] / 10.0) / fTestSphereVolume);
    G4double fAbsLowFluenceShell = fMeanLowEnergy * ((1.0e9 / (G4double)fNevents) * ((tarcRun->fLow_Fluence_Step_Shell[ij1] / 10.0) / fTestShellVol) / fBinWidth);

    G4double fAbsLowFluenceShellError = (tarcRun->fLow_Fluence_Step_Shell[ij1] > 0)
                 ? (std::pow(tarcRun->fLow_Fluence_Step_Shell[ij1], 0.5) / tarcRun->fLow_Fluence_Step_Shell[ij1]) * fAbsLowFluenceShell
                 : 0.0;
    G4double fAbsLowError = (tarcRun->fFlux_Low_Data[ij1] != 0.0)
                 ? (std::pow(tarcRun->fFlux_Low_Data[ij1], 0.5) / tarcRun->fFlux_Low_Data[ij1]) * fAbsLowFlux
                 : 0.0;
    fTARC_Integral += tarcRun->fFlux_Low_Data[ij1] ;
    fTARC_Integral_E += tarcRun->fFlux_Low_Data[ij1] * fMeanLowEnergy;

    fAnalysisManager->FillNtupleDColumn(4, 0, fMeanLowEnergy);
    fAnalysisManager->FillNtupleDColumn(4, 1, tarcRun->fFlux_Low_Data_in[ij1] * 100.0);
    fAnalysisManager->FillNtupleDColumn(4, 2, tarcRun->fFlux_Low_Syst_Err[ij1] * 100.0);
    fAnalysisManager->FillNtupleDColumn(4, 3, fAbsLowFlux);
    fAnalysisManager->FillNtupleDColumn(4, 4, fAbsLowFluxPerp);
    fAnalysisManager->FillNtupleDColumn(4, 5, fAbsLowFluence);
    fAnalysisManager->FillNtupleDColumn(4, 6, fAbsLowError);
    fAnalysisManager->FillNtupleDColumn(4, 7, tarcRun->fFlux_Low[ij1]);
    fAnalysisManager->FillNtupleDColumn(4, 8, tarcRun->fLow_Fluence_Step[ij1]);         // low_fluence_step
    fAnalysisManager->FillNtupleDColumn(4, 9, 0.0);        // abs low fluence cyl
    fAnalysisManager->FillNtupleDColumn(4, 10, 0.0);           // abs low fluence front
    fAnalysisManager->FillNtupleDColumn(4, 11, fAbsLowFluenceShell);           // abs low fluence shell
    fAnalysisManager->FillNtupleDColumn(4, 11, fAbsLowFluenceShellError);           // abs low fluence shell error
    fAnalysisManager->AddNtupleRow(4);
  }

  for (G4int ij1 = 0; ij1 < tarcRun->fMaxFluenceData; ij1++){
    G4double fMeanLithiumEnergy   = std::exp(0.5 * (std::log(tarcRun->fFlux_Lithium_Energy[ij1 + 1]) + std::log(tarcRun->fFlux_Lithium_Energy[ij1])));
    G4double fAbsLithiumFlux      = (tarcRun->fFlux_Lithium_Data[ij1] * (1.0e9 / (G4double)fNevents)) / fTestSphereSurfaceArea;
    //  G4double fBinWidth = fFlux_Low_Energy[ij1 + 1] - fFlux_Low_Energy[ij1];
    G4double fBinWidth = tarcRun->fFlux_Lithium_Energy[ij1 + 1] - tarcRun->fFlux_Lithium_Energy[ij1];
    G4double fAbsLithiumFluxPerp = fMeanLithiumEnergy * (((tarcRun->fCos_Lithium_Flux[ij1] * (1.0e9 / (G4double)fNevents)) / fTestSphereSurfaceArea) / fBinWidth);
    G4double fAbsLithiumFluence   =   (1.0e9 / (G4double)fNevents) * ((tarcRun->fLithium_Fluence_Step[ij1] / 10.0) / fTestSphereVolume);
    G4double fAbsLithiumFluenceShell = fMeanLithiumEnergy * ( (1.0e9 / (G4double)fNevents) * ((tarcRun->fLithium_Fluence_Step_Shell[ij1] / 10.0) / fTestShellVol) / fBinWidth);

    G4double fAbsErrorLi =(tarcRun->fFlux_Low_Data[ij1] != 0.0)
          ? (std::pow(tarcRun->fFlux_Lithium_Data[ij1], 0.5) / tarcRun->fFlux_Lithium_Data[ij1]) * fAbsLithiumFlux
          : 0.0;

    fTARC_lithium   += tarcRun->fFlux_Lithium_Data[ij1];
    fTARC_lithium_E += tarcRun->fFlux_Lithium_Data[ij1] * fMeanLithiumEnergy;

    fAnalysisManager->FillNtupleDColumn(5, 0, fMeanLithiumEnergy);
    fAnalysisManager->FillNtupleDColumn(5, 1, tarcRun->fFlux_Lithium_Data[ij1] * 100.0);
    fAnalysisManager->FillNtupleDColumn(5, 2, tarcRun->fFlux_Low_Syst_Err[ij1] * 100.0);
    fAnalysisManager->FillNtupleDColumn(5, 3, fAbsLithiumFlux);
    fAnalysisManager->FillNtupleDColumn(5, 4, fAbsLithiumFluxPerp);
    fAnalysisManager->FillNtupleDColumn(5, 5, fAbsLithiumFluence);
    fAnalysisManager->FillNtupleDColumn(5, 6, fAbsErrorLi);
    fAnalysisManager->FillNtupleDColumn(5, 7, tarcRun->fLithium_Flux[ij1] * 100.0);
    fAnalysisManager->FillNtupleDColumn(5, 8, tarcRun->fLithium_Fluence_Step[ij1]);
    fAnalysisManager->FillNtupleDColumn(5, 9, 0.0);
    fAnalysisManager->FillNtupleDColumn(5, 10, 0.0);
    fAnalysisManager->FillNtupleDColumn(5, 11, fAbsLithiumFluenceShell);
    fAnalysisManager->AddNtupleRow(5);
  }

  for (G4int i1 = 0; i1 < tarcRun->fMaxRadCount; i1++){
    for (G4int i2 = 0; i2 < tarcRun->fMaxRadCount; i2++) {
      fAnalysisManager->FillNtupleDColumn(13, 0, tarcRun->fRadList[i1]);
      fAnalysisManager->AddNtupleRow(13);
    }
  }
}

void G4TARCRunAction::RadialFluxHistogram(G4int fNevents, const G4TARCRun* aRun){
  auto fAnalysisManager = G4AnalysisManager::Instance();
  const G4TARCRun* tarcRun = static_cast<const G4TARCRun*>(aRun);

  for (G4int ijk1 = 0; ijk1 < fRefShellNumber; ijk1++) {  // ijk1 < fMaxTestFluxData; ijk1++){
    //G4cout << ijk1 << " / " << fRefShellNumber << " RO " << fOuterRadiusofShell[ijk1] << " RI " << fInnerRadiusofShell[ijk1] << "  " << G4endl;
    G4double shellVol = (4.0 / 3.0) * CLHEP::pi * (std::pow(fOuterRadiusofShell[ijk1], 3.0) - std::pow(fInnerRadiusofShell[ijk1], 3.0)); // mm^3
    G4double radL = 0.5 * (fOuterRadiusofShell[ijk1] + fInnerRadiusofShell[ijk1]);                                                      // mm
    for (G4int ijk2 = 0; ijk2 < tarcRun->fMaxRadCount; ijk2++){
      G4double fBinTmpWidth = tarcRun->fLithium_Radial_Energy_Upper[ijk2] - tarcRun->fLithium_Radial_Energy_Lower[ijk2];
      //G4cout << " ijk2 " << ijk2 << " R " << radL << "  F  " << fRadialFluenceStep[ijk1][ijk2] << G4endl;

      if (tarcRun->fRadialFluenceStep[ijk1][ijk2] != 0.0){
        G4double fAbsRadFluence = tarcRun->fLithium_Radial_Mean[ijk2]*(100.0 * (1.0e9 / (G4double)fNevents) * ((tarcRun->fRadialFluenceStep[ijk1][ijk2]) / shellVol) / fBinTmpWidth); // radii are in mm
        G4double fAbsRadFluenceTrueMean = tarcRun->fLithium_Radial_True_Mean[ijk2] * (100.0 * ( 1.0e9 / (G4double)fNevents) * tarcRun->fRadialFluenceStep[ijk1][ijk2] / shellVol / fBinTmpWidth);
        fAnalysisManager->FillNtupleDColumn(8, 0, radL);    // rad in mm
        fAnalysisManager->FillNtupleDColumn(8, 1, tarcRun->fLithium_Radial_Mean[ijk2]);
        fAnalysisManager->FillNtupleDColumn(8, 2, fAbsRadFluence);
        fAnalysisManager->FillNtupleDColumn(8, 3, tarcRun->fLithium_Radial_True_Mean[ijk2]);
        fAnalysisManager->FillNtupleDColumn(8, 4, fAbsRadFluenceTrueMean);
        fAnalysisManager->AddNtupleRow(8);
      }
    }
  }
}
