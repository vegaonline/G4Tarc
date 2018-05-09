#include "G4TARCHistoManager.hh"

G4TARCHistoManager *  G4TARCHistoManager::fHistoManager = 0;

G4TARCHistoManager *  G4TARCHistoManager::GetPointer() {
  if ( !fHistoManager ) {
    static G4TARCHistoManager manager;
    fHistoManager = &manager;
  }
  return fHistoManager;
}

G4TARCHistoManager::G4TARCHistoManager()
//  fHisto( 0 ),
//  fHistoBooked(false),
//  fPrimaryDef( 0 ),
//  fProton( 0 ),
// fNHisto( 25 ),
//  fNeutron( 0 ),
//  fEdepMax( fMaxEVal ),
//  fLength( fMaxLVal ),
//  fHLength(0.5 * fMaxLVal),
//  fPrimaryKineticEnergy ( fEVal0 ),
//  fVerbose( 1 ),
//  fNtuple_full(false),
//  fNBinsE( fMaxBin ),
//  fProtonIN(0),
//  fNCountTotal(0),
//  fVirtVol(0),
//  fNSlices( fMaxSlices )
{
    fHistoBooked = false;
    fEdepMax = fMaxEVal;
    fLength = fMaxLVal;
    fHLength = 0.5 * fMaxLVal;
    fPrimaryKineticEnergy = fEVal0;
    fVerbose = 1;
    fNtuple_full = false;
    fNBinsE = fMaxBin;
    fNSlices = fMaxSlices;
    //fHisto    = new G4TARCHisto();
    //fHisto->SetVerbose(fVerbose);
    //fNeutron = G4Neutron::Neutron();
    //fProton = G4Proton::Proton();
}

G4TARCHistoManager::~G4TARCHistoManager() {
  // delete fHisto;
  std::vector<G4int>().swap(fFluxTableList);
  std::vector<G4double>().swap(fExptEnergyBin);
  std::vector<G4double>().swap(fFluxRadTables);
  std::vector<G4double>().swap(fRadList);
  std::vector<G4double>().swap(fOuterRadiusofShell);
  std::vector<G4double>().swap(fInnerRadiusofShell);
  std::vector<G4double>().swap(fFlux);
  std::vector<G4double>().swap(fFlux_Energy);
  std::vector<G4double>().swap(fFlux_Data);
  std::vector<G4double>().swap(fFlux_Syst_Err);
  std::vector<G4double>().swap(fFlux_Energy_in);
  std::vector<G4double>().swap(fFlux_Data_in);
  std::vector<G4double>().swap(fFlux_Syst_Err_in);
  std::vector<G4double>().swap(fFlux_Low);
  std::vector<G4double>().swap(fFlux_Low_Radius);
  std::vector<G4double>().swap(fFlux_Low_Energy);
  std::vector<G4double>().swap(fFlux_Low_Energy_in);
  std::vector<G4double>().swap(fFlux_Low_Data);
  std::vector<G4double>().swap(fFlux_Low_Syst_Err);
  std::vector<G4double>().swap(fFlux_Lithium);

  std::vector<G4double>().swap(fFlux_Lithium_Radius);
  std::vector<G4double>().swap(fFlux_Lithium_Energy);
  std::vector<G4double>().swap(fFlux_Lithium_Energy_in);
  std::vector<G4double>().swap(fFlux_Lithium_Data);
  std::vector<G4double>().swap(fFlux_Lithium_Syst_Err);
  std::vector<G4double>().swap(fFluence1D);
  std::vector<G4double>().swap(fFluence_Radius);
  std::vector<G4double>().swap(fFluence_Energy);
  std::vector<G4double>().swap(fFluence_Data);
  std::vector<G4double>().swap(fFluence_Syst_Err);
  std::vector<G4double>().swap(fEflux_Data);
  std::vector<G4double>().swap(fFine_Energy);
  std::vector<G4double>().swap(fENflux);
  std::vector<G4double>().swap(fNeutflux);
  std::vector<G4double>().swap(fFluence_Spectrum);
  std::vector<G4double>().swap(fLithium_Radial_Energy_Lower);
  std::vector<G4double>().swap(fLithium_Radial_Energy_Upper);

  std::vector<G4double>().swap(fLithium_Radial_Mean);
  std::vector<G4double>().swap(fLithium_Fluence_Step_Shell);
  std::vector<G4double>().swap(fLithium_Fluence_Step);
  std::vector<G4double>().swap(fLow_Fluence_Step_Shell);
  std::vector<G4double>().swap(fFluence_Step_Shell);
  std::vector<G4double>().swap(fLithium_Flux);
  std::vector<G4double>().swap(fCos_Lithium_Flux);
  std::vector<G4double>().swap(fLow_Flux);
  std::vector<G4double>().swap(fCos_Low_Flux);
  std::vector<G4double>().swap(fCos_Flux);
  std::vector<G4double>().swap(fEFlux);
  std::vector<G4double>().swap(fFluence_Cyl);
  std::vector<G4double>().swap(fFluence_step);
  std::vector<G4double>().swap(fLow_Fluence_step);

  std::vector<std::vector<G4double> >().swap(fRadial_Fluence_Step);
  std::vector<std::vector<G4double> >().swap(fExptRadiiTables);
  std::vector<std::vector<G4double> >().swap(fExptFluenceTables);
  std::vector<std::vector<G4double> >().swap(fExptErrTables);
  std::vector<std::vector<G4double> >().swap(fExptEnergyTables);
  std::vector<std::vector<G4double> >().swap(fExptFluxTables);
  std::vector<std::vector<G4double> >().swap(fExptFluxErrTables);
  std::vector<std::vector<G4double> >().swap(fFlux_Radius);


}

void G4TARCHistoManager::DefineShellBlocks() {
  fHalfXBlockB           =     0.5 * 300 * mm;
  fHalfYBlockB           =     0.5 * 300 * mm;
  fHalfZBlockB           =     0.5 * 600 * mm;
  fHalfXVBox             =     0.5 * 150 * mm;
  fHalfYVBox             =     0.5 * 150 * mm;
  fHalfZVBox             =     0.5 * 300 * mm;
  fNewHalfZProt          =     0.5 * ((2.0 * fHalfZBlockB) / 3.0);
  fZposProt              = -fHalfZBlockB + fNewHalfZProt;
  fShellThickness        =     50.0 * mm;
  fMinInnerRadiusofShell =    10.0 * mm;
  fMaxOuterRadiusofShell =   1500.0 * mm;
  fInnerRadProtonShell   =     0.0 * mm;   //
  fOuterRadProtonShell   =   300.0 * mm;   // These two were thought as a spherical 4Pi measurement for Proton
  fShellNumber           = (G4int)((fMaxOuterRadiusofShell - fMinInnerRadiusofShell) / fShellThickness + 0.5);
  G4double tmp1          = fMaxOuterRadiusofShell;
  G4double tmp2          =     0.0;
  for (G4int ii = 0; ii < fShellNumber; ii++) {
    tmp2 = tmp1 - fShellThickness;
    fInnerRadiusofShell.push_back(tmp2);
    fOuterRadiusofShell.push_back(tmp1);
    tmp1 = tmp2;
  }
  fTestSphereRadius = 45.6 * cm;
  fTestSphereVolume = (4.0 / 3.0) * CLHEP::pi * (fTestSphereRadius * fTestSphereRadius * fTestSphereRadius);
  fTestSphereSurfaceArea = 4.0 * CLHEP::pi * (fTestSphereRadius * fTestSphereRadius);
  fTestShellVol          = (4.0 / 3.0) * CLHEP::pi * (std::pow(fMaxOuterRadiusofShell, 3.0) - std::pow(fMinInnerRadiusofShell, 3.0));
  fRadHole = 32.0 * mm;
  fLenCyl  = 150.0 * mm;
}


void G4TARCHistoManager::BookHistogram() {
  fHistoBooked = true;
  fAnalysisManager->SetFirstHistoId(1);
  fAnalysisManager->CreateH1("TH1","Gamma Edep /keV", 1000, 0.0, 1000*keV);
  fAnalysisManager->CreateH1("TH2","Neutron ener vs. 1/mom /eV", 100000, -1.0, 250.0);   // 100000, 0., 1000000.);
  fAnalysisManager->CreateH1("TH3","Electron Edep /keV", 10000, 0.0, 1000*keV);
  fAnalysisManager->CreateH1("TH4","Positron Edep /keV", 10000, 0.0, 1000*keV);
  fAnalysisManager->CreateH1("TH5","Other Edep /keV", 1000, 0.0, 1000*keV);
  fAnalysisManager->CreateH1("TH6","Particle Stack", 1000, 0.5, 12.5);
  fAnalysisManager->CreateH1("TH7","Neutrons/event", 30, 0.0, 30.0);
  fAnalysisManager->CreateH1("TH8","Protons/event", 30, 0.0, 30.0);
  fAnalysisManager->CreateH2("Neutron_Energy_Time","Neutron Energy vs. Time", 100, -4.0, 4.0, 100, -4.0, 7.0);
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


void G4TARCHistoManager::CreateTuples(){
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


void G4TARCHistoManager::ReadExperimentalDataFromFile(G4String& exptFileName){
  std::ifstream exptIN(exptFileName, std::ios::in);
  G4String lineIN;
  unsigned NCount = 0, restCount = 0, file0 = 0, iTableNum = 0;
  G4bool readPara = false;
  G4double v1, v2, v3, v4, v5, v6, v7, v8, v9, v10;
  G4bool isFlux = false;
  std::vector<G4double> tmpV1;
  std::vector<G4double> tmpV2;
  std::vector<G4double> tmpV3;
  fMaxFluxData = -999999.99;
  fMaxFluenceData = -999999.99;
  fMaxTestFluxData = -999999.99;
  //G4int fileCount = 0;
  while (getline(exptIN, lineIN)){
    lineIN = std::regex_replace(lineIN, std::regex("^ +| +$|( ) +"), "$1");
    if (lineIN.size() > 1 && lineIN.find("#", 0, 1) != std::string::npos){  // if the line starts with # sign
      std::size_t found1 = lineIN.find("Table");
      readPara = false;
      G4String tableNum = (lineIN.substr(found1 + 5, 3));
      tableNum = std::regex_replace(tableNum, std::regex("^ +| +$|( ) +"), "$1"); // stripping extra spaces
      iTableNum = std::atoi(tableNum);
      isFlux = (std::find( fFluxTableList.begin(), fFluxTableList.end(), iTableNum) != fFluxTableList.end());
      if (found1 != std::string::npos){
        file0 = (iTableNum == 0) ? 1 : 0;
      }
      lineIN="";
    } else if (lineIN.size() > 1 && lineIN.find("#", 0, 1) == std::string::npos){   // the line does not start with # symbol
      if (lineIN.find(";", 0, 1) != std::string::npos){
        NCount = atoi(lineIN.substr(1, lineIN.size()).c_str());
        fIFluxCountRef = (iTableNum == 40) ? NCount : 0;
        restCount = NCount;
        readPara = true;
        if (!isFlux) fMaxFluenceData = (std::max(fMaxFluenceData, (signed)NCount));
        if (isFlux && !fIFluxCountRef) fMaxFluxData = (std::max(fMaxFluxData, (signed)NCount));
        fMaxTestFluxData = fIFluxCountRef;
        //std::cout << "Table->" << iTableNum << " Data-> " << NCount
        //<< " Flux: " << isFlux << " Fluence: " << !(isFlux) << std::endl;
        continue;
      }
      if ( file0 && readPara){
        std::stringstream ss (lineIN);
        ss >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10;
        fExptEnergyBin.push_back(v1); fExptEnergyBin.push_back(v2); fExptEnergyBin.push_back(v3); fExptEnergyBin.push_back(v4);
        fExptEnergyBin.push_back(v5); fExptEnergyBin.push_back(v6); fExptEnergyBin.push_back(v7); fExptEnergyBin.push_back(v8);
        fExptEnergyBin.push_back(v9); fExptEnergyBin.push_back(v10);
        //G4cout << "filecount->" << fileCount << "   " << fExptEnergyBin.size() << G4endl;
        //for (unsigned ijk = 0 ; ijk < fExptEnergyBin.size(); ijk++) std::cout << fExptEnergyBin[ijk] << "  ";
        //std::cout << std::endl;
        file0 = 0;
        readPara = false;
      }
      if (!file0 && readPara){
        std::istringstream sdummy;
        std::stringstream ss;
        sdummy.str(lineIN); // sdummy is used for number of components which deletes sdumy. So ss is required
        ss.str(lineIN);
        v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = v9 = 0.0;
        G4int wcount = std::distance(std::istream_iterator<std::string>(sdummy), std::istream_iterator<std::string>());
        if (wcount == 3) {
          ss >> v1 >> v2  >> v3;
          tmpV1.push_back(v1);  tmpV2.push_back(v2);  tmpV3.push_back(v3);
          --restCount;
        } else if (wcount == 6) {
          ss >> v1 >> v2 >> v3 >> v4 >> v5 >> v6;
          tmpV1.push_back(v1);  tmpV2.push_back(v2);  tmpV3.push_back(v3);
          tmpV1.push_back(v4);  tmpV2.push_back(v5);  tmpV3.push_back(v6);
          restCount -= 2;
        }else if (wcount == 9){
          ss >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9;
          tmpV1.push_back(v1);  tmpV2.push_back(v2);  tmpV3.push_back(v3);
          tmpV1.push_back(v4);  tmpV2.push_back(v5);  tmpV3.push_back(v6);
          tmpV1.push_back(v7);  tmpV2.push_back(v8);  tmpV3.push_back(v9);
          restCount -= 3;
        }
      }
      if (!isFlux && (tmpV1.size()) == NCount){
        if (iTableNum!=0) ++fMaxFluenceTable;
        fExptRadiiTables.push_back(tmpV1);
        fExptFluenceTables.push_back(tmpV2);
        fExptErrTables.push_back(tmpV3);
        std::vector<G4double>().swap(tmpV1);
        std::vector<G4double>().swap(tmpV2);
        std::vector<G4double>().swap(tmpV3);
      } else if (isFlux && (tmpV1.size()) == NCount){
        fExptEnergyTables.push_back(tmpV1);
        fExptFluxTables.push_back(tmpV2);
        fExptFluxErrTables.push_back(tmpV3);
        std::vector<G4double>().swap(tmpV1);
        std::vector<G4double>().swap(tmpV2);
        std::vector<G4double>().swap(tmpV3);
      }
    }
    lineIN="";
  }
  exptIN.close();
  std::vector<G4double>().swap(tmpV1);
  std::vector<G4double>().swap(tmpV2);
  std::vector<G4double>().swap(tmpV3);


  fMaxRadCount = fExptRadiiTables[10].size();
  //  fMaxRadCount = 10  fMaxTestFluxData = 21;
  //  fMaxFluxData  = 95 fMaxFluenceData  = 102

  for (G4int ii = 0; ii < fMaxRadCount; ii++){
    fFluxRadTables.push_back(fExptRadiiTables[10][ii]);
  }


  for (std::size_t ijk2 = 0; ijk2 < fExptRadiiTables[3].size(); ijk2++){
    G4double rVal = fExptRadiiTables[3][ijk2] * 10.0; // in mm
    fRadList.push_back(rVal);
  }

  for (std::size_t ijk1 = 0; ijk1 < fFluxRadTables.size(); ijk1++){
    fFluxRadTables[ijk1] *= 10.0;
  }

  for (std::size_t ij1 = 0; ij1 < fExptEnergyTables[0].size(); ij1++){
    fFlux_Low_Energy.push_back(fExptEnergyTables[0][ij1]);
    fFlux_Low_Data.push_back(fExptFluxTables[0][ij1]);
    fFlux_Low_Syst_Err.push_back(fExptFluxErrTables[0][ij1]);
  }

  for (std::size_t ij1 = 0; ij1 < fExptEnergyTables[1].size(); ij1++){
    fFlux_Lithium_Energy.push_back(fExptEnergyTables[1][ij1]);
    fFlux_Lithium_Data.push_back(fExptFluxTables[1][ij1]);
    fFlux_Lithium_Syst_Err.push_back(fExptFluxErrTables[1][ij1]);
  }

  for (std::size_t ij1 = 0; ij1 < fExptEnergyTables[2].size(); ij1++){
    fFlux_Energy.push_back(fExptEnergyTables[2][ij1]);
    fFlux_Data.push_back(fExptFluxTables[2][ij1]);
    fFlux_Syst_Err.push_back(fExptFluxErrTables[2][ij1]);
  }

  fFlux_Energy_in = fFlux_Energy;
  fFlux_Data_in   = fFlux_Data;
  fFlux_Syst_Err_in = fFlux_Syst_Err;
  fFlux_Low_Energy_in = fFlux_Low_Energy;
  fFlux_Lithium_Energy_in = fFlux_Lithium_Energy;
}


void G4TARCHistoManager::FillRadialExperimentalData(){

  for (G4int ij1 = 0; ij1 < 8; ij1++) {  //  fExptRadiiTables.size(); ij1++){  0~ 41 to 7 ~ 48
    for (std::size_t ij2 = 0; ij2 < fExptRadiiTables[ij1].size(); ij2++){    //   fExptRadiiTables[ij1].size(); ij2++){
      fAnalysisManager->FillNtupleDColumn(9, 0, fExptRadiiTables[ij1][ij2]);
      fAnalysisManager->FillNtupleDColumn(9, 1, fExptEnergyBin[ij1]);
      fAnalysisManager->FillNtupleDColumn(9, 2, fExptFluenceTables[ij1][ij2]);
      fAnalysisManager->FillNtupleDColumn(9, 3, fExptErrTables[ij1][ij2]);
      //G4cout << "9 " << fExptRadiiTables[ij1][ij2] << "  " << fExptEnergyBin[ij1] << "   "
      //       << fExptFluenceTables[ij1][ij2] << "   " << fExptErrTables[ij1][ij2] << G4endl;
      fAnalysisManager->AddNtupleRow(9);
    }
  }


  for (std::size_t ij1 = 8; ij1 < 16; ij1++){ // 8 ~ 49 to 16 ~ 57
    G4int ijE = ij1 - 7;
    for (std::size_t ij2 = 0; ij2 < fExptRadiiTables[ij1].size(); ij2++){    //   fExptRadiiTables[ij1].size(); ij2++){
      fAnalysisManager->FillNtupleDColumn(10, 0, fExptRadiiTables[ij1][ij2]);
      fAnalysisManager->FillNtupleDColumn(10, 1, fExptEnergyBin[ijE]);
      fAnalysisManager->FillNtupleDColumn(10, 2, fExptFluenceTables[ij1][ij2]);
      fAnalysisManager->FillNtupleDColumn(10, 3, fExptErrTables[ij1][ij2]);
      //G4cout << "10 " << fExptRadiiTables[ij1][ij2] << "  " << fExptEnergyBin[ij1] << "   "
      //       << fExptFluenceTables[ij1][ij2] << "   " << fExptErrTables[ij1][ij2] << G4endl;
      fAnalysisManager->AddNtupleRow(10);
    }
  }

  G4cout << "Experimental data filling complete." << G4endl;
}

void G4TARCHistoManager::BeginOfRun() {
  fAnalysisManager = G4AnalysisManager::Instance();
  G4String path = getenv("dateStr");
  fAnalysisFileName = path + "/" + fAnalysisFileName;
  fAnalysisManager->OpenFile(fAnalysisFileName);

  if (!fHistoBooked) {
    BookHistogram();
    CreateTuples();
  }

  DefineShellBlocks();
  ReadExperimentalDataFromFile(fExptlDataFileName);
  FillRadialExperimentalData();


  fNtuple_full = false;
  fAbsX0 = fAbsY0 = fAbsZ0 = 0.0; // 0.5 * fLength;
  fNevt       = 0;
  fNelec      = 0;
  fNposit     = 0;
  fNgam       = 0;
  fNstep      = 0;
  fNprot_leak = 0;
  fNmuons     = 0;
  fNions      = 0;
  fNdeut      = 0;
  fNalpha     = 0;
  fNneutron   = 0;
  fNproton    = 0;
  fNaproton   = 0;
  fNneu_forw  = 0;
  fNneu_leak  = 0;
  fNneu_back  = 0;
  fNinelastic = 0;
  fNelastic   = 0;
  fNzero      = 0;
  fNeutron_check = 0;
  fNeutronStack = 0;
  fNumber_newTrack = 0;
  fGamma_flux = fNeutron_flux = fElectron_flux = fPiminus_flux = fPiPlus_flux = fPizero_flux = fPositron_flux = 0;
  fProton_flux = fMuon_flux = fOther_flux = fNEUTRON_fluence= 0;
  fEdepSum    = 0.0;
  fEdepSum2   = 0.0;
  fNstepEnergy = fMaxEVal / (G4double)fMaxEBin;  // max 8 GeV considered
  fVirtualDia = 3300.0*mm;
  fVirtVol    = (4.0 / 3.0) * CLHEP::pi * (0.5 * fVirtualDia) * (0.5 * fVirtualDia) * (0.5 * fVirtualDia);
  fRange      = fMaxLVal;
  fRho        = 5.0 * mm;  //-------------------------->  Diagnose with GDML file  ********
  fLMax       = fMaxLVal;
  fLBin       = fRange / (fLMax / mm);
  fTcut       = 1.0 * CLHEP::MeV;   // 100.0 tried with 300, 100 and result does not change.
  fNeutronInit = fNeutronSum = fNeutronBreed = 0.0;
  fEsecond    = G4DataVector(fMaxBin, 0.0); //  fNumMax, 0.);
  fMsecond    = G4DataVector(fMaxBin, 0.0); //  fNumMax, 0.);
  // n-spectra, E, T and ET
  fTmin       = 0.01 * CLHEP::eV;
  fTmax       = fMaxEVal;   //  ***** I want to give value ~ 2.0 * incident beam energy ******* // CHECK
  // fTmax    = 1000. * CLHEP::MeV;
  fTimeMin   = 1.0 * nanosecond;
  fTimeMax   = 1.0 * second;   //1.0 * millisecond;
  fNbin      = fMaxBin; // fNumMax; // 60; // 100;  // 1000; // 10; //

  fNmax = 0;
  fTotal_flux = 0.0;
  fTARC_Integral = 0.0; fTARC_Integral_E = 0.0; fTARC_lithium = 0.0;
  fTARC_lithium_E = 0.0; fTARC_helium = 0.0; fTARC_helium_E = 0.0;
  fExiting_Flux = 0;
  fExiting_check_Flux = 0;
  fExiting_Energy = 0.0;

  fNEsecond  = G4PhysicsLogVector(fTmin,fTmax,fNbin);
  fNTsecond  = G4PhysicsLogVector(fTimeMin,fTimeMax,fNbin);
  fNSecondSum1  = G4DataVector(fNbin, 0.0);
  fNSecondSum2  = G4DataVector(fNbin, 0.0);
  fNSecondSum3  = G4DataVector(fNbin, 0.0);
  fNETsum       = G4DataVector(fNbin * fNbin, 0.0);
  fLocal_Energy_Integral = G4DataVector(4, 0.0);

  startEnergy = 0.0;
  flag = false;
  number_generations = 0;
  fractional_bin_width = 0.2;

  fFluence_Spectrum.resize(1000, 0.0);
  fFluence.resize(fMaxFluenceTable, std::vector<G4double>(fMaxEBin, 0.0));
  fFluence_Step_Shell.resize(fMaxTestFluxData, 0.0);
  fFluence_Cyl.resize(fMaxTestFluxData, 0.0);
  fFluence_step.resize(fMaxTestFluxData, 0.0);

  fLow_Fluence_step.resize(fMaxFluenceData, 0.0);
  fLow_Fluence_Step_Shell.resize(fMaxFluenceData, 0.0);

  fRadial_Fluence_Step.resize(fMaxTestFluxData, std::vector<G4double>(fMaxRadCount, 0.0));

  fFlux.resize(fMaxTestFluxData, 0.0);
  fFlux_Radius.resize(fMaxRadCount, std::vector<G4double>(fMaxRadCount, 0.0));
  fEFlux.resize(fMaxTestFluxData, 0.0);
  fLow_Flux.resize(fMaxFluxData, 0.0);
  fENflux.resize(4, 0.0);
  fNeutflux.resize(4, 0.0);

  fCos_Lithium_Flux.resize(fMaxFluxData, 0.0);
  fCos_Low_Flux.resize(fMaxFluxData, 0.0);
  fCos_Flux.resize(fMaxFluxData, 0.0);

  fLithium_Flux.resize(fMaxFluxData, 0.0);
  fLithium_Radial_Mean.resize(fMaxRadCount, 0.0);
  fLithium_Radial_Energy_Lower.resize(fMaxRadCount, 0.0);
  fLithium_Radial_Energy_Upper.resize(fMaxRadCount, 0.0);
  fLithium_Fluence_Step.resize(fMaxFluenceData, 0.0);
  fLithium_Fluence_Step_Shell.resize(fMaxFluenceData, 0.0);


  fET.resize(fNbin, std::vector<G4double>(fNbin, 0.0));
  fEdNdE.resize(fNbin, std::vector<G4double>(fNbin, 0.0));


  fGunParticleX   = G4DataVector(fLMax, 0.0);
  fGunParticleY   = G4DataVector(fLMax, 0.0);
  fGunParticleZ   = G4DataVector(fLMax, 0.0);
  fGunParticleTLX = G4DataVector(fLMax, 0.0);
  fGunParticleTLY = G4DataVector(fLMax, 0.0);
  fGunParticleTLZ = G4DataVector(fLMax, 0.0);
  fGunParticleDep = G4DataVector(fLMax, 0.0);
  fGunParticleRho = G4DataVector(fLMax, 0.0);
  fRangeVector    = G4DataVector(fLMax, 0.0);
  fRhoVector      = G4DataVector(fLMax, 0.0);
  fDeltaVector    = G4DataVector(fLMax, 0.0);
  fNEfluxBin      = G4DataVector(fMaxBin, 0.0);
  for( G4int ii = 0; ii < fLMax; ii++ ) {
    fRangeVector[ii] = ii * fLBin - 0.5 * fRange; //
    fRhoVector[ii]   = ii * fRho / fLBin;
    fDeltaVector[ii] = ii * fTmax / fLBin;
  }
  fTkin = fTcut; // fEVal0 // Does it affect?
  fEbin = fTkin / fMaxBin;

  G4double EDelta = (fMaxEVal / fMaxBin);
  G4double MDelta = (fEVal0 / fMaxBin);
  for( G4int ii = 0; ii < fMaxBin; ii++ )
  {
    fEsecond[ii]  = ii * EDelta;
    // fNEsecond[ii] = ii * fEbin;
    fMsecond[ii]  = ii * MDelta;
  }
  if( fVerbose > 0 ) G4cout << "G4TARCHistoManager: Histograms are booked and run has been started:" << G4endl;
}


void G4TARCHistoManager::CreateNeutronFluxHisto(){
  fAbsolute_TotalFlux = (fTotal_flux * 1.0e9 / (G4double)fNevt) / fTestSphereSurfaceArea;

  for (G4int ij1 = 0; ij1 < fMaxTestFluxData; ij1++){
    G4double mean_energy   = 0.5 * (fFlux_Energy[ij1 + 1] + fFlux_Energy[ij1]);
    G4double abs_flux      = (fFlux[ij1] * ( 1.0e9 / (G4double)fNevt)) / fTestSphereSurfaceArea;
    G4double bin_tmp_width = (fFlux_Energy[ij1 + 1] - fFlux_Energy[ij1]);
    G4double abs_flux_perp = mean_energy * (((fCos_Flux[ij1] * (1.0e9 / (G4double)fNevt)) / fTestSphereSurfaceArea) / bin_tmp_width);
    G4double abs_eflux     = (fEFlux[ij1] * ( 1.0e9 / (G4double)fNevt)) /fTestSphereSurfaceArea;
    G4double abs_fluence   = mean_energy * (100.0 * (1.0e9 / (G4double)fNevt) * (fFluence_step[ij1] / fTestSphereVolume) / bin_tmp_width);

    G4double abs_error = 0.0;
    if (fFlux[ij1] != 0.0) abs_error = (std::pow(fFlux[ij1], 0.5) / fFlux[ij1]) * abs_flux;

    fTARC_helium   += fFlux_Data[ij1];
    fTARC_helium_E += fFlux_Data[ij1] * mean_energy;

    fLocal_Energy_Integral[2] += abs_flux * mean_energy;
    fLocal_Energy_Integral[3] += fFlux_Data[ij1] * mean_energy;

    fAnalysisManager->FillNtupleDColumn(3, 0, mean_energy);
    fAnalysisManager->FillNtupleDColumn(3, 1, fFlux_Data[ij1]);
    fAnalysisManager->FillNtupleDColumn(3, 2, fFlux_Syst_Err[ij1]);
    fAnalysisManager->FillNtupleDColumn(3, 3, abs_flux);
    fAnalysisManager->FillNtupleDColumn(3, 4, abs_flux_perp);
    fAnalysisManager->FillNtupleDColumn(3, 5, abs_fluence);
    fAnalysisManager->FillNtupleDColumn(3, 6, abs_error);
    fAnalysisManager->FillNtupleDColumn(3, 7, fFlux[ij1]);
    fAnalysisManager->FillNtupleDColumn(3, 8, 0.0);         // eflux_data
    fAnalysisManager->FillNtupleDColumn(3, 9, abs_eflux);
    fAnalysisManager->FillNtupleDColumn(3, 10, fFluence_step[ij1]);
    fAnalysisManager->FillNtupleDColumn(3, 12, 0.0);           // abs fluence cyl
    /*
    fAnalysisManager->FillNtupleDColumn(3, 13, 0.0);           // abs fluence front
    fAnalysisManager->FillNtupleDColumn(3, 14, 0.0);           // abs fluence shell
    */
    fAnalysisManager->AddNtupleRow(3);
  }

  for (G4int ij1 = 0; ij1 < fMaxFluenceData; ij1++){
    G4double mean_low_energy   = std::exp(0.5 * (std::log10(fFlux_Low_Energy[ij1]) + std::log10(fFlux_Low_Energy[ij1])));
    G4double abs_low_flux      = (fFlux_Low_Data[ij1] * ( 1.0e9 / (G4double)fNevt)) / fTestSphereSurfaceArea;
    G4double bin_tmp_width2    = fFlux_Low_Energy[ij1 + 1] - fFlux_Low_Energy[ij1];
    G4double abs_low_flux_perp = mean_low_energy * (((fCos_Low_Flux[ij1] *
      ( 1.0e9 / (G4double)fNevt)) / fTestSphereSurfaceArea) / bin_tmp_width2);
    G4double abs_low_fluence   = 100.0 * ( 1.0e9 / (G4double)fNevt) * ((fLow_Fluence_step[ij1]) / fTestSphereVolume);


    G4double abs_error = 0.0;
    if (fFlux_Low_Data[ij1] != 0.0) abs_error = (std::pow(fFlux_Low_Data[ij1], 0.5) / fFlux_Low_Data[ij1]) * abs_low_flux;
    fTARC_Integral += fFlux_Low_Data[ij1];
    fTARC_Integral_E += fFlux_Low_Data[ij1] * mean_low_energy;

    fAnalysisManager->FillNtupleDColumn(4, 0, mean_low_energy);
    fAnalysisManager->FillNtupleDColumn(4, 1, fFlux_Low_Data[ij1]);
    fAnalysisManager->FillNtupleDColumn(4, 2, fFlux_Low_Syst_Err[ij1]);
    fAnalysisManager->FillNtupleDColumn(4, 3, abs_low_flux);
    fAnalysisManager->FillNtupleDColumn(4, 4, abs_low_flux_perp);
    fAnalysisManager->FillNtupleDColumn(4, 5, abs_low_fluence);
    fAnalysisManager->FillNtupleDColumn(4, 6, abs_error);
    fAnalysisManager->FillNtupleDColumn(4, 7, fLow_Flux[ij1]);
    /*
    fAnalysisManager->FillNtupleDColumn(4, 9, 0.0);         // low_fluence_step
    fAnalysisManager->FillNtupleDColumn(4, 10, 0.0);        // abs low fluence cyl
    fAnalysisManager->FillNtupleDColumn(4, 11, 0.0);           // abs low fluence front
    fAnalysisManager->FillNtupleDColumn(4, 12, 0.0);           // abs low fluence shell
    */
    fAnalysisManager->AddNtupleRow(4);
  }

  for (G4int ij1 = 0; ij1 < fMaxFluenceData; ij1++){
    G4double mean_lithium_energy   = std::exp(0.5 * (std::log10(fFlux_Lithium_Energy[ij1 + 1]) + std::log10(fFlux_Lithium_Energy[ij1])));
    G4double abs_lithium_flux      = (fFlux_Lithium_Data[ij1] * (1.0e9 / (G4double)fNevt)) / fTestSphereSurfaceArea;
    G4double bin_tmp_width_lithium = fFlux_Low_Energy[ij1 + 1] - fFlux_Low_Energy[ij1];
    G4double abs_lithium_flux_perp = mean_lithium_energy * (((fCos_Lithium_Flux[ij1] * (1.0e9 / (G4double)fNevt))
      / fTestSphereSurfaceArea) / bin_tmp_width_lithium);
    G4double abs_lithium_fluence   = 100.0 * (1.0e9 / (G4double)fNevt) * ((fLithium_Fluence_Step[ij1]) / fTestSphereVolume);
    //G4double abs_lithium_fluence_shell = mean_lithium_energy * (100.0 * (1.0e9 / (G4double)fNevt)
    //  * ((fLithium_Fluence_Step_Shell[ij1]) / fTestShellVol) / bin_tmp_width_lithium);

    G4double abs_error_li = 0.0;
    if (fFlux_Low_Data[ij1] != 0.0) abs_error_li = (std::pow(fFlux_Lithium_Data[ij1], 0.5) / fFlux_Lithium_Data[ij1]) * abs_lithium_flux;

    fTARC_lithium   += fFlux_Lithium_Data[ij1];
    fTARC_lithium_E += fFlux_Lithium_Data[ij1] * mean_lithium_energy;

    fAnalysisManager->FillNtupleDColumn(5, 0, mean_lithium_energy);
    fAnalysisManager->FillNtupleDColumn(5, 1, fFlux_Lithium_Data[ij1]);
    fAnalysisManager->FillNtupleDColumn(5, 2, fFlux_Low_Syst_Err[ij1]);
    fAnalysisManager->FillNtupleDColumn(5, 3, abs_lithium_flux);
    fAnalysisManager->FillNtupleDColumn(5, 4, abs_lithium_flux_perp);
    fAnalysisManager->FillNtupleDColumn(5, 5, abs_lithium_fluence);
    fAnalysisManager->FillNtupleDColumn(5, 6, abs_error_li);
    fAnalysisManager->FillNtupleDColumn(5, 7, fLithium_Flux[ij1]);
    fAnalysisManager->FillNtupleDColumn(5, 8, fLithium_Fluence_Step[ij1]);
    fAnalysisManager->AddNtupleRow(5);
  }
}

void G4TARCHistoManager::CreateRadialFluxHisto(){
  for (G4int ijk1 = 0; ijk1 < fMaxTestFluxData; ijk1++) {  // ijk1 < fMaxTestFluxData; ijk1++){
    G4double shellVol = (4.0 / 3.0) * CLHEP::pi * (std::pow(fOuterRadiusofShell[ijk1], 3.0) - std::pow(fInnerRadiusofShell[ijk1], 3.0));
    G4double radL = 0.5 * (fOuterRadiusofShell[ijk1] + fInnerRadiusofShell[ijk1]);
    for (G4int ijk2 = 0; ijk2 < fMaxRadCount; ijk2++){
      G4double bin_tmp_width = fLithium_Radial_Energy_Upper[ijk2] - fLithium_Radial_Energy_Lower[ijk2];
      if (fRadial_Fluence_Step[ijk1][ijk2] != 0.0){
        G4double fAbs_Rad_Fluence = fLithium_Radial_Mean[ijk2]*(100.0 * (1.0e9/(G4double)fNevt) *
        ((fRadial_Fluence_Step[ijk1][ijk2]) / shellVol) / bin_tmp_width);
        fAnalysisManager->FillNtupleDColumn(8, 0, radL);
        fAnalysisManager->FillNtupleDColumn(8, 1, fLithium_Radial_Mean[ijk2]);
        fAnalysisManager->FillNtupleDColumn(8, 2, fAbs_Rad_Fluence);
        fAnalysisManager->FillNtupleDColumn(8, 3, 0.0);
        fAnalysisManager->FillNtupleDColumn(8, 4, 0.0);
        fAnalysisManager->AddNtupleRow(8);
      }
    }
  }
}

void G4TARCHistoManager::EndOfRun() {
  G4cout << "G4TARCHistoManager ; End of run actions are started" << G4endl;
  G4cout << "fNevt = " << fNevt << G4endl;
  G4cout << "EndOfRun(), fEdepSum = " << fEdepSum << G4endl;
  G4cout << "======================================================================" << G4endl;

  StartProcessing();
  CreateNeutronFluxHisto();
  CreateRadialFluxHisto();

  G4double x = ( G4double )fNevt;
  G4double perN = 1.0;
  if (fNevt > 0){
    x = perN = 1.0 / x;
  }
  TrackRun(x);  // track and get leaks
  NeutronRun(x); // neutron leak data and spectra
  GunParticleRun(x);  // beam distribution in target.
  fAnalysisManager->Write();
  fAnalysisManager->CloseFile();
}

void G4TARCHistoManager::BeginOfEvent(G4int nEvt) {
  fEdepEvt    = 0.0;
  fEdepEM     = 0.0;
  fEdepProton = 0.0;
  fEdepPI     = 0.0;
  fNsecondary = 0.0;
  fRangeSum   = G4ThreeVector(0.0, 0.0, 0.0);
  fStepSum    = 0.0;
  fDeltaSum   = 0.0;
  SetEventID(nEvt);
}

void G4TARCHistoManager::EndOfEvent() {
  fEdepSum  += fEdepEvt;
  fEdepSum2 += fEdepEvt * fEdepEvt;
  fNevt ++;
  if (fNsecondary > 0) fNinelastic++;
  //  fHisto->Fill(21,fEdepEvt/fPrimaryKineticEnergy,1.0);
  //  fHisto->Fill(22,fEdepEM/fPrimaryKineticEnergy,1.0);
  //  fHisto->Fill(23,fEdepPI/fPrimaryKineticEnergy,1.0);
  //  fHisto->Fill(24,fEdepProton/fPrimaryKineticEnergy,1.0);

  fAnalysisManager->FillH1(7, fNeutronStack);
  WriteEventRange( fRangeSum, fStepSum, fDeltaSum );
}


// Normalization for no hadron interaction in the target. EM and hadron elastic shuld be inactivated
void G4TARCHistoManager::AddNzero(const G4Track* myTrack) {
  G4double cosTheta = myTrack->GetDynamicParticle()->GetMomentumDirection().x();
  // change x, y z
  if (cosTheta > 0.999) {
    fNzero++;
  }else {
    fNelastic++;
  }
}


void G4TARCHistoManager::GunParticleDistribution( const G4Step* myStep) {
  G4double length = myStep->GetStepLength();
  G4double dEdx   = myStep->GetTotalEnergyDeposit();

  G4ThreeVector pGun1 = myStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector pGun2 = myStep->GetPostStepPoint()->GetPosition();

  fRangeSum += (pGun2 - pGun1);
  fStepSum += length;
  fDeltaSum += dEdx;
}

void G4TARCHistoManager::AddTargetStep(const G4Step* myStep) {
  fNstep++;

  G4double fEdep = myStep->GetTotalEnergyDeposit();

  if (fVerbose > 1) {
    G4cout << "TargetSD::ProcessHits: beta1 =" << myStep->GetPreStepPoint()->GetVelocity() / c_light
           << " beta2 = "                       << myStep->GetPostStepPoint()->GetVelocity() / c_light
           << " weight = "                      <<  myStep->GetTrack()->GetWeight()
           << G4endl;
  }
  //if (fEdep >= DBL_MIN) {
  const G4Track* myTrack = myStep->GetTrack();
  G4ThreeVector pos = 0.5 * ( myStep->GetPreStepPoint()->GetPosition()
                          + myStep->GetPostStepPoint()->GetPosition()
                          );

  //G4double x = pos.x(); //- fAbsX0;
  //G4double y = pos.y(); //- fAbsY0;
  G4double z = pos.z();// - fAbsZ0;

  //Scoring.
  fEdepEvt += fEdep;

  //G4cout << " Energy for histogram-> " << fEdep/keV << " keV" << G4endl;
  //  fHisto->Fill(0, z, fEdep/MeV);
  const G4ParticleDefinition* pd = myTrack->GetDefinition();

  if (pd == G4Gamma::Gamma() || pd == G4Electron::Electron() || pd == G4Positron::Positron()) {
    fEdepEM += fEdep;
  } else if (pd == G4Proton::Proton() || pd == G4AntiProton::AntiProton()){
    fEdepProton += fEdep;
  } else if (pd == G4PionPlus::PionPlus() || pd == G4PionMinus::PionMinus()){
    fEdepPI += fEdep;
  }

  if (fVerbose > 1) G4cout << "HistoManager::AddEnergy: E(keV) = " << fEdep / keV
                           << "; z(mm) = " << z / mm
                           << "; step(mm) = " << myStep->GetStepLength() / mm
                           << " by " << pd->GetParticleName()
                           << " E(MeV) = " << myTrack->GetKineticEnergy() / MeV
                           << G4endl;
  //}
}


void G4TARCHistoManager::ScoreNewTrack( const G4Track* myTrack) {
  fNumber_newTrack++;
  const G4ParticleDefinition* pd = myTrack->GetDefinition();
  G4String name = pd->GetParticleName();
  G4double ke = myTrack->GetKineticEnergy();
  //G4double TKin = myTrack->GetKineticEnergy();
  //if (myTrack->GetTrackID() == 1)
  //  G4cout << name << "  energy = " << ke << "  ID---> " << myTrack->GetTrackID() << G4endl;

// For Primary Track
  if (myTrack->GetTrackID() == 1) {
    fNevt++;
    if (ke/MeV > 0){
      G4double enerMean = GetGPSEnergy();
      enerMean = (0.5 * (ke + enerMean));
      //SetGPSEnergyIN(enerMean);  // already set in primaryGenerator from GPS
    }
    //fPrimaryKineticEnergy = ke;
    fPrimaryDef = pd;
    G4ThreeVector dir = myTrack->GetMomentumDirection();

    if (fVerbose > 1) G4cout << "### Primary "   << name << " KineticE(MeV) = " << ke / MeV
                             << "; Mass(MeV) = " << pd->GetPDGMass() / MeV
                             << "; Pos(mm) = "   << myTrack->GetPosition() / mm
                             << "; Dir = "       << myTrack->GetMomentumDirection()
                             << G4endl;
  } else {
// For Secondary tracks.
    if (fVerbose > 1) G4cout << "=== Secondary " << name << " KineticE(MeV) = " << ke / MeV
                             << "; Mass(MeV) = " << pd->GetPDGMass() / MeV
                             << "; Pos(mm) = "   << myTrack->GetPosition() / mm
                             << "; Dir = "       << myTrack->GetMomentumDirection()
                             << G4endl;
    ke = std::log10(ke / eV);

    if (pd == G4Gamma::GammaDefinition()){
      fNgam++;
      //const G4VProcess* procG = myTrack->GetCreatorProcess();         //  check
      //  fHisto->Fill(1, ke, 1.0);
    }else if (pd == G4Electron::ElectronDefinition()){
      fNelec++;
      // // fHisto->Fill(2, ke, 1.0);
    }else if (pd == G4Positron::PositronDefinition()) {
      fNposit++;
      //  fHisto->Fill(3, ke, 1.0);
    } else if (pd == G4Proton::ProtonDefinition()){
      fNproton++;
      // // fHisto->Fill(4, ke, 1.0);
    } else if (pd == G4Neutron::NeutronDefinition()){//&& TKin < fTcut){  // <----- CHECK
      fNneutron++;
      fNeutronStack++;
      // fEventAction->AddNeutronStack();
      //fHisto->Fill(5, ke, 1.0);
      //G4cout << ke << G4endl;
      fAnalysisManager->FillNtupleDColumn(15, 0, log10(ke));
      fAnalysisManager->AddNtupleRow(15);
    } else if (pd == G4AntiProton::AntiProtonDefinition()){
      fNaproton++;
    } else if ( pd == G4PionPlus::PionPlusDefinition() ) {
      fNpions++;
      // fHisto->Fill(6, ke, 1.0);
      //  fHisto->Fill(19, ke, 1.0);
    } else if ( pd == G4PionMinus::PionMinusDefinition()) {
      fNpions++;
      //  fHisto->Fill(6, ke, 1.0);
      //  fHisto->Fill(20, ke, 1.0);
    } else if ( pd == G4PionZero::PionZeroDefinition()) {
      fNpi0++;
      //  fHisto->Fill(7, ke, 1.0);
    } else if ( pd == G4KaonPlus::KaonPlusDefinition() ||  pd == G4KaonMinus::KaonMinusDefinition()) {
      fNkaons++;
      //  fHisto->Fill(8, ke, 1.0);
    } else if ( pd == G4KaonZeroShort::KaonZeroShortDefinition() ||  pd == G4KaonZeroLong::KaonZeroLongDefinition()) {
      fNkaons++;
      //  fHisto->Fill(9, ke, 1.0);
    } else if (pd == G4Deuteron::DeuteronDefinition() || pd == G4Triton::TritonDefinition()){
      fNdeut++;
      // // fHisto->Fill(10, ke, 1.0);
    } else if (pd == G4He3::He3() || pd == G4Alpha::Alpha()) {
      fNalpha++;
      //  fHisto->Fill(11, ke, 1.0);
    } else if (pd->GetParticleType() == "nucleus"){
      fNions++;
      //  fHisto->Fill(12, ke, 1.0);
    } else if (pd == G4MuonPlus::MuonPlusDefinition() || pd == G4MuonMinus::MuonMinusDefinition()){
      fNmuons++;
      //  fHisto->Fill(13, ke, 1.0);
    }
  }
}

void G4TARCHistoManager::AddLeakingParticle(const G4Track* myTrack) {
  const G4ParticleDefinition* pd           = myTrack->GetDefinition();
  G4String pname  = pd->GetParticleName();
  //const G4DynamicParticle*    dynSecondary = myTrack->GetDynamicParticle();
  //G4double                    Tkin         = dynSecondary->GetKineticEnergy();

  if (myTrack->GetKineticEnergy()/MeV <= 0.0) {
    //G4double mass   = pd->GetPDGMass();
    return;
  }

  //G4double en       = std::log10(myTrack->GetKineticEnergy()/eV);
  G4ThreeVector pos = myTrack->GetPosition();
  G4ThreeVector dir = myTrack->GetMomentumDirection();

  G4double xx       = pos.x();
  G4double yy       = pos.y();
  G4double zz       = pos.z();
  G4bool isLeaking  = false;


  if ( (
         // (xx > -fAbsX0 && dir.x() > 0.0) // - was added later originally > +fAbsX0
         // ||   (yy > -fAbsY0 && dir.y() > 0.0)
         // ||
         (zz > fAbsZ0 && dir.z() > 0.0)     //  (zz > -fHLength && dir.z() > 0.0)
       )
  ) { // forward
      isLeaking = true;
      if (pname == "neutron") { // && Tkin < fTcut) {
        ++fNneu_forw;
        //  fHisto->Fill(15, en,  1.0);
      }
  }
  if ( (
         // (xx < fAbsX0 && dir.x() < 0.0)
         // ||   (yy < fAbsY0 && dir.y() < 0.0)
         // ||
           (zz < -fAbsZ0 && dir.z() < 0.0) //(zz < fHLength && dir.z() < 0.0)
       )
  ) { // backward
    isLeaking = true;
    if (pname == "neutron") { // && Tkin < fTcut)   {
       ++fNneu_back;
       //  fHisto->Fill(16, en,  1.0);
     }
  }
  if ( (
         // (std::abs(xx) <= fAbsX0 && (zz * dir.z()  + yy * dir.y()) > 0.0 )
         // ||   (std::abs(yy) <= fAbsY0 && (xx * dir.x()  + zz * dir.z()) > 0.0 )
         // ||
         //(std::abs(zz) <= -fHLength && (xx * dir.x()  + yy * dir.y()) > 0.0 )
         (std::abs(xx) <= fAbsX0 || std::abs(yy) <= fAbsY0)
         && ((xx * dir.x()  + yy * dir.y()) > 0.0 )
       )
  ) { // side
    isLeaking = true;
    if (pd == fNeutron){ // && Tkin < fTcut)  {
      ++fNneu_leak;
      // fHisto->Fill(14, en,  1.0);
    }
  }

  if (isLeaking) {
    if (pd == G4Proton::Proton()) { // && myTrack->GetTrackID() == 1){
       ++fNprot_leak;
       // fHisto->Fill(17, en,  1.0);
     }
    if (pd == G4PionPlus::PionPlus() || pd == G4PionMinus::PionMinus()) {
      ++fNPionleak;
      // fHisto->Fill(18, en,  1.0);
    }
  }
}


// pA->nX low energy data neutron final state
void G4TARCHistoManager::NeutFinalState(const G4Track* myTrack){   //}, const G4Step* myStep) {
  fNsecondary++;
  G4double TKin = 0.0; //, Tspectra = 0.0, stepLen, time;
  //G4int ix;
  G4int ii, jj;

  if (myTrack->GetDynamicParticle()->GetDefinition()->GetParticleName() == "neutron"){
	if ((myTrack->GetLogicalVolumeAtVertex()->GetName() == "blockB_log")
       || (myTrack->GetLogicalVolumeAtVertex()->GetName() == "blockA_log")) fNeutronSum+= 1.0;

    // future plan for neutron breed at breeder with X
    TKin = myTrack->GetDynamicParticle()->GetKineticEnergy();
    //stepLen = myStep->GetStepLength();
    //time = myTrack->GetGlobalTime();
    for (ii = 0; ii < fNbin; ii++) {
      if (TKin <= fNEsecond.GetLowEdgeEnergy(ii)){
        jj = (ii == 0) ? ii : ii - 1;
        break;
      }
    }
    jj = (ii == fNbin) ? fNbin-1 : jj;
  }
}


void G4TARCHistoManager::TargetProfile(const G4Track* myTrack){ //}, const G4Step* myStep) {
  //G4double TKinE = 0.0, TSpectr = 0.0, stepLen, time, xn;
  G4double rn;
  G4int ix = 0;

 if (myTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "neutron") {
    if (
	   ((myTrack->GetLogicalVolumeAtVertex()->GetName() == "blockB_log") || (myTrack->GetLogicalVolumeAtVertex()->GetName() == "blockA_log"))
			&& (myTrack->GetParentID() == 1)
			) {
			fNeutronInit += 1.0;
		}
	}

  if (myTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "neutron") {
    if (
            ((myTrack->GetLogicalVolumeAtVertex()->GetName() == "blockB_log")
            || (myTrack->GetLogicalVolumeAtVertex()->GetName() == "blockA_log"))
           //  && (  (myTrack->GetVertexPosition().x() >= -35.0 && myTrack->GetVertexPosition().x() <= 35.0)
           //  || (myTrack->GetVertexPosition().y() >= -35.0 && myTrack->GetVertexPosition().y() <= 35.0)
           //  || (myTrack->GetVertexPosition().z() >= -1500.0 && myTrack->GetVertexPosition().z() <= -299.0)
           // )
           && 
		   (myTrack->GetGlobalTime() <= 10.0 * nanosecond)
           && (myTrack->GetParentID() == 1)
     ) {
       rn = myTrack->GetVertexPosition().x() + 0.5 * fRange;
       ix = G4int(rn / fLBin + 0.5);
       //  if (ix >= 0 && ix < fLMax) fGunParticleX[ix] += 1.0;
       if (ix >= 0 && ix < fVirtualDia) fGunParticleX[ix] += 1.0;

       rn = myTrack->GetVertexPosition().y() + 0.5 * fRange;
       ix = G4int(rn / fLBin + 0.5);
       //  if (ix >= 0 && ix < fLMax) fGunParticleY[ix] += 1.0;
       if (ix >= 0 && ix < fVirtualDia) fGunParticleY[ix] += 1.0;

       rn = myTrack->GetVertexPosition().z() + 0.5 * fRange;
       ix = G4int(rn / fLBin + 0.5);
       // if (ix >= 0 && ix < fLMax) fGunParticleZ[ix] += 1.0;
       if (ix >= 0 && ix < fVirtualDia) fGunParticleZ[ix] += 1.0;
     }
  }
  // if (fNeutronInit) G4cout << "N init " << fNeutronInit << G4endl;
}

void G4TARCHistoManager::AddEnergyTimeHole(const G4Track* myTrack, const G4Step* myStep) {
  G4double KE, myTime;//, myStepLength;
  //size_t ii, jj, kk;

  std::vector<G4double> tmp;

  if (myTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "neutron"){
    KE = myTrack->GetDynamicParticle()->GetKineticEnergy();
    //myStepLength = myStep->GetStepLength();
    myTime = myTrack->GetGlobalTime();

    G4TouchableHandle touch1 = myStep->GetPreStepPoint()->GetTouchableHandle();
    std::string LVname1 = touch1->GetVolume()->GetLogicalVolume()->GetName();
    std::size_t pos1 = LVname1.find("_");
    G4int ijk1= std::atoi(LVname1.substr(4, pos1 - 4).c_str());
    //G4TouchableHandle touch2 = myStep->GetPostStepPoint()->GetTouchableHandle();
    //std::string LVname2 = touch2->GetVolume()->GetLogicalVolume()->GetName();
    //std::size_t pos2 = LVname2.find("_");
    //G4int ijk2= std::atoi(LVname2.substr(4, pos2 - 4).c_str());

    G4double eval = KE / eV;
    G4int fluxEBin = eval / (fNstepEnergy);
    G4cout << "stepE " << fNstepEnergy << G4endl;
    fluxEBin = (fluxEBin >= fMaxEBin) ? fMaxEBin - 1 : fluxEBin;

    // insert LV, time, energy to fETVirtual

    tmp.push_back(G4double(ijk1));
    tmp.push_back(myTime / microsecond);
    tmp.push_back(eval);

    fETVirtual.push_back(tmp); // needs sorting on ijk1 first and then on myTime
    std::sort(fETVirtual.begin(), fETVirtual.end(),
      [](const std::vector< G4double >& a, const std::vector< G4double >& b){
        if (a[0] == b[0])
          return a[1] < b[1];
        else
          return a[0] < b[0];
        } );

    //G4cout << ijk1 << " E= " << eval << " Ebin= " << fluxEBin << " step= " << myStepLength << G4endl;

    std::vector<G4double>().swap(tmp);
    //fFluence[ijk1][fluxEBin] += myStepLength;
    //fFluence[ijk1][fluxEBin] /= fTotVolVBox;
  }
}

void G4TARCHistoManager::AddEnergyTime(const G4Track* myTrack){ //, const G4Step* myStep) {
  G4double Tkin = 0.0, myTime, myGlobalTime;   //, myStepLength;
  G4double Qenergy, Qtime; // , eVal;
  G4int ii, jj, eii = 0, tjj = 0;
  std::vector<G4double> tmp;

  Tkin = myTrack->GetDynamicParticle()->GetKineticEnergy();
  //eVal = Tkin / eV;
//  myStepLength = myStep->GetStepLength();

  if (myTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "neutron") {
    fTrackID = myTrack->GetTrackID();
    myTime = myTrack->GetProperTime();
    myGlobalTime = myTrack->GetGlobalTime();
    for (ii = 0; ii < fNbin; ii++) {
      if (Tkin <= fNEsecond.GetLowEdgeEnergy(ii)) {
        eii = (ii == 0) ? ii : ii - 1;
        //Qenergy = (Tkin/eV != 0.0) ? (Tkin/eV) : 0.0;  // Tkin; //
        Qenergy = Tkin / eV;
        break;
      }
    }
    if (ii == fNbin) eii = fNbin -1 ;
    for (jj = 0; jj < fNbin; jj++) {
      if (myTime <= fNTsecond.GetLowEdgeEnergy(jj)) {
        Qtime = (myGlobalTime/microsecond);
        tjj = (jj == 0) ?  jj : jj - 1;
        break;
      }
    }
    //G4cout << "jj : " << jj << "  tjj: " << tjj << "  eii: " << eii << G4endl;
    if (jj == fNbin) tjj = fNbin - 1;
    fET[tjj][eii] += 1.0;

    tmp.push_back(Qtime);
    tmp.push_back(Qenergy);
    fNSpectra.push_back(tmp);
    std::vector<G4double>().swap(tmp);
    fOldTrackID = fTrackID;
  }

}

void G4TARCHistoManager::WriteEventRange(G4ThreeVector rsum, G4double lsum, G4double dsum) {
  G4int ix;
  G4double x, y, z, rho, xbin, ybin, zbin, rbin, dbin;    //, dEdx, length;

  x = rsum.x();
  y = rsum.y();
  z = rsum.z();

  rho = std::sqrt(x * x + y * y);

  xbin = ybin = zbin = fRange / fMaxBin;   // 100 can be changed checking the binning
  rbin               = fRho / fMaxBin;
  dbin               = fTmax  / fMaxBin;

  ix = (G4int)(x / (xbin) + 0.5);
  if (ix > 0 && ix < fMaxBin) fGunParticleX[ix] += 1.0;

  ix = (G4int)(y / (ybin) + 0.5);
  if (ix > 0 && ix < fMaxBin) fGunParticleY[ix] += 1.0;

  ix = (G4int)(z / (zbin) + 0.5);
  if (ix > 0 && ix < fMaxBin) fGunParticleZ[ix] += 1.0;

  ix = (G4int)(lsum / xbin + 0.5);
  if (ix > 0 && ix < fMaxBin) fGunParticleTLX[ix] += 1.0;

  ix = (G4int)(lsum / ybin + 0.5);
  if (ix > 0 && ix < fMaxBin) fGunParticleTLY[ix] += 1.0;

  ix = (G4int)(lsum / zbin + 0.5);
  if (ix > 0 && ix < fMaxBin) fGunParticleTLZ[ix] += 1.0;

  ix = (G4int)(dsum / dbin + 0.5);
  if (ix > 0 && ix < fMaxBin) fGunParticleDep[ix] += 1.0;

  ix = (G4int)(rho / rbin + 0.5);
  if (ix > 0 && ix < fMaxBin) fGunParticleRho[ix] += 1.0;
}


template <typename T>
void G4TARCHistoManager::Check10s(T inVal, T &outVal, G4String &uniStr){
  G4int iresult = 0;
  div_t divresult = div((inVal / MeV), 1000);
  if (divresult.quot >= 1){
    iresult = 9;  // unit = Giga
  } else {
    divresult = div((inVal / keV), 1000);
    if (divresult.quot >=1){
      iresult = 6;   // Mega
    } else {
      divresult = div((inVal / eV), 1000);
      if (divresult.quot >= 1) {
        iresult = 3;   // kilo
      } else {
        iresult = 1;
      }
    }
  }
  switch(iresult){
    case 9: outVal = inVal / GeV; uniStr = " GeV"; break;
    case 6: outVal = inVal / MeV; uniStr = " MeV"; break;
    case 3: outVal = inVal / keV; uniStr = " keV"; break;
    case 1: outVal = inVal; uniStr = " eV"; break;
  }
}

void G4TARCHistoManager::TrackRun(G4double x) {
  std::ofstream trackout("track.dat", std::ios::out);
  G4double xStep         = x * (G4double)fNstep;
  G4double xElectron     = x * (G4double)fNelec;
  G4double xPositron     = x * (G4double)fNposit;
  G4double xGamma        = x * (G4double)fNgam;
  G4double xProton       = x * (G4double)fNproton;
  G4double xAntiProton   = x * (G4double)fNaproton;
  G4double xProtonLeak   = x * (G4double)fNprot_leak;
  G4double xNeutron      = x * (G4double)fNneutron;
  G4double xNeutronLeak  = x * (G4double)fNneu_leak;
  G4double xneuF         = x * (G4double)fNneu_forw;
  G4double xneuB         = x * (G4double)fNneu_back;
  G4double xMuons        = x * (G4double)fNmuons;
  G4double xIons         = x * (G4double)fNions;
  G4double xAlpha        = x * (G4double)fNalpha;
  G4double xDeut         = x * (G4double)fNdeut;
  G4double xp0           = x * (G4double)fNPionleak;

  fEdepSum  *= x;
  fEdepSum2 *= x;
  fEdepSum2 -= fEdepSum * fEdepSum;
  fEdepSum2  = (fEdepSum2 > 0.0) ? std::sqrt(fEdepSum2) : 0.0;

  trackout << "x = 1/Number_of_Events = "  << x << " for number of Events = " << fNevt   << G4endl;

  //div_t divresult = div((fEdepSum/MeV),1000);
  G4double trVal = 0.0;
  G4String enerUnit = "";
  Check10s(fEdepSum, trVal, enerUnit);
  trackout << "  Per event data:" << G4endl;
  trackout << std::setprecision(4) << "Energy Deposited <mean> " << trVal << enerUnit;
  Check10s(fEdepSum2, trVal, enerUnit);
  trackout << "  <rms> " << trVal << enerUnit <<  G4endl;
  Check10s(fPrimaryKineticEnergy, trVal, enerUnit);
  trackout << "Beam energy = " << trVal << enerUnit << G4endl;
  trackout << "                                                                        "  << G4endl;

  trackout << " Production in Target (per event):"        << "                                     "  << G4endl;
  trackout << std::setprecision(4) << "Average Number of steps: "       << xStep             << G4endl;
  trackout << std::setprecision(4) << "Average Number of Gammas: "      << xGamma            << G4endl;
  trackout << std::setprecision(4) << "Average Number of Electrons: "   << xElectron         << G4endl;
  trackout << std::setprecision(4) << "Average Number of Positrons: "   << xPositron         << G4endl;
  trackout << std::setprecision(4) << "Average Number of Protons: "     << xProton           << G4endl;
  trackout << std::setprecision(4) << "Average Number of AntiProton: "  << xAntiProton       << G4endl;
  trackout << std::setprecision(4) << "Average Number of Neutrons: "    << xNeutron          << G4endl;
  trackout << std::setprecision(4) << "Average Number of Muons: "       << xMuons            << G4endl;
  trackout << std::setprecision(4) << "Average Number of D + T: "       << xDeut             << G4endl;
  trackout << std::setprecision(4) << "Average Number of He3 + alpha: " << xAlpha            << G4endl;
  trackout << std::setprecision(4) << "Average Number of ions: "        << xIons             << G4endl;
  trackout << "                                                      "  << "              "  << G4endl;

  trackout << " Leakage from the system (per event): "                                            << G4endl;
  trackout << std::setprecision(4) << "Average Number of forward Neutrons: "      << xneuF        << G4endl;
  trackout << std::setprecision(4) << "Average Number of reflected Neutrons: "    << xneuB        << G4endl;
  //trackout << std::setprecision(4) << "Average Number of other leaked Neutrons " << xNeutronLeak << G4endl;
  trackout << std::setprecision(4) << "Average Number of total leaked Neutrons: " <<  xNeutronLeak
                                                                                    + xneuF
                                                                                    + xneuB      << G4endl;
  trackout << std::setprecision(4) << "Average Number of leaked Protons: "        << xProtonLeak  << G4endl;
  trackout << std::setprecision(4) << "Average Number of leaked Pions: "          << xp0          << G4endl;
  trackout <<                                                                                     G4endl;


  trackout << "==========================================================" << G4endl;
  trackout << "Exiting flux : " << fExiting_Flux << G4endl;
  Check10s(fExiting_Energy, trVal, enerUnit);
  trackout << "Total Exiting Energy : " << trVal << enerUnit << G4endl;
  Check10s(fAnalysisManager->GetH1(1)->mean(), trVal, enerUnit);
  trackout << " Gamma Edep <mean> : " << trVal << enerUnit;
  Check10s(fAnalysisManager->GetH1(1)->rms(), trVal, enerUnit);
   trackout << " rms: " << trVal << enerUnit << G4endl;
  trackout << "Neutron Lethargy <mean> : " << fAnalysisManager->GetH1(2)->mean() << " rms: " << fAnalysisManager->GetH1(2)->rms()  << G4endl;
  trackout << "==========================================================" << G4endl;

  G4double kEffective, rho, rat, react, perN=x;
  kEffective = (fNeutronInit!= 0.0) ? fNeutronSum / fNeutronInit : 0.0;
  rho        =  (kEffective != 0.0) ? (kEffective - 1.0) / kEffective : 0.0;  // reactivity :: deviation from criticality
  rat        = (kEffective != 0.0) ? std::log(kEffective) : 0.0;
  react      = rat / (1.0 + rat);

  trackout << " IMP Parameters : "    << G4endl;
  trackout << " Neutron_Init/p = "    << fNeutronInit* perN << ",  Neutron_Sum/p = " << fNeutronSum * perN << G4endl;
  trackout << " kEffective = "        << kEffective         << ", Rho = "            << rho                << G4endl;
  trackout << " Estimated reactivity = " << react                                    << G4endl             << G4endl;
  trackout << "==========================================================================================" << G4endl;
  trackout << G4endl;

}


void G4TARCHistoManager::NeutronRun(G4double x) {
  G4double perN = x;

  std::ofstream fETVirt("ETVirtual.dat", std::ios::out);
  //fETVirt << (fETVirtual.size() * fETVirtual[0].size()) << G4endl;
  for (std::size_t ii = 0; ii < fETVirtual.size(); ii++){
    for (std::size_t jj = 0; jj < fETVirtual[0].size(); jj++){
      fETVirt << fETVirtual[ii][jj] << "  ";
    }
    fETVirt<< G4endl;
  }
  fETVirt.close();

/*
  std::ofstream fFlu("FluenceData.dat", std::ios::out);
  std::cout << fFluence.size() << G4endl;
  for (std::size_t ii = 0; ii < fFluence.size(); ii++){
    for (std::size_t jj = 0; jj < fMaxEBin; jj++){
      fFlu << ii << "  " << jj << " " << fFluence[ii][jj] << G4endl;
    }
  }
  fFlu.close();
*/

  G4double ngSum = 0.0;
  std::ofstream nspec("neutronSpectra.dat", std::ios::out);
  nspec << fNSecondSum1.size() << G4endl;
  for (size_t k = 0; k < fNSecondSum1.size(); k++){
      nspec << fNEsecond.GetLowEdgeEnergy(k)/MeV << "  "  << perN * fNSecondSum1[k]
                                                 << "   " << perN * fNSecondSum2[k]
                                                 << "  "  << perN * fNSecondSum3[k] << G4endl;
      ngSum += fNSecondSum1[k];
  }
  G4cout << "integral S-spectrum per event = " << ngSum * perN << G4endl;

//--------------------------------------------------
  std::ofstream tenspectr("tenspectr.dat",std::ios::out);
  tenspectr<<fNETsum.size()<<G4endl;
  for( G4int k = 0; k < fNbin; k++ ){
    for( G4int j = 0; j < fNbin; j++ ){
      // tenspectr<<perN*fNETsum[k]<<G4endl;
      if (perN * fET[j][k] != 0.0)
        tenspectr << (G4int)k << "    " << (G4int)j << "    " << perN * fET[j][k] << G4endl;
    }
  }
//----------------------------------------------------
//------------> This is working but commented as it is taking time to process and write huge data at the end of the execution.
/*
  std::ofstream neutSpec("neutSpec.dat", std::ios::out);
  for (size_t ii = 0; ii < fNSpectra.size(); ii ++){
    neutSpec << fNSpectra[ii][0] << "   " << fNSpectra[ii][1] << G4endl;
  }
  neutSpec.close();
  */
//-------------------------------------------------
  std::ofstream teaxis("teaxis.dat",std::ios::out);
  teaxis << fNEsecond.GetVectorLength() << G4endl;
  for( size_t k = 0; k < fNEsecond.GetVectorLength(); k++ )
  {
    teaxis << fNEsecond.GetLowEdgeEnergy(k) << "        " << fNTsecond.GetLowEdgeEnergy(k) << G4endl;
  }

}


void G4TARCHistoManager::GunParticleRun(G4double x) {
  std::ofstream gunspectrum("gun.dat", std::ios::out);

  gunspectrum << fRangeVector.size() << G4endl;

  for (size_t k = 0; k < fRangeVector.size(); k++) {
    gunspectrum << fRangeVector[k] << "  " << x * fGunParticleX[k]
                                   << "  " << x * fGunParticleY[k]
                                   << "  " << x * fGunParticleZ[k]
                                   << "  " << x * fGunParticleTLX[k]
                                   << "  " << x * fGunParticleTLY[k]
                                   << "  " << x * fGunParticleTLZ[k]
                                   << "  " <<     fDeltaVector[k]
                                   << "  " << x * fGunParticleDep[k]
                                   << "  " <<     fRhoVector[k]
                                   << "  " << x * fGunParticleRho[k] << G4endl;
  }
}


void G4TARCHistoManager::StartProcessing(){
  std::size_t ArrayDim = fFlux_Energy.size();
  G4double mean_energy = 0.0;
  G4int j = 0;
  G4cout << "Entered in StartProcessing. " << G4endl;
  for (std::size_t ii1 = 0; ii1 < ArrayDim; ii1++){
    mean_energy = 0.5 * (fFlux_Energy[ii1] + fFlux_Energy_in[ii1 + 1]);
    fFine_Energy.push_back(fFlux_Energy_in[ii1]);
    fEflux_Data.push_back(mean_energy * fFlux_Data_in[ii1]);
    fEflux_Integral += fFlux_Data_in[ii1] * mean_energy/1.0e6;
    ++j;
    fFine_Energy[j] = mean_energy;
  }

  //for (G4int k = 0; k < ArrayDim; k++){G4cout << " Flux Energy: " << fFine_Energy[k] << G4endl;  }

  G4double bin_width = (std::log(1.0e5) - std::log(0.01)) / 100.0;
  //G4double energy_start = 0.01;

  //fFlux_Low_Energy[0] = energy_start;        //
  //fFlux_Lithium_Energy[0] = energy_start;    //  C H E CK   as data has been read from the file

  G4int radial_index = 0;

  for(G4int ii1 = 0; ii1 < fMaxFluenceData; ii1++){
    fFlux_Low_Energy[ii1 + 1] = std::exp(bin_width + std::log(fFlux_Low_Energy[ii1]));
    fFlux_Lithium_Energy[ii1 + 1] = std::exp(bin_width + std::log(fFlux_Lithium_Energy[ii1]));
    mean_energy = std::exp(0.5 * (std::log(fFlux_Low_Energy[ii1 + 1]) + std::log(fFlux_Low_Energy[ii1])));
    //G4double Lithium_Mean_Energy = std::exp(0.5 * (std::log(fFlux_Lithium_Energy[ii1 + 1]) + std::log(fFlux_Lithium_Energy[ii1])));
    //G4cout << " Lithium Energy: " << fFlux_Lithium_Energy[ii1] << "  radila index: " << radial_index << G4endl;
    //G4cout << " Energy (rindex): " << fExptEnergyBin[radial_index] << G4endl;

    if (fFlux_Lithium_Energy[ii1] < fExptEnergyBin[radial_index]
      && fFlux_Lithium_Energy[ii1 + 1] > fExptEnergyBin[radial_index]){
        fLithium_Radial_Energy_Lower.push_back(fFlux_Lithium_Energy[ii1]);
        fLithium_Radial_Energy_Upper.push_back(fFlux_Lithium_Energy[ii1 + 1]);
        radial_index++;
    }

    //G4double bug_ratio = mean_energy / (fFlux_Low_Energy[ii1 + 1] - fFlux_Low_Energy[ii1]);
    //G4cout << "bug ratio: " << bug_ratio << "  mean_energy: " << mean_energy
    //       << " lowE1: " << fFlux_Low_Energy[ii1] << " lowE2: " << fFlux_Low_Energy[ii1] << G4endl;
  }
}

void G4TARCHistoManager::NeutronEnergyTime(G4double thisE, G4double thisT, G4double E0){
  //G4cout << fNtuple_full << " NEUT:----> thisE : " << thisE << "  thisT : " << thisT << G4endl;
  G4double tempT = thisT / microsecond;
  G4double tempE = thisE / eV;
  G4double tempE0 = E0 / eV;
  if (tempT > 0.0 && tempE > 0.0) fAnalysisManager->FillH2(1, log10(tempT), log10(tempE), 1.0);
  //if (fNtuple_full) {
    fAnalysisManager->FillNtupleDColumn(1, 0, tempE);
    fAnalysisManager->FillNtupleDColumn(1, 1, tempT);
    fAnalysisManager->FillNtupleDColumn(1, 2, tempE0);
    fAnalysisManager->AddNtupleRow(1);
  //}
}

void G4TARCHistoManager::otherEnergyTime(G4double thisE, G4double thisT, G4double E0){
  //G4cout << fNtuple_full << "  OTHERS-------> thisE : " << thisE << "  thisT : " << thisT << G4endl;
  G4double tempT = thisT / microsecond;
  G4double tempE = thisE / eV;
  G4double tempE0 = E0 / eV;
  if (tempT > 0.0 && tempE > 0.0) fAnalysisManager->FillH2(2, log10(tempT), log10(tempE), 1.0);
  //if (fNtuple_full) {
    fAnalysisManager->FillNtupleDColumn(14, 0, tempE);
    fAnalysisManager->FillNtupleDColumn(14, 1, tempT);
    fAnalysisManager->FillNtupleDColumn(14, 2, tempE0);
    fAnalysisManager->AddNtupleRow(14);
  //}
}

void G4TARCHistoManager::exitingTally(G4bool exiting_flag, G4double energyL){
  //G4cout << exiting_flag << "  energyL:  " << energyL << " exitingF: " << fExiting_Flux << "  exitingE: " << fExiting_Energy << G4endl;
  if(exiting_flag) {
    //G4TARCHistoManager::AddExitingFlux(energyL);
    AddExitingFlux(energyL);
    fAnalysisManager->FillNtupleDColumn(2, 0, energyL / MeV);
    fAnalysisManager->AddNtupleRow(2);
  }
}

void G4TARCHistoManager::analysePS(G4double partEnergy, G4String particleName
  //, G4double partTime
  , G4double partMomentum
  //, G4double zMomentum
){
  //G4int Iparticle = -9;
  if(particleName == "gamma") {
    fAnalysisManager->FillH1(1, partEnergy / keV);
  // hGammaEdep->fill(energy/keV);
    //Iparticle = 1;
  }else if(particleName == "neutron") {
  //  testMax1 = std::max(testMax1, partEnergy); testMin1 = std::min(testMin1, partEnergy);
  //  testMax2 = std::max(testMax2, 1.0 / partMomentum); testMin2 = std::min(testMin2, 1.0/partMomentum);

    //G4cout << " filling neutron lethargy with: energy = " << partEnergy / eV << " and momentum = " << 1.0 / partMomentum << G4endl;
    //fAnalysisManager->FillH1(2, partEnergy / eV, 1.0 / partMomentum);
    fAnalysisManager->FillH1(2, partEnergy / eV, 1.0 / partMomentum);
    // hNeutronEdep->fill(energy/eV,1/momentum);  // fill(x,weight)
    //Iparticle = 2;
    if(partEnergy / MeV < 2.0) {
      fNeutflux[0]++;
      fENflux[0] += partEnergy;
    } else if(partEnergy / MeV > 2.0 && partEnergy / MeV < 20.0) {
      fNeutflux[1]++;
      fENflux[1] += partEnergy;
    } else if(partEnergy / MeV > 20.0) {
      fNeutflux[2]++;
      fENflux[2] += partEnergy;
    }
    if(partEnergy / MeV > 1000.0) {
      fNeutflux[3]++;
      fENflux[3] += partEnergy;
    }
  }
  if(particleName == "e-") {
    fAnalysisManager->FillH1(3, partEnergy / keV);
    // hElectronEdep->fill(energy/keV);  // fill(x,weight)
    //Iparticle = 3;
  }
  if(particleName == "e+") {
    fAnalysisManager->FillH1(4, partEnergy / keV);
    // hPositronEdep->fill(energy/keV);  // fill(x,weight)
    //Iparticle = 4;
  }
  if(particleName == "other") {
    fAnalysisManager->FillH1(5, partEnergy / keV);
    // hOtherEdep->fill(energy/keV);  // fill(x,weight)
    //Iparticle = 5;
  }
}


void G4TARCHistoManager::ProcessStepping(const G4Step* myStep){
  //StartProcessing();
  G4int StepNo = myStep->GetTrack()->GetCurrentStepNumber();
  G4double partEnergy = myStep->GetPreStepPoint()->GetKineticEnergy();
  G4double partTime   = myStep->GetPreStepPoint()->GetGlobalTime();
  G4double partMomentum = myStep->GetPreStepPoint()->GetMomentum().mag();
  //G4double zMomentum = myStep->GetPreStepPoint()->GetMomentum().z();
  G4double angle = myStep->GetPreStepPoint()->GetMomentum().angle(myStep->GetPreStepPoint()->GetPosition());
  G4double cosAngle = std::abs(cos(angle));
  G4ParticleDefinition* particleType = myStep->GetTrack()->GetDefinition();
  G4String particleName = particleType->GetParticleName();
  G4double primEnergy = 0.0;

  analysePS(partEnergy, particleName, partMomentum);      //  , partTime, partMomentum, zMomentum);
  if (StepNo == 1 && particleName == "neutron"){
    startEnergy = partEnergy;
    startTime = partTime;
    flag = true;
  }

  G4int thisTrackID = myStep->GetTrack()->GetTrackID();
  G4int parentTrackID = myStep->GetTrack()->GetParentID();

  if (StepNo == 1 && thisTrackID == 1){
    parent_energy.clear();
    parent_particle.clear();
    parent_particleID.clear();
    number_generations = 0;
  }

  parent_energy[thisTrackID] = partEnergy;
  parent_particle[thisTrackID] = particleName;
  parent_particleID[thisTrackID] = parentTrackID;

  G4bool reduced_tally = false;

  if (thisTrackID == 1 && StepNo == 1) primEnergy = partEnergy;
  if (thisTrackID != 1 && StepNo == 1) {
    G4int tempID = thisTrackID;
    while(parent_particleID[tempID] != 1){
      tempID = parent_particleID[tempID];
      ++number_generations;
    }
    if (parent_particle[parentTrackID] == "neutron" && particleName == "neutron"){
      reduced_tally = true;
      parent_particle.erase(parentTrackID);
    }
    analyseSecondaries (partEnergy, particleName, partTime, partMomentum, parentTrackID, primEnergy,
      parent_energy[parentTrackID], parent_particle[parentTrackID], reduced_tally, number_generations);
  }

  if (particleName == "neutron"){
    NeutronEnergyTime(partEnergy, partTime, startEnergy);
  } else {
    if (particleName == "Pb207" || particleName == "Pb208")  otherEnergyTime(partEnergy, partTime, startEnergy);
  }
  G4double radiusPre = myStep->GetPreStepPoint()->GetPosition().mag();
  G4double radiusPost = myStep->GetPostStepPoint()->GetPosition().mag();
  //G4double zPos = myStep->GetPreStepPoint()->GetPosition().z();
  G4double StepLength = myStep->GetStepLength();
  G4String vol = myStep->GetTrack()->GetVolume()->GetName();


  G4TouchableHistory* thePreTouchable  = (G4TouchableHistory*) (myStep->GetPreStepPoint()->GetTouchable());
  G4TouchableHistory* thePostTouchable = (G4TouchableHistory*) (myStep->GetPostStepPoint()->GetTouchable());

  if (particleName == "neutron"){
    G4String PreVolGrichine = thePreTouchable->GetVolume()->GetName();
    if (myStep->GetTrack()->GetNextVolume()){
      G4String PreVol = thePreTouchable->GetVolume()->GetName();
      G4String PostVol = thePostTouchable->GetVolume()->GetName();
      if (PreVol == "lab_phys" && PostVol == "world_log_PV"){
        exitingTally(true, partEnergy);
      }
      if (PostVol == "world_log_PV"){
        exitingTallyCheck(true);
      }
    }
  }

  G4bool pre_inside = false;
  G4bool post_inside = false;
  //G4double radius = -999999.9;
  if  ((radiusPre <= (fMaxOuterRadiusofShell + fMyTol))
   &&  (radiusPre >= (fMinInnerRadiusofShell - fMyTol)) && particleName == "neutron") pre_inside = true;
  if  ((radiusPost <= (fMaxOuterRadiusofShell + fMyTol))
   &&  (radiusPost >= (fMinInnerRadiusofShell - fMyTol)) && particleName == "neutron") post_inside = true;
  if (pre_inside && post_inside){
    analyseNeutronShellFluence(partEnergy, StepLength);
  }

  for (G4int ishell = 0; ishell < fShellNumber; ishell++){
    G4bool pre_inside_radial = false;
    G4bool post_inside_radial = false;
    G4double radOut = fOuterRadiusofShell[ishell];
    G4double radIn  = fInnerRadiusofShell[ishell];
    if ((radiusPre <= (radOut + fMyRadTol))
     && (radiusPre >= (radIn - fMyRadTol)) && (particleName == "neutron")) pre_inside_radial = true;
    if ((radiusPost <= (radOut + fMyRadTol))
     && (radiusPost >= (radIn - fMyRadTol)) && (particleName == "neutron")) post_inside_radial = true;
    if (pre_inside_radial && post_inside_radial){
      analyseNeutronRadialFluence(partEnergy, StepLength, ishell);  // partTime, StepLength, ishell);
    }
  }

  for (std::size_t ii = 0; ii < fFluxRadTables.size(); ++ii){
    if ( (radiusPre < fFluxRadTables[ii] && radiusPost > fFluxRadTables[ii])
       ||(radiusPre > fFluxRadTables[ii] && radiusPost < fFluxRadTables[ii])){
         G4double radiusL = fFluxRadTables[ii];

         analyseNeutronFlux(partEnergy, thisTrackID, radiusL, cosAngle, particleName);
    }
  }
}

void G4TARCHistoManager::analyseSecondaries(G4double energyL, G4String nameL, G4double timeL, G4double momentumL,
  G4int ParentIDL, G4double primaryEnergyL, G4double parentEnergyL, G4String parentParticleL, G4bool reduced_fluxL,
  G4int number_generationsL){
  G4int Iparticle=-9;
  G4double temp_time = timeL / microsecond;
  G4double temp_energy = energyL / eV;
  G4double temp_momentum = momentumL;
  if(nameL == "gamma") {
    fGamma_flux++;
    Iparticle = 1;
  } else if(nameL == "neutron") {
    if(!reduced_fluxL) fNeutron_check++;
    fNeutron_flux++;
    Iparticle = 2;
    if(reduced_fluxL) Iparticle = -2;
    if(energyL > 0.1*eV && energyL < 10.0*keV) fNEUTRON_fluence++;
  } else if(nameL == "e-") {
    fElectron_flux++;
    Iparticle = 3;
    //    return;
  } else if(nameL == "pi-") {
    fPiminus_flux++;
    Iparticle = 4;
  } else if(nameL == "pi+") {
    fPiPlus_flux++;
    Iparticle = 5;
  } else if(nameL == "pi0") {
    fPizero_flux++;
    Iparticle = 6;
  } else if(nameL == "e+") {
    fPositron_flux++;
    Iparticle = 7;
  } else if(nameL == "proton") {
    fProton_flux++;
    Iparticle = 8;
  } else if(nameL == "mu-") {
    fMuon_flux++;
    Iparticle = 9;
  } else if(nameL == "mu+") {
    fMuon_flux++;
    Iparticle = 10;
  } else {
    fOther_flux++;
    Iparticle = 99;
    return;
  }
  G4int iParent = 0;
  if (parentParticleL == "gamma")        iParent = 1;
  else if (parentParticleL == "neutron") iParent = 2;
  else if(reduced_fluxL)                 iParent = -2;
  else if (parentParticleL == "e-")      iParent = 3;
  else if (parentParticleL == "pi-")     iParent = 4;
  else if (parentParticleL == "pi+")     iParent = 5;
  else if (parentParticleL == "pi0")     iParent = 6;
  else if (parentParticleL == "e+")      iParent = 7;
  else if (parentParticleL == "proton")  iParent = 8;
  else if (parentParticleL == "proton")  iParent = 8;
  else if (parentParticleL == "mu-")     iParent = 9;
  else if (parentParticleL == "mu+")     iParent = 10;


  //if(fNtuple_full) {
    fAnalysisManager->FillNtupleDColumn(0,0, temp_energy);
    fAnalysisManager->FillNtupleDColumn(0,1, temp_time);
    fAnalysisManager->FillNtupleIColumn(0,2, Iparticle);
    fAnalysisManager->FillNtupleDColumn(0,3, temp_momentum);
    fAnalysisManager->FillNtupleIColumn(0,4, ParentIDL);
    fAnalysisManager->FillNtupleDColumn(0,5, primaryEnergyL);
    fAnalysisManager->FillNtupleIColumn(0,6, iParent);
    fAnalysisManager->FillNtupleDColumn(0,7, parentEnergyL);
    fAnalysisManager->FillNtupleIColumn(0,8, number_generationsL);
    fAnalysisManager->FillNtupleIColumn(0,9, fNevent_id);
    fAnalysisManager->AddNtupleRow(0);
  //}
}


void G4TARCHistoManager::analyseNeutronRadialFluence(G4double partEnergyL, //G4double partTimeL,
G4double StepLengthL, G4int ishellL){
  //G4cout << "radial index : " << ishellL << " fShellNumber: " << fShellNumber << G4endl;
  if (ishellL < 0 || ishellL > fShellNumber) G4cout << "WARNING! radial index is wrong !!!!!!!" << G4endl;
  G4double tempEnergy = partEnergyL / eV;
  if (tempEnergy < fLithium_Radial_Energy_Upper[9] && tempEnergy > fLithium_Radial_Energy_Lower[0]){
    for (G4int i = 0 ; i < fMaxRadCount; i++) {
      if (tempEnergy > fLithium_Radial_Energy_Lower[i] && tempEnergy > fLithium_Radial_Energy_Upper[i]){
        fRadial_Fluence_Step[ishellL][i] += StepLengthL / mm;
      }
    }
  }
}

//void G4TARCHistoManager::analyseNeutronFlux(G4double n_EnergyL, G4double partTimeL, G4double startEnergyL,
//  G4int thisTrackIDL, G4int parentTrackIDL, G4double zMomentumL, G4double startTimeL, G4double radiusL,
//  G4double zPosL, G4double parent_energyL, G4String parent_particleL, G4double cosAngleL,
//  G4int number_generationsL, G4String particleNameL, G4bool reduced_tallyL)
void G4TARCHistoManager::analyseNeutronFlux(G4double n_EnergyL, G4int thisTrackIDL, G4double radiusL,
  G4double cosAngleL, G4String particleNameL)
  //G4double zPosL,G4double cosAngleL, G4String particleNameL)
  {
    if (particleNameL == "neutron"){
      if (thisTrackIDL == fOldTrackID && std::abs(radiusL - 456.0*mm) < 0.1){
        ++fDuplicate_neutrons;
      } else {
        fDuplicate_neutrons = 0;
      }
      fOldTrackID = thisTrackIDL;
    }
    G4double tempEnergy = n_EnergyL / eV;
    if (particleNameL == "neutron"){
      for (G4int ii1 = 0; ii1 < fMaxRadCount; ii1++){
        if (std::abs(radiusL - fFluxRadTables[ii1]) < 0.1){
          for (G4int ii2 = 0; ii2 < fMaxRadCount; ii2++){
            if (std::abs(n_EnergyL - fExptEnergyBin[ii2]) < fractional_bin_width * fExptEnergyBin[ii2])
              fFlux_Radius[ii1][ii2] += 1.0 / std::abs(cosAngleL);
          }
        }
      }
      if (std::abs(radiusL - 456.0*mm) < 0.1){
        if (n_EnergyL > 0.345*eV && n_EnergyL < 1.0e5*eV){
          fTARC_Integral   += 1.0 / std::abs(cosAngleL);
          fTARC_Integral_E += n_EnergyL / eV * 1.0 / std::abs(cosAngleL);
        }
        G4int nVal = (G4int)((2.0 + std::log10(n_EnergyL / eV)) / 0.09);
        if (nVal < 0) nVal = 0;
        if (nVal < 100) {
          fFluence_Spectrum[nVal] += 1.0 /std::abs(cosAngleL);
          if (nVal > fNmax) fNmax = nVal;
        }

        if (tempEnergy > 0.0194 && tempEnergy < 1.0e5){
          fTARC_lithium   += 1.0 / std::abs(cosAngleL);
          fTARC_lithium_E += tempEnergy * 1.0 / std::abs(cosAngleL);
        }
        if (tempEnergy > 59500.0 && tempEnergy < 1825092.0){
          fTARC_helium   += 1.0 / std::abs(cosAngleL);
          fTARC_helium_E += tempEnergy * 1.0 / std::abs(cosAngleL);
        }
      }

      //G4cout << "radiusL----------> " << radiusL << G4endl;
      std::size_t LithiumMax = fFlux_Lithium_Energy.size();
      //G4cout << "LithiumMax:-> " << LithiumMax << "   FluxLithium: -> " << fFlux_Lithium.size() << G4endl;
      if (std::abs(radiusL - 456.0 * mm) < 0.1){
        fTotal_flux++;
        if (n_EnergyL / eV < fFlux_Energy[0]) fLocal_Energy_Integral[0] += n_EnergyL / eV;
        if (tempEnergy < fFlux_Lithium_Energy[LithiumMax-1]){
          for (std::size_t ijk1 = 0; ijk1 < LithiumMax; ijk1++){
            if(tempEnergy > fFlux_Lithium_Energy[ijk1] && tempEnergy < fFlux_Lithium_Energy[ijk1 + 1]){
              fLithium_Flux[ijk1]++;
              fCos_Lithium_Flux[ijk1] = 1.0 / std::abs(cosAngleL);
            }
          }
        }
        if (tempEnergy < fFlux_Low_Energy[fMaxFluxData]){
          for (G4int ijk1 = 0; ijk1 < fMaxFluxData; ijk1++){
            if (tempEnergy >fFlux_Low_Energy[ijk1] && tempEnergy < fFlux_Low_Energy[ijk1 + 1]){
              fLow_Flux[ijk1]++;
              fCos_Low_Flux[ijk1] += 1.0 / std::abs(cosAngleL);
            }
          }
        } else if ( tempEnergy > fFlux_Energy[0]) {
            for(G4int ijk1 = 0; ijk1 < fMaxTestFluxData; ijk1++){
              if (tempEnergy > fFlux_Energy[ijk1] && tempEnergy < fFlux_Energy[ijk1 + 1]){
                fLocal_Energy_Integral[1] += n_EnergyL / eV;
                fFlux[ijk1]++;
                if (cosAngleL != 0.0) fCos_Flux[ijk1] += 1.0 / std::abs(cosAngleL);
                fEFlux[ijk1] += n_EnergyL / eV;
              }
            }
        }
      }
    }

}


void G4TARCHistoManager::analyseNeutronShellFluence(G4double energyL, G4double StepLengthL){
    G4double tempE = energyL / eV;
    std::size_t lastTagFluence = fFlux_Lithium_Energy.size();
    //G4cout << "lastElementFluence : " << lastTagFluence << "  analyseNeutronShellFluence in Histo" << G4endl;
    if (tempE < fFlux_Lithium_Energy[lastTagFluence - 1]){
      for (std::size_t ii1 = 0; ii1 < lastTagFluence; ii1++){
        if (tempE > fFlux_Lithium_Energy[ii1] && tempE < fFlux_Lithium_Energy[ii1 + 1]) fLithium_Fluence_Step_Shell[ii1] += StepLengthL / mm;
      }
    }
    std::size_t lastTagLowEFlux = fFlux_Low_Energy.size();
    //G4cout << "lastElementLowEFlux : " << lastTagLowEFlux << "  analyseNeutronShellFluence in Histo" << G4endl;
    if (tempE < fFlux_Low_Energy[lastTagLowEFlux - 1]){
      for (std::size_t ii1 = 0; ii1 < lastTagLowEFlux; ii1++){
        if (tempE > fFlux_Low_Energy[ii1] && tempE < fFlux_Low_Energy[ii1 + 1]) fLow_Fluence_Step_Shell[ii1] += StepLengthL / mm;
      }
    }
    std::size_t lastTagFluxE = fFlux_Energy.size();
    //G4cout << "lastElementFluxE : " << lastTagFluxE << "  analyseNeutronShellFluence in Histo" << G4endl;
    if (tempE > fFlux_Energy[0]){
      for (std::size_t ii1 = 0; ii1 < lastTagFluxE; ii1++){
        if (tempE > fFlux_Energy[ii1] && tempE < fFlux_Energy[ii1 + 1]) fFluence_Step_Shell[ii1] += StepLengthL / mm;
      }
    }
  }


//void G4TARCHistoManager::Fill(G4int id, G4double x, G4double w) {
  // fHisto->Fill(id, x, w);
//}


  // // // // // // // // // // // // // // // // //
/*
  void G4TARCHistoManager::SetVerbose(G4int val)
  {
    fVerbose = val;
    fHisto->SetVerbose(val);
  }
*/
