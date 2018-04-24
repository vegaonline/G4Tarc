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
:
  fHisto( 0 ),
  fHistoBooked(false),
  fPrimaryDef( 0 ),
  fProton( 0 ),
  //fNHisto( 25 ),
  fNeutron( 0 ),
  fEdepMax( fMaxEVal ),
  fLength( fMaxLVal ),
  fHLength(0.5 * fMaxLVal),
  fPrimaryKineticEnergy ( fEVal0 ),
  fVerbose( 0 ),
  fNBinsE( fMaxBin ),
  fProtonIN(0),
  fNCountTotal(0),
  fVirtVol(0),
  fNSlices( fMaxSlices ) {
    fHisto    = new G4TARCHisto();
    fHisto->SetVerbose(fVerbose);
    //fNeutron = G4Neutron::Neutron();
    //fProton = G4Proton::Proton();
}

G4TARCHistoManager::~G4TARCHistoManager() {
  // delete fHisto;
  std::vector<G4double>().swap(fExptEnergyBin);
  std::vector<std::vector<G4double> >().swap(fExptRadiiTables);
  std::vector<std::vector<G4double> >().swap(fExptFluenceTables);
  std::vector<std::vector<G4double> >().swap(fExptErrTables);
  std::vector<std::vector<G4double> >().swap(fExptEnergyTables);
  std::vector<std::vector<G4double> >().swap(fExptFluxTables);
  std::vector<std::vector<G4double> >().swap(fExptFluxErrTables);
  std::vector<G4double>().swap(fOuterRadiusofShell);
  std::vector<G4double>().swap(fInnerRadiusofShell);

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
  fShellThickness        =     2.0 * mm;
  fMinInnerRadiusofShell =    10.0 * mm;
  fMaxOuterRadiusofShell =  1500.0 * mm;
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
  fTestSphereVolume = (4.0 / 3.0) * pi * (fTestSphereRadius * fTestSphereRadius * fTestSphereRadius);
  fTestSphereSurfaceArea = 4.0 * pi * (fTestSphereRadius * fTestSphereRadius);
  fRadHole = 32.0 * mm;
  fLenCyl  = 150.0 * mm;
}


void G4TARCHistoManager::BookHistogram() {
  fHistoBooked = true;
  fAnalysisManager->SetFirstHistoId(1);
  fAnalysisManager->CreateH1("Protons_per_Event", "Protons/event", 50, 0.0, 100.0);
  fAnalysisManager->CreateH1("Neutrons_per_Event", "Neutrons/event", 50, 0.0, 100.0);
  fAnalysisManager->CreateH1("Neutron_Energy", "Neutron Energy vs 1/mom /eV", 100000, 0.0, 1000000.0);
  fAnalysisManager->CreateH1("proton_Energy", "Proton Energy Deposition/keV", 1000, 0.0, 1000.0 * keV);
  fAnalysisManager->CreateH1("Particle_Stack", "Particle Stack", 12, 0.5, 12.5);
  fAnalysisManager->CreateH1("log10En", "Log10 Energy (eV) of neutron", 100, -9.0, 9.0);
  fAnalysisManager->CreateH2("Neutron_Energy_Time", "Neutron Energy vs Time", 100, -1.0, 4.0, 100, -2.0, 7.0);
}


void G4TARCHistoManager::CreateTuples(){
  fAnalysisManager->CreateNtuple("G4TARC_Secondaries", "Secondary Particle Info");
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

  fAnalysisManager->CreateNtuple("G4TARC_Energy_Time", "Neutron Time");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("time");
  fAnalysisManager->CreateNtupleDColumn("primary");
  fAnalysisManager->FinishNtuple(); // ntupleID: 1

  fAnalysisManager->CreateNtuple("G4TARC_Exiting", "Neutrons Exiting");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->FinishNtuple(); // ntupleID: 2 - filled

  fAnalysisManager->CreateNtuple("G4TARC_Flux_4002", "Neutrons G4TARC flux");
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

  fAnalysisManager->CreateNtuple("G4TARC_Flux_4004", "Neutrons G4TARC flux");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("tarcflux");
  fAnalysisManager->CreateNtupleDColumn("errstat");
  fAnalysisManager->CreateNtupleDColumn("g4flux");
  fAnalysisManager->CreateNtupleDColumn("g4perp");
  fAnalysisManager->CreateNtupleDColumn("gfluence");
  fAnalysisManager->CreateNtupleDColumn("g4err");
  fAnalysisManager->CreateNtupleDColumn("rawflux");
  fAnalysisManager->CreateNtupleDColumn("gstep");
  fAnalysisManager->CreateNtupleDColumn("gfl_cyl");
  fAnalysisManager->CreateNtupleDColumn("g4front");
  // fAnalysisManager->CreateNtupleDColumn("g4_shell");
  fAnalysisManager->FinishNtuple(); // ntupleID: 4

  fAnalysisManager->CreateNtuple("G4TARC_Created_Neutrons", "Created Neutrons");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("time");
  fAnalysisManager->CreateNtupleIColumn("particle");
  fAnalysisManager->CreateNtupleDColumn("momentum");
  fAnalysisManager->CreateNtupleDColumn("zmom");
  fAnalysisManager->FinishNtuple(); // ntupleID: 5

  fAnalysisManager->CreateNtuple("G4TARC_Created_Neutrons", "Created Neutrons");
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
  fAnalysisManager->FinishNtuple(); // ntupleID: 6

  fAnalysisManager->CreateNtuple("G4TARC_Radial_Shell_Fluence", "Radial Shell Fluence");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("fluence");
  fAnalysisManager->CreateNtupleDColumn("true_e");
  fAnalysisManager->CreateNtupleDColumn("true_f");
  fAnalysisManager->FinishNtuple(); // ntupleID: 7

  //fAnalysisManager->CreateNtuple("G4TARC_Radial_Fluence_Exptl_Data", "Radial Fluence Expt. Data");
  fAnalysisManager->CreateNtuple("h8", "Radial Fluence Expt. Data");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("data");
  fAnalysisManager->CreateNtupleDColumn("error");
  fAnalysisManager->FinishNtuple(); // ntupleID: 8

  fAnalysisManager->CreateNtuple("G4TARC_Radial_Fluence_He3", "Radial Fluence He3");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("data");
  fAnalysisManager->CreateNtupleDColumn("error");
  fAnalysisManager->FinishNtuple(); // ntupleID: 9


  fAnalysisManager->CreateNtuple("G4TARC_3.5GeV_He3_Expt_Data", "Radial Fluence He3 Expt Data");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("data");
  fAnalysisManager->CreateNtupleDColumn("stat_err");
  fAnalysisManager->CreateNtupleDColumn("syst_err");
  fAnalysisManager->FinishNtuple(); // ntupleID: 10

  fAnalysisManager->CreateNtuple("G4TARC_Radial_Fluence_Li", "Radial Fluence Li");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("data");
  fAnalysisManager->CreateNtupleDColumn("stat_err");
  fAnalysisManager->CreateNtupleDColumn("syst_err");
  fAnalysisManager->FinishNtuple(); // ntupleID: 11

  fAnalysisManager->CreateNtuple("G4TARC_Radial_Fluence", "Radial Fluence");
  fAnalysisManager->CreateNtupleDColumn("radius");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("fluence");
  fAnalysisManager->CreateNtupleDColumn("he_data");
  fAnalysisManager->FinishNtuple(); // ntupleID: 12

  fAnalysisManager->CreateNtuple("G4TARC_Energy_Time", "OTHER Time");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->CreateNtupleDColumn("time");
  fAnalysisManager->CreateNtupleDColumn("primary");
  fAnalysisManager->FinishNtuple(); // ntupleID: 13

  fAnalysisManager->CreateNtuple("log10(En)", "log10(En)");
  fAnalysisManager->CreateNtupleDColumn("energy");
  fAnalysisManager->FinishNtuple(); // ntupleID: 14 for neutron

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
  std::cout << " TARC Experimental data file reading report:-" << std::endl;
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
    } else {   // the line does not start with # symbol
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
        //for (unsigned ijk = 0 ; ijk < fExptEnergyBin.size(); ijk++) std::cout << fExptEnergyBin[ijk] << "  ";
        //std::cout << std::endl;
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
        std::vector<G4double>().swap(tmpV2);
        std::vector<G4double>().swap(tmpV1);
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
}

void G4TARCHistoManager::FillRadialExperimentalData(){
  for (G4int jj = 0; jj < 8; jj++){   // this is the row of 2D vector for read Data
    //G4cout << jj << "   " << fExptRadiiTables[jj].size() << G4endl;
    for (unsigned ii = 0; ii < fExptRadiiTables[jj].size(); ii++){
      fAnalysisManager->FillNtupleDColumn(8, 0, fExptRadiiTables[jj][ii]);
      fAnalysisManager->FillNtupleDColumn(8, 1, fExptEnergyBin[jj]);
      fAnalysisManager->FillNtupleDColumn(8, 2, fExptFluenceTables[jj][ii]);
      fAnalysisManager->FillNtupleDColumn(8, 3, fExptErrTables[jj][ii]);
      fAnalysisManager->AddNtupleRow(8);
    }
  }
  for (G4int jj = 8; jj < 17; jj++){
    //G4cout << jj << "   " << fExptRadiiTables[jj].size() << G4endl;
    for (unsigned ii = 0; ii < fExptRadiiTables[jj].size(); ii++){
      fAnalysisManager->FillNtupleDColumn(9, 0, fExptRadiiTables[jj][ii]);
      fAnalysisManager->FillNtupleDColumn(9, 1, fExptEnergyBin[jj]);
      fAnalysisManager->FillNtupleDColumn(9, 2, fExptFluenceTables[jj][ii]);
      fAnalysisManager->FillNtupleDColumn(9, 3, fExptErrTables[jj][ii]);
      fAnalysisManager->AddNtupleRow(9);
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
  G4cout << " File is being read from Data file" << G4endl;
  ReadExperimentalDataFromFile(fExptlDataFileName);

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

  fTARC_Integral = 0.0; fTARC_Integral_E = 0.0; fTARC_lithium = 0.0; fTARC_lithium_E = 0.0; fTARC_helium = 0.0; fTARC_helium_E = 0.0;

  fNEsecond  = G4PhysicsLogVector(fTmin,fTmax,fNbin);
  fNTsecond  = G4PhysicsLogVector(fTimeMin,fTimeMax,fNbin);
  fNSecondSum1  = G4DataVector(fNbin, 0.0);
  fNSecondSum2  = G4DataVector(fNbin, 0.0);
  fNSecondSum3  = G4DataVector(fNbin, 0.0);
  fNETsum       = G4DataVector(fNbin * fNbin, 0.0);


  for( G4int ii = 0; ii < fNbin; ii++ )  {
    std::vector<G4double> temp;
    for( G4int jj = 0; jj < fNbin; jj++ ) temp.push_back(0.0);
    fET.push_back(temp);
    fEdNdE.push_back(temp);
  }

  for (G4int ii=0; ii < fMaxFluenceTable; ii ++) {
    std::vector<G4double> temp;
    for (G4int jj = 0; jj < fMaxEBin; jj++) temp.push_back(0.0);
    fFluence.push_back(temp);
  }

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
  if( fVerbose > 0 )
    G4cout << "G4TARCHistoManager: Histograms are booked and run has been started"
           << G4endl;
}


void G4TARCHistoManager::CreateNeutronFluxHisto(){
  fAbsolute_TotalFlux = (fHisto->fTotal_flux * 1.0e9 / (G4double)fNevt) / fTestSphereSurfaceArea;
  G4cout << "Total absolute flux = " << fAbsolute_TotalFlux << G4endl;

}

void G4TARCHistoManager::EndOfRun() {
  G4cout << "G4TARCHistoManager ; End of run actions are started" << G4endl;
  G4cout << "fNevt = " << fNevt << G4endl;
  G4cout << "EndOfRun(), fEdepSum = " << fEdepSum << G4endl;
  G4cout << "======================================================================" << G4endl;

  FillRadialExperimentalData();

  fHisto->SetGeomParam(fShellNumber, fShellThickness, fInnerRadiusofShell, fOuterRadiusofShell);
  fHisto->SetExptDataParam(fMaxFluxData, fMaxFluenceData, fMaxTestFluxData, fExptEnergyBin,
                                fExptEnergyTables, fExptFluxTables, fExptFluxErrTables);
  //fHisto->StartProcessing(fAnalysisManager);

  //CreateNeutronFluxHisto();

  G4double x = ( G4double )fNevt;
  G4double perN = 1.0;
  if (fNevt > 0){
    x = perN = 1.0 / x;
  }
  TrackRun(x);  // track and get leaks
  NeutronRun(x); // neutron leak data and spectra

  GunParticleRun(x);  // beam distribution in target.

  //****************** This is for G4ParticleHPThermalScattering
  //fNHisto = 1;
/*
  for (G4int i = 0; i < fNHisto; i++)
    fHisto->ScaleH1(i, x);           // Normalize Histogram

  fHisto->Save();
  */
  fAnalysisManager->Write();
  fAnalysisManager->CloseFile();
}

void G4TARCHistoManager::BeginOfEvent() {
  fEdepEvt    = 0.0;
  fEdepEM     = 0.0;
  fEdepProton = 0.0;
  fEdepPI     = 0.0;
  fNsecondary = 0.0;
  fRangeSum   = G4ThreeVector(0.0, 0.0, 0.0);
  fStepSum    = 0.0;
  fDeltaSum   = 0.0;
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

  WriteEventRange( fRangeSum, fStepSum, fDeltaSum );
}


// Normalization for no hadron interaction in the target. EM and hadron elastic shuld be inactivated
void G4TARCHistoManager::AddNzero(const G4Track* myTrack, const G4Step* myStep ) {
  G4double cosTheta = myTrack->GetDynamicParticle()->GetMomentumDirection().x();
  // change x, y z
  if (cosTheta > 0.999) {
    fNzero++;
  }else {
    fNelastic++;
  }
}


void G4TARCHistoManager::GunParticleDistribution( const G4Track* myTrack, const G4Step* myStep) {
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

  G4double x = pos.x(); //- fAbsX0;
  G4double y = pos.y(); //- fAbsY0;
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
  const G4ParticleDefinition* pd = myTrack->GetDefinition();
  G4String name = pd->GetParticleName();
  G4double ke = myTrack->GetKineticEnergy();
  G4double TKin = myTrack->GetKineticEnergy();
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
      // fEventAction->AddNeutronStack();
      //fHisto->Fill(5, ke, 1.0);
      //G4cout << ke << G4endl;

      fAnalysisManager->FillNtupleDColumn(14, 0, log10(ke));
      fAnalysisManager->AddNtupleRow(14);
    } else if (pd == G4AntiProton::AntiProtonDefinition()){
      fNaproton++;
    } else if ( pd == G4PionPlus::PionPlusDefinition() ) {
      fNcpions++;
      // fHisto->Fill(6, ke, 1.0);
      //  fHisto->Fill(19, ke, 1.0);
    } else if ( pd == G4PionMinus::PionMinusDefinition()) {
      fNcpions++;
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
  const G4DynamicParticle*    dynSecondary = myTrack->GetDynamicParticle();
  G4double                    Tkin         = dynSecondary->GetKineticEnergy();

  if (myTrack->GetKineticEnergy()/MeV <= 0.0) {
    G4String pname  = pd->GetParticleName();
    G4double mass   = pd->GetPDGMass();
    return;
  }

  G4double en       = std::log10(myTrack->GetKineticEnergy()/eV);
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
         (zz > -fHLength && dir.z() > 0.0)
       )
  ) { // forward
      isLeaking = true;
      if (pd == fNeutron) { // && Tkin < fTcut) {
        ++fNneu_forw;
        //  fHisto->Fill(15, en,  1.0);
      }
  }
  if ( (
         // (xx < fAbsX0 && dir.x() < 0.0)
         // ||   (yy < fAbsY0 && dir.y() < 0.0)
         // ||
           (zz < fHLength && dir.z() < 0.0)
       )
  ) { // backward
    isLeaking = true;
    if (pd == fNeutron) { // && Tkin < fTcut)   {
       ++fNneu_back;
       //  fHisto->Fill(16, en,  1.0);
     }
  }
  if ( (
         // (std::abs(xx) <= fAbsX0 && (zz * dir.z()  + yy * dir.y()) > 0.0 )
         // ||   (std::abs(yy) <= fAbsY0 && (xx * dir.x()  + zz * dir.z()) > 0.0 )
         // ||
         //(std::abs(zz) <= -fHLength && (xx * dir.x()  + yy * dir.y()) > 0.0 )
         (std::abs(xx) > fHLength || std::abs(yy) > fHLength)
         && ((xx * dir.x()  + yy * dir.y()) > 0.0 )
       )
  ) { // side
    isLeaking = true;
    if (pd == fNeutron){ // && Tkin < fTcut)  {
      ++fNneu_leak;
      G4cout << "-------------------->   side " << fNneu_leak << G4endl;
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
void G4TARCHistoManager::NeutFinalState(const G4Track* myTrack, const G4Step* myStep) {
  fNsecondary++;
  G4double TKin = 0.0, Tspectra = 0.0, stepLen, time;
  G4int ix;
  size_t ii, jj;

  if (myTrack->GetDynamicParticle()->GetDefinition()->GetParticleName() == "neutron"){
    fNeutronSum+= 1.0;
    // future plan for neutron breed at breeder with X
    TKin = myTrack->GetDynamicParticle()->GetKineticEnergy();
    //stepLen = myStep->GetStepLength();
    time = myTrack->GetGlobalTime();
    for (ii = 0; ii < fNbin; ii++) {
      if (TKin <= fNEsecond.GetLowEdgeEnergy(ii)){
        jj = (ii == 0) ? ii : ii - 1;
        break;
      }
    }
    jj = (ii == fNbin) ? fNbin-1 : jj;
  }
}


void G4TARCHistoManager::TargetProfile(const G4Track* myTrack, const G4Step* myStep) {
  G4double TKinE = 0.0, TSpectr = 0.0, stepLen, time, xn;
  G4int ix = 0;

  if (myTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "neutron") {
    if (
           (myTrack->GetLogicalVolumeAtVertex()->GetName() == "blockB_log")
           && (  myTrack->GetVertexPosition().x() >= -35.0 && myTrack->GetVertexPosition().x() <= 35.0
           || myTrack->GetVertexPosition().y() >= -35.0 && myTrack->GetVertexPosition().y() <= 35.0
           || myTrack->GetVertexPosition().z() >= -1500.0 && myTrack->GetVertexPosition().z() <= -299.0
           )
           && (myTrack->GetGlobalTime()/nanosecond <= 10.0)
           && (myTrack->GetParentID() == 1)
     ) {
       fNeutronInit += 1.0;
       xn = myTrack->GetVertexPosition().x() + 0.5 * fRange;
       ix = G4int(xn / fLBin + 0.5);
       if (ix >= 0 && ix < fLMax) fGunParticleX[ix] += 1.0;

       xn = myTrack->GetVertexPosition().y() + 0.5 * fRange;
       ix = G4int(xn / fLBin + 0.5);
       if (ix >= 0 && ix < fLMax) fGunParticleY[ix] += 1.0;

       xn = myTrack->GetVertexPosition().z() + 0.5 * fRange;
       ix = G4int(xn / fLBin + 0.5);
       if (ix >= 0 && ix < fLMax) fGunParticleZ[ix] += 1.0;
     }
  }
}

void G4TARCHistoManager::AddEnergyTimeHole(const G4Track* myTrack, const G4Step* myStep) {
  G4double KE, myTime, myStepLength;
  size_t ii, jj, kk;

  std::vector<G4double> tmp;

  if (myTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "neutron"){
    KE = myTrack->GetDynamicParticle()->GetKineticEnergy();
    myStepLength = myStep->GetStepLength();
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

void G4TARCHistoManager::AddEnergyTime(const G4Track* myTrack, const G4Step* myStep) {
  G4double Tkin = 0.0, myTime, myGlobalTime, myStepLength;
  G4double Qenergy, Qtime, eVal;
  size_t ii, jj, eii, tjj;
  std::vector<G4double> tmp;

  if (myTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "neutron") {
    Tkin = myTrack->GetDynamicParticle()->GetKineticEnergy();
    eVal = Tkin / eV;
    myStepLength = myStep->GetStepLength();
    myTime = myTrack->GetProperTime();
    myGlobalTime = myTrack->GetGlobalTime();

    G4int fluxBin = eVal / fStepE;
    fluxBin = (fluxBin >= fMaxBin) ? fMaxBin - 1 : fluxBin;
    fNEfluxBin[fluxBin] += eVal;

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
    if (jj == fNbin) tjj = fNbin - 1;
    fET[tjj][eii] += 1.0;

    tmp.push_back(Qtime);
    tmp.push_back(Qenergy);
    fNSpectra.push_back(tmp);
    std::vector<G4double>().swap(tmp);
  }
}

void G4TARCHistoManager::WriteEventRange(G4ThreeVector rsum, G4double lsum, G4double dsum) {
  G4int ix;
  G4double x, y, z, rho, xbin, ybin, zbin, rbin, dbin, dEdx, length;

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

  div_t divresult = div((fEdepSum/MeV),1000);
  if (divresult.quot >= 1){
      trackout << " Total Eenergy Deposited = " << fEdepSum/GeV << " GeV"<< G4endl;
  } else {
      trackout << " Total Eenergy Deposited = " << fEdepSum/MeV << " MeV"<< G4endl;
  }

  fEdepSum  *= x;
  fEdepSum2 *= x;
  fEdepSum2 -= fEdepSum * fEdepSum;
  fEdepSum2  = (fEdepSum2 > 0.0) ? std::sqrt(fEdepSum2) : 0.0;

  trackout << "x = 1/fNevt = "               << x                                          << G4endl;

  divresult = div((fEdepSum/MeV),1000);
  if (divresult.quot >= 1){
      trackout << "x * Total Energy deposited = " << fEdepSum/GeV  << " GeV"         << G4endl;
  }else {
      trackout << "x * Total Energy deposited = " << fEdepSum/MeV  << " MeV"         << G4endl;
  }

  if (fPrimaryDef){
    trackout << "Beam Particle "                << fPrimaryDef->GetParticleName()    << " having ";
    divresult = div((fPrimaryKineticEnergy/MeV),1000);
    if (divresult.quot >= 1){
        trackout << "beam energy = " << fPrimaryKineticEnergy/GeV  << " GeV"         << G4endl;
    }else {
        trackout << "beam Energy = " << fPrimaryKineticEnergy/MeV  << " MeV"         << G4endl;
    }
  }
  trackout << "Number of Events = "             << fNevt                                    << G4endl;
  trackout << "                                                                        "  << G4endl;

  trackout << " Production in Target:"        << "                                     "  << G4endl;

  divresult = div((fEdepSum/MeV),1000);
  if (divresult.quot >= 1){
    trackout << std::setprecision(4) << "Mean Energy deposited = " << fEdepSum/GeV  << " GeV "
             << " RMS " << fEdepSum2/MeV << " MeV" "    "<< G4endl;
  }else {
    trackout << std::setprecision(4) << "Mean Energy deposited = " << fEdepSum/MeV  << " MeV "
             << " RMS = " << fEdepSum2/MeV << " MeV" "    "<< G4endl;
  }

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
  trackout << "                                                      " << "              "  << G4endl;

  trackout << " Leakage from the system: "                                                       << G4endl;
  trackout << std::setprecision(4) << "Average Number of forward Neutrons: "      << xneuF        << G4endl;
  trackout << std::setprecision(4) << "Average Number of reflected Neutrons: "    << xneuB        << G4endl;
  //trackout << std::setprecision(4) << "Average Number of other leaked Neutrons " << xNeutronLeak << G4endl;
  trackout << std::setprecision(4) << "Average Number of total leaked Neutrons: " <<  xNeutronLeak
                                                                                    + xneuF
                                                                                    + xneuB      << G4endl;
  trackout << std::setprecision(4) << "Average Number of leaked Protons: "        << xProtonLeak  << G4endl;
  trackout << std::setprecision(4) << "Average Number of leaked Pions: "          << xp0          << G4endl;
  trackout <<                                                                                     G4endl;



  G4double kEffective, rho, rat, react, perN=x;
  kEffective = (fNeutronInit!= 0.0) ? fNeutronSum / fNeutronInit : 0.0;
  rho        =  (kEffective != 0.0) ? (kEffective - 1.0) / kEffective : 0.0;  // reactivity :: deviation from criticality
  rat        = (kEffective != 0.0) ? std::log(kEffective) : 0.0;
  react      = rat / (1.0 + rat);

  trackout << " IMP Parameters : "    << G4endl;
  trackout << " Neutron_Init/p = "    << fNeutronInit* perN << ",  Neutron_Sum/p = " << fNeutronSum * perN << G4endl;
  trackout << " kEffective = "        << kEffective         << ", Rho = "            << rho                << G4endl;
  trackout << " Estimated reactivity = " << react                                       << G4endl             << G4endl;
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


  std::ofstream fFlu("FluenceData.dat", std::ios::out);
  std::cout << fFluence.size() << G4endl;
  for (std::size_t ii = 0; ii < fFluence.size(); ii++){
    for (std::size_t jj = 0; jj < fMaxEBin; jj++){
      fFlu << ii << "  " << jj << " " << fFluence[ii][jj] << G4endl;
    }
  }
  fFlu.close();

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
  for( size_t k = 0; k < fNbin; k++ ){
    for( size_t j = 0; j < fNbin; j++ ){
      // tenspectr<<perN*fNETsum[k]<<G4endl;
      if (perN * fET[j][k] != 0.0)
        tenspectr << (G4int)k << "    " << (G4int)j << "    " << perN * fET[j][k] << G4endl;
    }
  }
//----------------------------------------------------
  std::ofstream neutSpec("neutSpec.dat", std::ios::out);
  for (size_t ii = 0; ii < fNSpectra.size(); ii ++){
    neutSpec << fNSpectra[ii][0] << "   " << fNSpectra[ii][1] << G4endl;
  }
  neutSpec.close();
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

void G4TARCHistoManager::Fill(G4int id, G4double x, G4double w) {
  // fHisto->Fill(id, x, w);
}


  // // // // // // // // // // // // // // // // //
/*
  void G4TARCHistoManager::SetVerbose(G4int val)
  {
    fVerbose = val;
    fHisto->SetVerbose(val);
  }
*/
