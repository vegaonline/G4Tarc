/**********************************************************
 * @file      G4TARCParallelWorld.cxx
 * @author    Abhijit Bhattacharyya
 * @brief     Parallelization processes
 *********************************************************/
#include "G4TARCParallelWorld.hh"

G4ThreadLocal G4bool G4TARCParallelWorld::fSDConstructed = false;


G4TARCParallelWorld::G4TARCParallelWorld( G4String& tarcParallelWorld)
  : G4VUserParallelWorld( tarcParallelWorld), fConstructed(false), fSerial(false) {
  DefineShellsBlocks();
}

G4TARCParallelWorld::~G4TARCParallelWorld() {
  std::vector<G4LogicalVolume*>().swap(fLVvector);
  std::vector<G4double>().swap(fInnerRadiusofShell);
  std::vector<G4double>().swap(fOuterRadiusofShell);
}

// Here defining concentric shells of finite thickness
void G4TARCParallelWorld::DefineShellsBlocks() {
  fHalfXBlockB           =     0.5 * 300 * mm;
  fHalfYBlockB           =     0.5 * 300 * mm;
  fHalfZBlockB           =     0.5 * 600 * mm;
  fHalfXVBox             =     0.5 * 150 * mm;
  fHalfYVBox             =     0.5 * 150 * mm;
  fHalfZVBox             =     0.5 * 300 * mm;
  fNewHalfZProt          =     0.5 * ((2.0 * fHalfZBlockB) / 3.0);
  fZposProt              = -fHalfZBlockB + fNewHalfZProt;
  fDiaMaxSphere          =  3300.0 * mm;
  fRadMaxSphere          =     0.5 * fDiaMaxSphere;
  fShellThickness        =     5.0 * cm;
  fRefShellThickness        =     2.0 * mm;
  fRefShellNumber = fRadiusReference.size();
  fRefShellOuterRad    = 457 * mm;
  fRefShellInnerRad   = fRefShellOuterRad - fRefShellThickness;
  fMinInnerRadiusofShell =    10.0 * mm;
  fMaxOuterRadiusofShell =  1500.0 * mm;
  fInnerRadProtonShell   =     0.0 * mm;   //
  fOuterRadProtonShell   =   300.0 * mm;   // These two were thought as a spherical 4Pi measurement for Proton
  // fShellNumber           = (G4int)((fMaxOuterRadiusofShell - fMinInnerRadiusofShell) / fShellThickness + 0.5);
  fRefShellNumber = fRadiusReference.size();

  for (G4int ii = 0; ii < fRefShellNumber; ii++) {
    G4double tmp = fRadiusReference[ii];
    fOuterRadiusofShell.push_back(tmp);
    tmp -=fRefShellThickness;
    fInnerRadiusofShell.push_back(tmp);
  }
  std::cout << " Define Blocks in ParallelWorld  initialised." << std::endl;
}


void G4TARCParallelWorld::Construct() {
  if (fConstructed) return;
  // A dummy material is used to fill the volulmes of the readout geometry.  Can we use Pb etc?
  G4Material* dummyMat = nullptr;
  G4String nameTmp;

  G4Colour col1 (1.0, 0.0, 0.0);  // red
  G4Colour col2 (1.0, 0.7, 0.0);
  G4Colour col3 (0.0, 1.0, 1.0);  //cyan
  G4VisAttributes* fAtt1 = new G4VisAttributes(col1);
  G4VisAttributes* fAtt2 = new G4VisAttributes(col2);
  G4VisAttributes* fAtt3 = new G4VisAttributes(col3);
  fAtt1->SetVisibility(true);
  fAtt2->SetVisibility(true);
  fAtt3->SetVisibility(true);
  fAtt3->SetForceSolid(true);

  // declaring the readout world
  // First the whole World
  ghostWorld = GetWorld();
  ghostWorldLog = ghostWorld->GetLogicalVolume();
  fLVvector.push_back(ghostWorldLog);
  fPVolumeStore.AddPVolume(G4GeometryCell(*ghostWorld, 0));

/*  This is being commented to run it smoothly without error and then it would be fixed.

  // Make a rectangular plate for proton hits and spallation point
  // blockB copy 50 equivalent : Z length is 300 mm for blockB. Here we consider 100 mm.
  G4Box* fTestBlockB50 = new G4Box("BoxPHit", fHalfXBlockB, fHalfYBlockB, fNewHalfZProt);
  fVBoxLogProton = new G4LogicalVolume(fTestBlockB50, dummyMat, "VB50Proton_log");
  fLVvector.push_back(fVBoxLogProton);
  nameTmp = GetCellName(0);
  // Here originally copy 50 was at 0,0,0 with zlength 300. Since we reduced to 100 and wish to
  // stick lowest Z bouondary to be at same position, origin of this block is shifted to -1000
  // so that
  fVBoxPVProton = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fZposProt), nameTmp, fVBoxLogProton, ghostWorld, false, 0);
  fVBoxLogProton->SetVisAttributes(fAtt1);
  G4GeometryCell fProtonCell(*fVBoxPVProton, 0);
  fPVolumeStore.AddPVolume(fProtonCell);

*/

  // Make the thinner shell bunches
  for (G4int i = 0; i < fRefShellNumber; i++) {
    G4Sphere* radShellSphere = new G4Sphere("shellSphere", fInnerRadiusofShell[i], fOuterRadiusofShell[i],
                            0.0*deg, 360.0*deg, 0.0*deg, 180.0*deg);
    nameTmp = "radial_shell_log_" + std::to_string(i);
    fShellLog = new G4LogicalVolume(radShellSphere, dummyMat, nameTmp);
    fLVvector.push_back(fShellLog);
    nameTmp = GetCellName(i);
    fShellPhys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), nameTmp, fShellLog, ghostWorld, false, i);
    fShellLog->SetVisAttributes(fAtt3);
    G4GeometryCell radCell(*fShellPhys, i);
    fPVolumeStore.AddPVolume(radCell);
  }
  fConstructed = true;
}


G4String G4TARCParallelWorld::GetCellName(G4int i) {
  G4String tmp = "";
  tmp = (i < 10) ? "0" + std::to_string(i) : std::to_string(i);
  return "cell_" + tmp;
}

G4GeometryCell G4TARCParallelWorld::GetGeometryCell(G4int i) {
  G4String name(GetCellName(i));
  const G4VPhysicalVolume* tmpPV =fPVolumeStore.GetPVolume(name);
  if (tmpPV) {
    return G4GeometryCell(*tmpPV, 0);
  } else {
    return G4GeometryCell(*ghostWorld, -2);
  }
}


void G4TARCParallelWorld::ConstructSD() {
  if (!fSDConstructed) {
    auto SDman = G4SDManager::GetSDMpointer();
    G4String fltName, fParticleName, psName, filterName;

    fNeutronFilter = new G4SDParticleFilter(filterName="neutronFilter", fParticleName = "neutron");
    fGammaFilter = new G4SDParticleFilter(filterName="gammaFilter", fParticleName = "gamma");
    fNeutronFilter = new G4SDParticleFilter(filterName="neutronFilter", fParticleName = "neutron");
    fElectronFilter = new G4SDParticleFilter(filterName="ElectronFilter", fParticleName = "e-");
    fPositronFilter = new G4SDParticleFilter(filterName="PositronFilter", fParticleName = "e+");
    /*
    fEPFilter = new G4SDParticleFilter(filterName="EPFilter");
    fEPFilter->add(fParticleName="e-");
    fEPFilter->add(fParticleName="e+");
    */
    fProtonFilter = new G4SDParticleFilter(filterName="ProtonFilter", fParticleName = "proton");
    fAntiProtonFilter = new G4SDParticleFilter(filterName="AntiProtonFilter", fParticleName = "anti_proton");
    fMuPlusFilter = new G4SDParticleFilter(filterName="MuPFilter", fParticleName = "mu+");
    fMuMinusFilter = new G4SDParticleFilter(filterName="MuMFilter", fParticleName = "mu-");
    fMuonFilter = new G4SDParticleFilter(filterName="MuonFilter");
    fMuonFilter->add(fParticleName="mu+");
    fMuonFilter->add(fParticleName="mu-");
    fPiPlusFilter = new G4SDParticleFilter(filterName="PiPFilter", fParticleName = "pi+");
    fPiMinusFilter = new G4SDParticleFilter(filterName="PiMFilter", fParticleName = "pi-");
    fPiZeroFilter = new G4SDParticleFilter(filterName="Pi0Filter", fParticleName = "pi0");
    fPionFilter = new G4SDParticleFilter(filterName="PionFilter");
    fPionFilter->add(fParticleName="pi+");
    fPionFilter->add(fParticleName="pi-");
    fPionFilter->add(fParticleName="pi0");
    fDeuteronFilter = new G4SDParticleFilter(filterName="DeuteronFilter", fParticleName = "deuteron");
    fTritonFilter = new G4SDParticleFilter(filterName="TritonFilter", fParticleName = "triton");
    fAlphaFilter = new G4SDParticleFilter(filterName="AlphaFilter", fParticleName = "alpha");
    fHe3Filter = new G4SDParticleFilter(filterName="He3Filter", fParticleName = "He3");

    G4MultiFunctionalDetector* fTARCNeutSD = new G4MultiFunctionalDetector("TARCNeut");
    SDman->AddNewDetector(fTARCNeutSD);
    G4SDParticleFilter* fNeutronFilter = new G4SDParticleFilter("NeutronFilter", fParticleName = "neutron");
    fTARCNeutSD->SetFilter(fNeutronFilter);

    G4MultiFunctionalDetector* fTARCNeutSRCSD = new G4MultiFunctionalDetector("TARCNeutSRC");
    SDman->AddNewDetector(fTARCNeutSRCSD);
    fTARCNeutSRCSD->SetFilter(fNeutronFilter);



    /*
    // Here first defining two detectors for neutron and proton study
    //----------------------------- SD 01 Neutron ----------------------------------------
    // Create a MultiFunction Detector for Neutron study in General
    TARCSDName = "TARCNeutronSD";
    fTARCNeutronDet = new G4MultiFunctionalDetector(TARCSDName);
    SDman->AddNewDetector(fTARCNeutronDet);
    // Then create filters to study
    neutronFilter = new G4SDParticleFilter(fltName="neutronFilter", fParticleName="neutron");
    fTARCNeutronDet->SetFilter(neutronFilter);

    //----------------------------- SD 02 Proton ----------------------------------------
    // Create a MultiFunction Detector for Proton study in General
    TARCSDName = "TARCProtonSD";
    fTARCProtonDet = new G4MultiFunctionalDetector(TARCSDName);
    SDman->AddNewDetector(fTARCProtonDet);
    // Then create filters to study
    protonFilter = new G4SDParticleFilter(fltName="protonFilter", fParticleName="proton");
    fTARCProtonDet->SetFilter(protonFilter);
    */

    G4String sslog = "sampleSphere_log";
    G4String sslog2 = "sampleSphere_log2";

    //------------- Next we need to add concerned LV to these detectors ----------------------
    for (std::vector<G4LogicalVolume*>::iterator it = fLVvector.begin(); it != fLVvector.end(); ++it) {
      if (((*it)->GetName().find(sslog) !=std::string::npos) ||  ((*it)->GetName().find(sslog2)!=std::string::npos)) {
        SetSensitiveDetector( (*it)->GetName(), fTARCNeutSRCSD);
      } else {
        SetSensitiveDetector( (*it)->GetName(), fTARCNeutSD);
      }
    }

    G4PSNofCollision*   scorer0 = new G4PSNofCollision(psName="Collisions");
    fTARCNeutSD->RegisterPrimitive(scorer0);
    fTARCNeutSRCSD->RegisterPrimitive(scorer0);

    G4PSNofCollision*   scorer1 = new G4PSNofCollision(psName="CollWeight");
    scorer1->Weighted(true);
    fTARCNeutSD->RegisterPrimitive(scorer1);
    fTARCNeutSRCSD->RegisterPrimitive(scorer1);

    G4PSPopulation*   scorer2N = new G4PSPopulation(psName="Population");
    fTARCNeutSD->RegisterPrimitive(scorer2N);
    fTARCNeutSRCSD->RegisterPrimitive(scorer2N);

    G4PSTrackCounter* scorer3 = new G4PSTrackCounter(psName="Track_Enter",fCurrent_In);
    fTARCNeutSD->RegisterPrimitive(scorer3);
    fTARCNeutSRCSD->RegisterPrimitive(scorer3);

    G4PSTrackLength* scorer4 = new G4PSTrackLength(psName="Track_Length");
    fTARCNeutSD->RegisterPrimitive(scorer4);
    fTARCNeutSRCSD->RegisterPrimitive(scorer4);

    G4PSTrackLength* scorer5 = new G4PSTrackLength(psName="Track_Length_Weighted");
    scorer5->Weighted(true);
    fTARCNeutSD->RegisterPrimitive(scorer5);
    fTARCNeutSRCSD->RegisterPrimitive(scorer5);

    G4PSTrackLength* scorer6 = new G4PSTrackLength(psName="Track_Length_Weighted.KE");
    scorer6->Weighted(true);
    scorer6->MultiplyKineticEnergy(true);
    fTARCNeutSD->RegisterPrimitive(scorer6);
    fTARCNeutSRCSD->RegisterPrimitive(scorer6);

    G4PSTrackLength* scorer7 = new G4PSTrackLength(psName="Track_Length_Weighted_By_Velocity");
    scorer7->Weighted(true);
    scorer7->DivideByVelocity(true);
    fTARCNeutSD->RegisterPrimitive(scorer7);
    fTARCNeutSRCSD->RegisterPrimitive(scorer7);

    G4PSTrackLength* scorer8 = new G4PSTrackLength(psName="Track_Length_Weighted_KE_By_Velocity");
    scorer8->Weighted(true);
    scorer8->MultiplyKineticEnergy(true);
    scorer8->DivideByVelocity(true);
    fTARCNeutSD->RegisterPrimitive(scorer8);
    fTARCNeutSRCSD->RegisterPrimitive(scorer8);
    
    fSDConstructed = true;
  }
}


void G4TARCParallelWorld::SetSerialGeometry(G4bool ser) {
  if (fSerial = ser) return;
  fSerial = ser;
  if (!fConstructed) return;
}

G4VIStore* G4TARCParallelWorld::CreateImportanceStore(){
  G4cout << " G4TARCParallelWorld:: Creating Importance Store " << G4endl;
  if (!fPVolumeStore.Size())  {
    G4Exception("G4TARCParallelWorld::CreateImportanceStore"
               ,"Testing...",RunMustBeAborted
               ,"no physical volumes created yet!");
  }
  // creating and filling the importance store
  //  G4IStore *istore = new G4IStore(*fWorldVolume);
  G4IStore *istore = G4IStore::GetInstance(GetName());
  G4GeometryCell gWorldVolumeCell(GetWorldVolumeAddress(), 0);
  G4double imp = 1;
  istore->AddImportanceGeometryCell(imp, gWorldVolumeCell);
  // set importance values and create scorers
  //  G4int number_shells = analysis->GetNumberShells();
  // G4int cell(fShellNumber);
  for (G4int cell = 0; cell < fRefShellNumber; cell++) {
    G4GeometryCell gCell = GetGeometryCell(cell);
    imp = 1;
    istore->AddImportanceGeometryCell(imp, gCell.GetPhysicalVolume(), cell);
  }
  return istore;
}
