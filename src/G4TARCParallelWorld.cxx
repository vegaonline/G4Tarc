/**********************************************************
 * @file      G4TARCParallelWorld.cxx
 * @author    Abhijit Bhattacharyya
 * @brief     Parallelization processes
 *********************************************************/
#include "G4TARCParallelWorld.hh"



G4TARCParallelWorld::G4TARCParallelWorld( G4String& tarcParallelWorld)
: G4VUserParallelWorld( tarcParallelWorld), fConstructed(false) {
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
  std::cout << " Define Blocks initialised." << std::endl;
}


void G4TARCParallelWorld::Construct() {
  if (fConstructed) return;
  DefineShellsBlocks();
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


  // Make the thinner shell bunches
  for (G4int i = 0; i < fShellNumber; i++) {
    G4Sphere* radShellSphere = new G4Sphere("shellSphere",
                            fInnerRadiusofShell[i], fOuterRadiusofShell[i],
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
  auto SDman = G4SDManager::GetSDMpointer();
  G4String TARCSDName, fltName, particleName, psName;
  G4SDParticleFilter* neutronFilter;
  G4SDParticleFilter* protonFilter;

  // Here first defining two detectors for neutron and proton study
  //----------------------------- SD 01 Neutron ----------------------------------------
  // Create a MultiFunction Detector for Neutron study in General
  TARCSDName = "TARCNeutronSD";
  fTARCNeutronDet = new G4MultiFunctionalDetector(TARCSDName);
  SDman->AddNewDetector(fTARCNeutronDet);
  // Then create filters to study
  neutronFilter = new G4SDParticleFilter(fltName="neutronFilter", particleName="neutron");
  fTARCNeutronDet->SetFilter(neutronFilter);

  //----------------------------- SD 02 Proton ----------------------------------------
  // Create a MultiFunction Detector for Proton study in General
  TARCSDName = "TARCProtonSD";
  fTARCProtonDet = new G4MultiFunctionalDetector(TARCSDName);
  SDman->AddNewDetector(fTARCProtonDet);
  // Then create filters to study
  protonFilter = new G4SDParticleFilter(fltName="protonFilter", particleName="proton");
  fTARCProtonDet->SetFilter(protonFilter);

  //------------- Next we need to add concerned LV to these detectors ----------------------
    for (std::vector<G4LogicalVolume*>::iterator it = fLVvector.begin(); it != fLVvector.end(); ++it) {
      if ( (*it)->GetName().find("Proton")!=std::string::npos){
        SetSensitiveDetector( (*it)->GetName(), fTARCProtonDet);
      } else{
        SetSensitiveDetector( (*it)->GetName(), fTARCNeutronDet);
      }
    }

    G4PSNofCollision*   scorer0 = new G4PSNofCollision(psName="Collisions");
    fTARCNeutronDet->RegisterPrimitive(scorer0);

    G4PSNofCollision*   scorer1 = new G4PSNofCollision(psName="CollWeight");
    scorer1->Weighted(true);
    fTARCNeutronDet->RegisterPrimitive(scorer1);


    G4PSPopulation*   scorer2N = new G4PSPopulation(psName="Population");
    fTARCNeutronDet->RegisterPrimitive(scorer2N);

    G4PSPopulation*   scorer2P = new G4PSPopulation(psName="Population");
    fTARCProtonDet->RegisterPrimitive(scorer2P);

    G4PSTrackCounter* scorer3 = new G4PSTrackCounter(psName="Track_Enter",fCurrent_In);
    fTARCNeutronDet->RegisterPrimitive(scorer3);
    G4PSTrackLength* scorer4 = new G4PSTrackLength(psName="Track_Length");
    fTARCNeutronDet->RegisterPrimitive(scorer4);

    G4PSTrackLength* scorer5 = new G4PSTrackLength(psName="Track_Length_Weighted");
    scorer5->Weighted(true);
    fTARCNeutronDet->RegisterPrimitive(scorer5);

    G4PSTrackLength* scorer6 = new G4PSTrackLength(psName="Track_Length_Weighted.KE");
    scorer6->Weighted(true);
    scorer6->MultiplyKineticEnergy(true);
    fTARCNeutronDet->RegisterPrimitive(scorer6);

    G4PSTrackLength* scorer7 = new G4PSTrackLength(psName="Track_Length_Weighted_By_Velocity");
    scorer7->Weighted(true);
    scorer7->DivideByVelocity(true);
    fTARCNeutronDet->RegisterPrimitive(scorer7);

    G4PSTrackLength* scorer8 = new G4PSTrackLength(psName="Track_Length_Weighted_KE_By_Velocity");
    scorer8->Weighted(true);
    scorer8->MultiplyKineticEnergy(true);
    scorer8->DivideByVelocity(true);
    fTARCNeutronDet->RegisterPrimitive(scorer8);
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
  istore->AddImportanceGeometryCell(1, gWorldVolumeCell);
  // set importance values and create scorers
  //  G4int number_shells = analysis->GetNumberShells();
  // G4int cell(fShellNumber);
  for (G4int cell = 0; cell < fShellNumber; cell++) {
    G4GeometryCell gCell = GetGeometryCell(cell);
    /*
    G4cout << " adding cell: " << cell << " replica: "
           << gCell.GetReplicaNumber() << " name: "
           << gCell.GetPhysicalVolume().GetName() << G4endl;
    */
    imp = 1;
    istore->AddImportanceGeometryCell(imp, gCell.GetPhysicalVolume(), cell);
  }
  return istore;
}
