/*************************************************************************
 * @file      G4TARCDetectorConstruction.cxx
 * @author    Abhijit Bhattacharyya
 * @brief     This file is for Detector Geometry setup
 *************************************************************************/

#include "G4TARCDetectorConstruction.hh"

G4TARCDetectorConstruction::G4TARCDetectorConstruction()
: G4VUserDetectorConstruction(), fBeamBlock(0), fABlocks(0),
   fBBlocks(0), fCBlocks(0), fSampleTubes(0), fSampleSpheres(0),
  // fAllLead(0),
    fDetectorMessenger(0)   {
  fDetectorMessenger = new G4TARCDetectorMessenger(this);
}

G4TARCDetectorConstruction::G4TARCDetectorConstruction( G4String gdmlFileName ) //const G4GDMLParser& parser )
: G4VUserDetectorConstruction(), fDetectorMessenger(0), fLogiWorld(0),
  fBeamBlock(0), fABlocks(0), fBBlocks(0), fCBlocks(0), fSampleTubes(0),
  fSampleSpheres(0) //, fAllLead(0) {
{
    if (!fFileLoaded) SetReadFile(gdmlFileName);
    fDetectorMessenger = new G4TARCDetectorMessenger(this);
}

G4TARCDetectorConstruction::~G4TARCDetectorConstruction() {
  // if I uncomment then it crashes at the end after G4Kernel quitting.
  if (fDetectorMessenger) delete fDetectorMessenger;
}


void G4TARCDetectorConstruction::SetReadFile( const G4String& gdmlFile ) {
  fGdmlFileNAME = gdmlFile;
}


G4VPhysicalVolume* G4TARCDetectorConstruction::Construct() {
  if (!fFileLoaded) {
    fParser.Read(fGdmlFileNAME, true);
    fFileLoaded = true;
    fWorldPhysVol = fParser.GetWorldVolume();
  }
    //*************** testing logical volumes ****************
    fLVS = G4LogicalVolumeStore::GetInstance();
    int icount = 0;
    G4VSolid *thisBlock;
    G4ThreeVector pMin;
    G4ThreeVector pMax;

    for (fLVciter = fLVS->begin(); fLVciter !=fLVS->end(); fLVciter++) {
      G4cout << icount << "    " << (*fLVciter)->GetName() << "  " ;
      if ((*fLVciter)->GetName() == "blockB_log")
      {
        thisBlock = (*fLVciter)->GetSolid();
        thisBlock->BoundingLimits(pMin, pMax);
        setBoundsToProtonTarget(pMin, pMax);
        G4cout << "\n -------------> " << thisBlock->GetName() << " is bounded by "
               << "  pMin:  " << pbTargetMin << "  pMax: " << pbTargetMax
               << G4endl;
        G4cout << "\n  " << (*fLVciter)->GetName() << ":-->  "<< (*fLVciter)->GetNoDaughters()<< "   ";
        for (G4int iij = 0 ; iij < (*fLVciter)->GetNoDaughters(); iij++){
          G4cout << " Daughter volume = " << (*fLVciter)->GetDaughter(iij)->GetName()
                 << "  Logical of Daughter : "
                 <<  (*fLVciter)->GetDaughter(iij)->GetLogicalVolume()->GetName()
                 << G4endl;
        }

        G4GDMLAuxListType auxInfo = fParser.GetVolumeAuxiliaryInformation(*fLVciter);

        if (auxInfo.size()>0) {
          G4cout << "Auxiliary Information is found for Logical Volume :  "
			           << (*fLVciter)->GetName() << G4endl;
          print_aux(&auxInfo, "SS");
        } else {
          G4cout << "There is no auxilliary information" << G4endl;
        }
      }
      ++icount;
    }
    G4cout << G4endl;


  //G4cout << *(G4Material::GetMaterialTable() ) << G4endl;
  //G4cout << "Physical Volume Name: " << fWorldPhysVol->GetName()
  //       << "  Logical Volume under this Physical volume: " << fWorldPhysVol->GetLogicalVolume()->GetName()
  //       << G4endl;


  //set visualization attributes to world.
  auto sampleSphereVisAttRed1 = new G4VisAttributes(G4Colour(1.0, 0.2, 0.1));
  auto sampleTubeVisAttGreen1 = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
  auto blockAVisAttBlue1      = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  auto blockBVisAttBlue2      = new G4VisAttributes(G4Colour(0.1, 0.0, 0.8));
  auto blockCVisAttBlue3      = new G4VisAttributes(G4Colour(0.2, 0.0, 0.6));
  auto beamblockAVisAttRed    = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  auto leadVisAttRed1         = new G4VisAttributes(G4Colour(0.1, 0.2, 0.2));

  G4VPhysicalVolume* labPhy  = fWorldPhysVol->GetLogicalVolume()->GetDaughter(0);
  fLAB = labPhy->GetLogicalVolume();
  G4VPhysicalVolume* leadPhy = labPhy->GetLogicalVolume()->GetDaughter(0);
  //fAllLead = leadPhy->GetLogicalVolume();
  G4int noDaughters = leadPhy->GetLogicalVolume()->GetNoDaughters();
  for (int icount = 0; icount < noDaughters; icount++) {
    G4VPhysicalVolume* phys = leadPhy->GetLogicalVolume()->GetDaughter(icount);
    if (phys->GetLogicalVolume()->GetName() == "beamblockB_log" && !fLogiBeam){
      fBeamBlock = phys->GetLogicalVolume();
      fBeamBlock->SetVisAttributes(&beamblockAVisAttRed);
      fLogiBeam = true;
    }
    if (phys->GetLogicalVolume()->GetName() == "blockA_log" && !fLogiA){
      fABlocks = phys->GetLogicalVolume();
      fABlocks->SetVisAttributes(&blockAVisAttBlue1);
      fLogiA = true;
    }
    if (phys->GetLogicalVolume()->GetName() == "blockB_log" && !fLogiB){
      fBBlocks = phys->GetLogicalVolume();
      fBBlocks->SetVisAttributes(&blockBVisAttBlue2);
      fLogiB = true;
    }
    if (phys->GetLogicalVolume()->GetName() == "blockC_log" && !fLogiC){
      fCBlocks = phys->GetLogicalVolume();
      fCBlocks->SetVisAttributes(&blockCVisAttBlue3);
      fLogiC = true;
    }
    if (phys->GetLogicalVolume()->GetName() == "sampleTube_log" && !fLogiTube){
      fSampleTubes = phys->GetLogicalVolume();
      fSampleTubes->SetVisAttributes(&sampleTubeVisAttGreen1);
      fLogiTube = true;
    }
    if (phys->GetLogicalVolume()->GetName() == "sampleSphere_log" && !fLogiSphere){
      fSampleSpheres = phys->GetLogicalVolume();
      fSampleSpheres->SetVisAttributes(&sampleSphereVisAttRed1);
      fLogiSphere = true;
    }
  }
  fWorldPhysVol->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::GetInvisible());
  return fWorldPhysVol;
}


void G4TARCDetectorConstruction::CheckGDMLGeomTree(const G4VPhysicalVolume* masterVol){
  int k1 = masterVol->GetLogicalVolume()->GetNoDaughters();
  int icnt = 0;
  if (k1){
    for (int ii1 = 0; ii1 < k1; ii1++){
      G4VPhysicalVolume* phys0 = fWorldPhysVol->GetLogicalVolume()->GetDaughter(ii1);
      int k2 = phys0->GetLogicalVolume()->GetNoDaughters();
      G4cout << ii1 << "   " << phys0->GetLogicalVolume()->GetName() << " with " << k2 << " daughters." << G4endl;
      ++icnt;
      for (int ii2 = 0; ii2 < k2; ii2++){
        G4VPhysicalVolume* phys1 = phys0->GetLogicalVolume()->GetDaughter(ii2);
        int k3 = phys1->GetLogicalVolume()->GetNoDaughters();
        G4cout << ii2 << "   " << phys1->GetLogicalVolume()->GetName() << " with " << k3 << " daughters." << G4endl;
        ++icnt;
        for (int ii3 = 0; ii3 < k3; ii3++) {
          G4VPhysicalVolume* phys2 = phys1->GetLogicalVolume()->GetDaughter(ii3);
          int k4 = phys2->GetLogicalVolume()->GetNoDaughters();
          G4cout << ii3 << "   " << phys2->GetLogicalVolume()->GetName() << " with " << k4 << " daughters." << G4endl;
          ++icnt;
        }
      }
    }
  }
}

// Sensitive Detectors
void G4TARCDetectorConstruction::ConstructSDandField() {
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  static G4ThreadLocal G4bool initialized = false;
  if (!initialized){
    G4TARCCheckVolumeSD* fCheckVolumeSD = new G4TARCCheckVolumeSD("checkSD");
    (G4SDManager::GetSDMpointer())->AddNewDetector( fCheckVolumeSD);
    //fAllLead->SetSensitiveDetector(fCheckVolumeSD);
    // Let us consider that beamBlock is target initially
    fBeamBlock->SetSensitiveDetector(fCheckVolumeSD);
    fABlocks->SetSensitiveDetector(fCheckVolumeSD);
    //fBBlocks->SetSensitiveDetector(fCheckVolumeSD);
    fCBlocks->SetSensitiveDetector(fCheckVolumeSD);

    G4TARCTargetSD* fTargetSD = new G4TARCTargetSD("targetSD");
    (G4SDManager::GetSDMpointer())->AddNewDetector(fTargetSD);
    fBBlocks->SetSensitiveDetector(fTargetSD);
    //fBeamBlock->SetSensitiveDetector(fCheckVolumeSD);
    // We shall use sample tube later with rare earth for fission.
    //fSampleTubes->SetSensitiveDetector(fTargetSD);
    //fSampleSpheres->SetSensitiveDetector(fTargetSD);

    initialized = true;
  }
}

void G4TARCDetectorConstruction::print_aux(const G4GDMLAuxListType* auxInfoList, G4String prepend="|") {
	for(std::vector<G4GDMLAuxStructType>::const_iterator iaux = auxInfoList->begin(); iaux != auxInfoList->end(); iaux++ ) {
		G4String str=iaux->type;
		G4String val=iaux->value;
		G4String unit=iaux->unit;

		G4cout << prepend << str << " : " << val  << " " << unit << G4endl;

		if (iaux->auxList) print_aux(iaux->auxList, prepend + "|");
	}
	return;
}

void G4TARCDetectorConstruction::setBoundsToProtonTarget( const G4ThreeVector &pMin,  const G4ThreeVector &pMax){
  pbTargetMin = pMin;
  pbTargetMax = pMax;
}
