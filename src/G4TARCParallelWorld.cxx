/**********************************************************
 * @file      G4TARCParallelWorld.cxx
 * @author    Abhijit Bhattacharyya
 * @brief     Parallelization processes
 *********************************************************/
#include "G4TARCParallelWorld.hh"



G4TARCParallelWorld::G4TARCParallelWorld( G4String& tarcParallelWorld)
: G4VUserParallelWorld( tarcParallelWorld), fConstructed(false) {

}

G4TARCParallelWorld::~G4TARCParallelWorld() {}

void G4TARCParallelWorld::Construct() {
  if (fConstructed) return;
  // A dummy material is used to fill the volulmes of the readout geometry.
  G4Material* dummyMat = nullptr;

  // declaring the readout world
  ghostWorld = GetWorld();
  ghostWorldLog = ghostWorld->GetLogicalVolume();
  fConstructed = true;

  // target block : blockB copy 50 at (0, 0, 0) 30 cms X 30 cms X 60 cms
  LeadTargetS = new G4Box("LeadTarget", fHalfXblockB, fHalfYblockB, fHalfZblockB);
  LeadTargetLV = new G4LogicalVolume (LeadTargetS, dummyMat, "target_log", 0, 0, 0);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), LeadTargetLV, "target_phys", ghostWorldLog, false, 0);

  // Vritual box around holes. Let is place around hole 4 (say)
  VirtualBoxS =  new G4Box("VBox4holes", fHalfXVBox, fHalfYVBox, fHalfZVBox);
  VirtualBox1LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox1_log", 0, 0, 0);
  VB1PV = new G4PVPlacement(0, G4ThreeVector(  1050.0,     0.0, fZposVBox), VirtualBox1LV,  "vbox1_phys", ghostWorldLog, false,  0);
  VirtualBox2LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox2_log", 0, 0, 0);
  VB2PV = new G4PVPlacement(0, G4ThreeVector(   600.0,   300.0, fZposVBox), VirtualBox2LV,  "vbox2_phys", ghostWorldLog, false,  0);
  VirtualBox3LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox3_log", 0, 0, 0);
  VB3PV = new G4PVPlacement(0, G4ThreeVector(   150.0,     0.0, fZposVBox), VirtualBox3LV,  "vbox3_phys", ghostWorldLog, false,  0);
  VirtualBox4LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox4_log", 0, 0, 0);
  VB4PV = new G4PVPlacement(0, G4ThreeVector(     0.0, -1500.0, fZposVBox), VirtualBox4LV,  "vbox4_phys", ghostWorldLog, false,  0);
  VirtualBox5LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox5_log", 0, 0, 0);
  VB5PV = new G4PVPlacement(0, G4ThreeVector(     0.0,  -600.0, fZposVBox), VirtualBox5LV,  "vbox5_phys", ghostWorldLog, false,  0);
  VirtualBox6LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox6_log", 0, 0, 0);
  VB6PV = new G4PVPlacement(0, G4ThreeVector(     0.0,   600.0, fZposVBox), VirtualBox6LV,  "vbox6_phys", ghostWorldLog, false,  0);
  VirtualBox7LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox7_log", 0, 0, 0);
  VB7PV = new G4PVPlacement(0, G4ThreeVector(     0.0,   900.0, fZposVBox), VirtualBox7LV,  "vbox7_phys", ghostWorldLog, false,  0);
  VirtualBox8LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox8_log", 0, 0, 0);
  VB8PV = new G4PVPlacement(0, G4ThreeVector(     0.0,  1200.0, fZposVBox), VirtualBox8LV,  "vbox8_phys", ghostWorldLog, false,  0);
  VirtualBox9LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox9_log", 0, 0, 0);
  VB9PV = new G4PVPlacement(0, G4ThreeVector(     0.0,  1500.0, fZposVBox), VirtualBox9LV,  "vbox9_phys", ghostWorldLog, false,  0);
  VirtualBox10LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox10_log", 0, 0, 0);
  VB10PV = new G4PVPlacement(0, G4ThreeVector( -450.0,     0.0, fZposVBox), VirtualBox10LV,  "vbox10_phys", ghostWorldLog, false, 0);
  VirtualBox11LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox11_log", 0, 0, 0);
  VB11PV = new G4PVPlacement(0, G4ThreeVector( -600.0,   300.0, fZposVBox), VirtualBox11LV,  "vbox11_phys", ghostWorldLog, false, 0);
  VirtualBox12LV = new G4LogicalVolume(VirtualBoxS, dummyMat, "vBox12_log", 0, 0, 0);
  VB12PV = new G4PVPlacement(0, G4ThreeVector(-1050.0,     0.0, fZposVBox), VirtualBox12LV,  "vbox12_phys", ghostWorldLog, false,  0);

  // Create virtual sphere with max dia for total neutron flux within it
  // VirtualSphereS = new G4Orb("VirtualSphereS", fRadMaxSphere );
  // VirtualSphereLV = new G4LogicalVolume(VirtualSphereS, dummyMat, "vsph_log", 0, 0);
  // new G4PVPlacement(0, G4ThreeVector(0, 0, 0), VirtualSphereLV, "vsph_phys", ghostWorldLog, false, 0);

}

void G4TARCParallelWorld::ConstructSD() {
  G4String leadTargetName = "/tarc/leadTargetSD";
  G4TARCleadTargSD* leadTargSD = new G4TARCleadTargSD(leadTargetName);
  G4SDManager::GetSDMpointer()->AddNewDetector(leadTargSD);
  LeadTargetLV->SetSensitiveDetector(leadTargSD);


  G4String vBoxwithHoleName = "VirtualBox";
  G4TARCVirtualSD* virtSD = new G4TARCVirtualSD(vBoxwithHoleName);
  G4SDManager::GetSDMpointer()->AddNewDetector(virtSD);
  VirtualBox1LV->SetSensitiveDetector(virtSD);
  VirtualBox2LV->SetSensitiveDetector(virtSD);
  VirtualBox3LV->SetSensitiveDetector(virtSD);
  VirtualBox4LV->SetSensitiveDetector(virtSD);
  VirtualBox5LV->SetSensitiveDetector(virtSD);
  VirtualBox6LV->SetSensitiveDetector(virtSD);
  VirtualBox7LV->SetSensitiveDetector(virtSD);
  VirtualBox8LV->SetSensitiveDetector(virtSD);
  VirtualBox9LV->SetSensitiveDetector(virtSD);
  VirtualBox10LV->SetSensitiveDetector(virtSD);
  VirtualBox11LV->SetSensitiveDetector(virtSD);
  VirtualBox12LV->SetSensitiveDetector(virtSD);

/*
  G4String virtualdetName = "VirtualSphere";
  G4TARCvirtualSD* virtSD = new G4TARCvirtualSD(virtualdetName);
  G4SDManager::GetSDMpointer()->AddNewDetector(virtSD);
  VirtualSphereLV->SetSensitiveDetector(virtSD);
*/
}
