/*************************************************************************
 * @file      G4TARCDetectorConstruction.cxx
 * @author    Abhijit Bhattacharyya
 * @brief     This file is for Detector Geometry setup
 *************************************************************************/

#include "G4TARCDetectorConstruction.hh"

G4TARCDetectorConstruction::G4TARCDetectorConstruction()
: G4VUserDetectorConstruction(),
  //fBeamBlock(0), fABlocks(0),   fBBlocks(0), fCBlocks(0), fSampleTubes(0), fSampleSpheres(0), fAllLead(0),
    fDetectorMessenger(0)   {
  fDetectorMessenger = new G4TARCDetectorMessenger(this);
}

G4TARCDetectorConstruction::G4TARCDetectorConstruction( G4String gdmlFileName ) //const G4GDMLParser& parser )
: G4VUserDetectorConstruction(), fDetectorMessenger(0)
  //, fLogiWorld(0),  fBeamBlock(0), fABlocks(0), fBBlocks(0), fCBlocks(0), fSampleTubes(0),  fSampleSpheres(0) //, fAllLead(0) {
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

  // Declaring SD for whole volume
  fLVS = G4LogicalVolumeStore::GetInstance();
  for(fLVciter = fLVS->begin(); fLVciter != fLVS->end(); fLVciter++) {
    G4String LVName = (*fLVciter)->GetName();
    std::size_t found1 = LVName.find("blockA_log");
    std::size_t found2 = LVName.find("blockB_log");
    std::size_t found3 = LVName.find("blockC_log");
    std::size_t found4 = LVName.find("Sphere_log");
    std::size_t found5 = LVName.find("Tube_log");
    if ((found1 != std::string::npos) || (found2 != std::string::npos)
     || (found3 != std::string::npos) || (found4 != std::string::npos)
     || (found5 != std::string::npos)
    ){
      fLVvectorMini.push_back(*fLVciter);
      //G4cout << (*fLVciter)->GetName() << G4endl;
    }
  }
  return fWorldPhysVol;
}


// Sensitive Detectors
void G4TARCDetectorConstruction::ConstructSDandField() {
  //auto sdManager = G4SDManager::GetSDMpointer();
  //sdManager->SetVerboseLevel(1);
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  static G4ThreadLocal G4bool initialized = false;

  if (!initialized){
    G4TARCVolumeSD* fAllBlocksSD = new G4TARCVolumeSD("AllVolSD");
    (G4SDManager::GetSDMpointer())->AddNewDetector(fAllBlocksSD);

    G4TARCTargetSD* fBlockBSD = new G4TARCTargetSD("BVolSD");
    (G4SDManager::GetSDMpointer())->AddNewDetector(fBlockBSD);

    for (std::vector<G4LogicalVolume*>::iterator it = fLVvectorMini.begin(); it != fLVvectorMini.end(); ++it) {
      if ((*it)->GetName().find("blockB_log")!=std::string::npos){
        SetSensitiveDetector( (*it)->GetName(), fBlockBSD);
      } else {
        SetSensitiveDetector( (*it)->GetName(), fAllBlocksSD);
      }
    }
    initialized = true;
  }
}
