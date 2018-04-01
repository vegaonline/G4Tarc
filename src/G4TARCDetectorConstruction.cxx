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
  return fWorldPhysVol;
}


// Sensitive Detectors
void G4TARCDetectorConstruction::ConstructSDandField() {
  auto sdManager = G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  sdManager->SetVerboseLevel(1);
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
