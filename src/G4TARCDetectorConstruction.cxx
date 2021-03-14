/*************************************************************************
 * @file      G4TARCDetectorConstruction.cxx
 * @author    Abhijit Bhattacharyya
 * @brief     This file is for Detector Geometry setup
 *************************************************************************/

#include "G4TARCDetectorConstruction.hh"

G4ThreadLocal G4bool G4TARCDetectorConstruction::fConstructedSDandField = false;

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
  if (fDetectorMessenger) delete fDetectorMessenger;
}


void G4TARCDetectorConstruction::SetReadFile( const G4String& gdmlFile ) {
  fGdmlFileNAME = gdmlFile;
}


G4VPhysicalVolume* G4TARCDetectorConstruction::Construct() {
  if (!fConstructed) {
    if (!fFileLoaded) {
      fParser.Read(fGdmlFileNAME, true);
      fFileLoaded = true;
      fWorldPhysVol = fParser.GetWorldVolume();
    }
	G4cout << "GDML loaded" << G4endl;
    fConstructed = true;

    addTransU();   // define extra material(s) for testSphere LV

    // Declaring sampleSphere volume with Tc99
    fLVS = G4LogicalVolumeStore::GetInstance();
    G4String sslog  = "sampleSphere_log";
    G4String sslog2 = "sampleSphere_log2";
    for(fLVciter = fLVS->begin(); fLVciter != fLVS->end(); fLVciter++) {
      G4String LVName = (*fLVciter)->GetName();
      G4cout << "LVName: " << LVName << G4endl;
      if ((LVName.find(sslog) != std::string::npos) || (LVName.find(sslog2) != std::string::npos)){
	      (*fLVciter)->SetMaterial(fTc99);
	      (*fLVciter)->SetVisAttributes(new G4Colour(1.0, 0.5, 0.5));
      }
      fLVvectorMini.push_back(*fLVciter);
    }
  }

  // G4cout <<*(G4Material::GetMaterialTable()) << G4endl;
  return fWorldPhysVol;
}


void G4TARCDetectorConstruction::addTransU(){
  G4int ncomponents;
  G4double atMass;
  G4int atNum, nNucleons;
  G4double density, abundance, massfraction;
  G4String name, symbol;

  //  99Tc
  atMass = 98.906254 * g / mole;
  atNum = 43;
  nNucleons = 99;
  density = 11.0 * g / cm3;
  G4Isotope* Tc99Iso = new G4Isotope("Tc", atNum, nNucleons, atMass);
  G4Element* Tc99Ele = new G4Element("Tecnicium-99", "99Tc", ncomponents=1);
  Tc99Ele->AddIsotope(Tc99Iso, abundance=100.0*perCent);
  fTc99 = new G4Material("Tc", density, ncomponents=1);
  fTc99->AddElement(Tc99Ele, massfraction=100.0*perCent);
}
