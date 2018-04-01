#include "G4TARCSteppingAction.hh"

G4TARCSteppingAction::G4TARCSteppingAction(G4TARCEventAction* anEvent)
: G4UserSteppingAction(), fEventAction(anEvent){
  DefineShellsBlocks();
  startEnergy = 0.0;
  flag = false;
  number_generations = 0;

  radii = G4DataVector(fShellNumber, 0.0);
  for (G4int irad = 0; irad < fShellNumber; irad++){
    radii[irad] = fInnershellRadius + irad * fShellThickness;
  }
}


// Here defining concentric shells of finite thickness
void G4TARCSteppingAction::DefineShellsBlocks() {
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
}
