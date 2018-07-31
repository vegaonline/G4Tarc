#include "G4TARCSteppingAction.hh"

G4TARCSteppingAction::G4TARCSteppingAction(G4TARCEventAction* anEvent)
: G4UserSteppingAction()//, fEventAction(anEvent), fHisto(0)
{
  fEventAction = anEvent;
  //fHisto = G4TARCHistoManager::GetPointer();
  fRefShellThickness = 2.0 * mm;
  fRefShellOuterRad = 456.0 * mm;
  fRefShellInnerRad = fRefShellOuterRad - fRefShellThickness;
  fMaxOuterRadiusofShell = 1500.0 * mm;
  fMinInnerRadiusofShell = 10.0 * mm;
  fEnergy0 = 0.0;
  fNumberGenerations = 0;
  fShellNumber = (G4int)((fMaxOuterRadiusofShell - fMinInnerRadiusofShell) / fShellThickness + 0.5);

}

void G4TARCSteppingAction::UserSteppingAction(const G4Step* myStep){
  //fHisto->ProcessStepping(myStep);
  ProcessStepping(myStep);
}


// User Stepping Action
void G4TARCSteppingAction::ProcessStepping(const G4Step* myStep){

  G4Track* myTrack = myStep->GetTrack();
  G4int StepNo = myTrack->GetCurrentStepNumber();
  G4double fParticleEnergy = myStep->GetPreStepPoint()->GetKineticEnergy();
  G4double fParticleTime   = myStep->GetPreStepPoint()->GetGlobalTime();
  G4double fParticleMomentum = myStep->GetPreStepPoint()->GetMomentum().mag();
  //G4double fZMomentum = myStep->GetPreStepPoint()->GetMomentum().z();
  G4double angle = myStep->GetPreStepPoint()->GetMomentum().angle(myStep->GetPreStepPoint()->GetPosition());
  G4double cosAngle = std::abs(cos(angle));
  //G4int thisTrackID = myStep->GetTrack()->GetTrackID();
  //G4int parentTrackID = myStep->GetTrack()->GetParentID();
  G4int thisTrackID = myTrack->GetTrackID();
  G4int parentTrackID = myTrack->GetParentID();

  G4ParticleDefinition* fParticleType = myTrack->GetDefinition();
  G4String fParticleName = fParticleType->GetParticleName();
  G4double primEnergy = 0.0;

  fEventAction->analysePS(fParticleEnergy, fParticleName, fParticleMomentum);      //  , fParticleTime, fParticleMomentum, zMomentum);


  if (StepNo == 1){
    if (thisTrackID == 1){
      fParentEnergy.clear();
      fParentParticle.clear();
      fParentParticleID.clear();
      fNumberGenerations = 0;
    }
    if ( fParticleName == "neutron"){
      fEnergy0 = fParticleEnergy;
      fTime0 = fParticleTime;
      flag = true;
    }
  }

  fParentEnergy[thisTrackID] = fParticleEnergy;
  fParentParticle[thisTrackID] = fParticleName;
  fParentParticleID[thisTrackID] = parentTrackID;

  G4bool reduced_tally = false;

  if (thisTrackID == 1 && StepNo == 1) primEnergy = fParticleEnergy;
  if (thisTrackID != 1 && StepNo == 1) {
    G4int tempID = thisTrackID;
    fNumberGenerations = 1;
    while(fParentParticleID[tempID] != 1){
      tempID = fParentParticleID[tempID];
      ++fNumberGenerations;
    }
    if (fParentParticle[parentTrackID] == "neutron" && fParticleName == "neutron"){
      reduced_tally = true;
      fParentParticle.erase(parentTrackID);
    }
    fEventAction->analyseSecondaries (fParticleEnergy, fParticleName, fParticleTime, fParticleMomentum, parentTrackID, primEnergy,
      fParentEnergy[parentTrackID], fParentParticle[parentTrackID], reduced_tally, fNumberGenerations);
  }

  if (fParticleName == "neutron"){
    fEventAction->NeutronEnergyTime(fParticleEnergy, fParticleTime, fEnergy0);
  } else {
    if (fParticleName == "Pb207" || fParticleName == "Pb208")  fEventAction->otherEnergyTime(fParticleEnergy, fParticleTime, fEnergy0);
  }
  G4double radiusPre = myStep->GetPreStepPoint()->GetPosition().mag();
  G4double radiusPost = myStep->GetPostStepPoint()->GetPosition().mag();
  //G4double zPos = myStep->GetPreStepPoint()->GetPosition().z();
  G4double StepLength = myStep->GetStepLength();
  G4String vol = myStep->GetTrack()->GetVolume()->GetName();

  G4TouchableHistory* thePreTouchable  = (G4TouchableHistory*) (myStep->GetPreStepPoint()->GetTouchable());
  G4TouchableHistory* thePostTouchable = (G4TouchableHistory*) (myStep->GetPostStepPoint()->GetTouchable());

G4cout << "----> HERE" << G4endl;

  if (fParticleName == "neutron"){
    if (myStep->GetTrack()->GetNextVolume()){
      G4String PreVol = thePreTouchable->GetVolume()->GetName();
      G4String PostVol = thePostTouchable->GetVolume()->GetName();

      if (PreVol == "lab_phys" && PostVol == "world_log_PV"){
        fEventAction->exitingTally(true, fParticleEnergy);
      }
      if (PostVol == "world_log_PV"){
        fEventAction->exitingTallyCheck(true);
      }
    }



    G4bool pre_inside = false;
    G4bool post_inside = false;

    if  ((radiusPre <= (fRefShellOuterRad + fMyTol)) && (radiusPre >= (fRefShellInnerRad - fMyTol)) ) pre_inside = true;
    if  ((radiusPost <= (fRefShellOuterRad + fMyTol)) && (radiusPost >= (fRefShellInnerRad - fMyTol))) post_inside = true;
    if (pre_inside && post_inside) fEventAction->analyseNeutronShellFluence(fParticleEnergy, StepLength);

    for (G4int ishell = 0; ishell < fRefShellNumber; ishell++){
      G4bool pre_inside_radial = false;
      G4bool post_inside_radial = false;
      G4double radOut = fOuterRadiusofShell[ishell];
      G4double radIn  = fInnerRadiusofShell[ishell];
      if ((radiusPre <= (radOut + fMyRadTol)) && (radiusPre >= (radIn - fMyRadTol))) pre_inside_radial = true;
      if ((radiusPost <= (radOut + fMyRadTol)) && (radiusPost >= (radIn - fMyRadTol)) ) post_inside_radial = true;
      if (pre_inside_radial && post_inside_radial) fEventAction->analyseNeutronRadialFluence(fParticleEnergy, StepLength, ishell);  // fParticleTime, StepLength, ishell);

    }

    if (vol == "sample_phys" || vol == "sampleTube_phys" || vol == "sample_phys2"){
      G4double radValue = fRefShellOuterRad;
      fEventAction->analyseNeutronFluence(fParticleEnergy, fParticleTime,  thisTrackID, radValue, StepLength,  parentTrackID, primEnergy,  fParticleName);
    }

    //for (std::size_t ii = 0; ii < fFluxRadTables.size(); ++ii){
    for (G4int ii = 0; ii < fMaxRadCount; ++ii){
      G4double radValue = fExptRadiiTables[8][ii];   // fFluxRadTables[ii] / 10.0;
      if ( (radiusPre < radValue && radiusPost > radValue) ||(radiusPre > radValue && radiusPost < radValue)){
          //G4cout << ii << " FluxRad  " << radValue << "  RadPre " << radiusPre <<  " radPost " << radiusPost << G4endl;
        G4double radiusL = radValue;
        fEventAction->analyseNeutronFlux(fParticleEnergy, thisTrackID, radiusL, cosAngle, fParticleName);
      }
    }
  }
  //G4cout << " Exiting ProcessStepping. " << G4endl;
}
