#include "G4TARCSteppingAction.hh"

G4TARCSteppingAction::G4TARCSteppingAction(G4TARCEventAction* anEvent)
: G4UserSteppingAction(), fEventAction(anEvent)// , fHisto(0)
{

  //fEventAction = anEvent;
  //fHisto = G4TARCHistoManager::GetPointer();

  fRefShellThickness = 2.0 * mm;
  fRefShellOuterRad = 456.0 * mm;
  fRefShellInnerRad = fRefShellOuterRad - fRefShellThickness;
  fMaxOuterRadiusofShell = 1500.0 * mm;
  fMinInnerRadiusofShell = 10.0 * mm;
  fEnergy0 = 0.0;
  fNumberGenerations = 0;
  fShellNumber = (G4int)((fMaxOuterRadiusofShell - fMinInnerRadiusofShell) / fShellThickness + 0.5);

  // G4TARCRun* thisRA = static_cast<G4TARCRun*> (G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  // thisRA->ReadExperimentalDataFromFile(fExptlDataFileName);

  //ReadExperimentalDataFromFile(fExptlDataFileName);

}

void G4TARCSteppingAction::UserSteppingAction(const G4Step* myStep){
  //fHisto->ProcessStepping(myStep);
  ProcessStepping(myStep);
}


// User Stepping Action
void G4TARCSteppingAction::ProcessStepping(const G4Step* myStep){
    //G4cout << ".................. Stepping started" << G4endl;

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


  } else {
    // G4cout << " Not entering AnalyseSecondary conditionally" << G4endl;
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

  if (fParticleName == "neutron"){
    G4String PreVol = thePreTouchable->GetVolume()->GetName();
    G4String PostVol = thePostTouchable->GetVolume()->GetName();
    G4cout << fParticleName << "       Prevol-> " << PreVol <<  "    PostVol " << PostVol << G4endl;
    if (myStep->GetTrack()->GetNextVolume()){
      if (PreVol == "lab_phys" && PostVol == "world_log_PV"){
        fEventAction->exitingTally(true, fParticleEnergy);
      } else if (PostVol == "world_log_PV"){
        fEventAction->exitingTallyCheck(true);
      }
    }

    G4bool pre_inside = false;
    G4bool post_inside = false;
    G4cout << "----> HERE" << G4endl;

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
  //G4cout << ".......................... Exiting ProcessStepping. " << G4endl;
}


void G4TARCSteppingAction::ReadExperimentalDataFromFile(G4String& exptFileName){
  fReadData = false;
  std::ifstream exptIN(exptFileName, std::ios::in);
  G4String lineIN;
  unsigned NCount = 0, restCount = 0, file0 = 0, iTableNum = 0;
  G4bool readPara = false;
  G4double v1, v2, v3, v4, v5, v6, v7, v8, v9, v10;
  G4bool isFlux = false;
  G4int is40 = 0;
  std::vector<G4double> tmpV1;
  std::vector<G4double> tmpV2;
  std::vector<G4double> tmpV3;
  fMaxFluxData = -999999.99;
  fMaxFluenceData = -999999.99;
  fMaxTestFluxData = -999999.99;
  //G4int fileCount = 0;
  while (getline(exptIN, lineIN)){
    lineIN = std::regex_replace(lineIN, std::regex("^ +| +$|( ) +"), "$1");
    if (lineIN.size() > 1 && lineIN.find("#", 0, 1) != std::string::npos){  // if the line starts with # sign
      std::size_t found1 = lineIN.find("Table");
      readPara = false;
      G4String tableNum = (lineIN.substr(found1 + 5, 3));
      tableNum = std::regex_replace(tableNum, std::regex("^ +| +$|( ) +"), "$1"); // stripping extra spaces
      iTableNum = std::atoi(tableNum);
      isFlux = ((std::find( fFluxTableList.begin(), fFluxTableList.end(), iTableNum) != fFluxTableList.end()) || (iTableNum == fMeanEnergyTable)) ;
      if (iTableNum == fMeanEnergyTable) ++is40;
      if (found1 != std::string::npos){
        file0 = (iTableNum == 0) ? 1 : 0;
      }
      lineIN="";
    } else if (lineIN.size() > 1 && lineIN.find("#", 0, 1) == std::string::npos){   // the line does not start with # symbol
      if (lineIN.find(";", 0, 1) != std::string::npos){
        NCount = atoi(lineIN.substr(1, lineIN.size()).c_str());
        fIFluxCountRef = (iTableNum == 40) ? NCount : 0;
        restCount = NCount;
        readPara = true;
        //if (!isFlux) fMaxFluenceData = (std::max(fMaxFluenceData, (signed)NCount));
        //if (isFlux && !fIFluxCountRef) fMaxFluxData = (std::max(fMaxFluxData, (signed)NCount));
        //fMaxTestFluxData = fIFluxCountRef;
        //std::cout << "Table->" << iTableNum << " Data-> " << NCount
        //<< " Flux: " << isFlux << " Fluence: " << !(isFlux) << std::endl;
        continue;
      }
      if ( file0 && readPara){
        std::stringstream ss (lineIN);
        ss >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10;
        /*
        fExptEnergyBin.push_back(v1); fExptEnergyBin.push_back(v2); fExptEnergyBin.push_back(v3); fExptEnergyBin.push_back(v4);
        fExptEnergyBin.push_back(v5); fExptEnergyBin.push_back(v6); fExptEnergyBin.push_back(v7); fExptEnergyBin.push_back(v8);
        fExptEnergyBin.push_back(v9); fExptEnergyBin.push_back(v10);
        //G4cout << "filecount->" << fileCount << "   " << fExptEnergyBin.size() << G4endl;
        //for (unsigned ijk = 0 ; ijk < fExptEnergyBin.size(); ijk++) std::cout << fExptEnergyBin[ijk] << "  ";
        //std::cout << std::endl;
        */
        file0 = 0;
        readPara = false;
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

          if (isFlux && is40 == 2){
            //fMeanEnergyT40List.push_back(v1);  fMeanEnergyT40List.push_back(v2);  fMeanEnergyT40List.push_back(v3);
          } else {
            tmpV1.push_back(v1);  tmpV2.push_back(v2);  tmpV3.push_back(v3);
          }
          --restCount;
        } else if (wcount == 6) {
          ss >> v1 >> v2 >> v3 >> v4 >> v5 >> v6;
          if (isFlux && is40 == 2){
            //fMeanEnergyT40List.push_back(v1);  fMeanEnergyT40List.push_back(v2);  fMeanEnergyT40List.push_back(v3);
            //fMeanEnergyT40List.push_back(v4);  fMeanEnergyT40List.push_back(v5);  fMeanEnergyT40List.push_back(v6);
          } else {
            tmpV1.push_back(v1);  tmpV2.push_back(v2);  tmpV3.push_back(v3);
            tmpV1.push_back(v4);  tmpV2.push_back(v5);  tmpV3.push_back(v6);
          }
          restCount -= 2;
        }else if (wcount == 9){
          ss >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9;
          if (isFlux && is40 == 2){
            //fMeanEnergyT40List.push_back(v1);  fMeanEnergyT40List.push_back(v2);  fMeanEnergyT40List.push_back(v3);
            //fMeanEnergyT40List.push_back(v4);  fMeanEnergyT40List.push_back(v5);  fMeanEnergyT40List.push_back(v6);
            //fMeanEnergyT40List.push_back(v7);  fMeanEnergyT40List.push_back(v8);  fMeanEnergyT40List.push_back(v8);
            //fMeanEnergyT40List.push_back(v9);
          } else {
            tmpV1.push_back(v1);  tmpV2.push_back(v2);  tmpV3.push_back(v3);
            tmpV1.push_back(v4);  tmpV2.push_back(v5);  tmpV3.push_back(v6);
            tmpV1.push_back(v7);  tmpV2.push_back(v8);  tmpV3.push_back(v9);
          }
          restCount -= 3;
        }
      }
      if (!isFlux && (tmpV1.size()) == NCount){
        //if (iTableNum!=0) ++fMaxFluenceTable;
        fExptRadiiTables.push_back(tmpV1);
        //fExptFluenceTables.push_back(tmpV2);
        //fExptErrTables.push_back(tmpV3);
        std::vector<G4double>().swap(tmpV1);
        //std::vector<G4double>().swap(tmpV2);
        //std::vector<G4double>().swap(tmpV3);
      } else if (isFlux && (tmpV1.size()) == NCount  && is40 !=2) {
        //fExptEnergyTables.push_back(tmpV1);
        //fExptFluxTables.push_back(tmpV2);
        //fExptFluxErrTables.push_back(tmpV3);
        //std::vector<G4double>().swap(tmpV1);
        //std::vector<G4double>().swap(tmpV2);
        //std::vector<G4double>().swap(tmpV3);
      }

    }
    lineIN="";
  }
  exptIN.close();
  std::vector<G4double>().swap(tmpV1);
  std::vector<G4double>().swap(tmpV2);
  std::vector<G4double>().swap(tmpV3);

  // Now sorting the vector fMeanEnergyT40List
  //std::sort(fMeanEnergyT40List.begin(), fMeanEnergyT40List.end(),  [] (G4double const& a, G4double const& b) { return a < b; });

  fMaxRadCount = fExptRadiiTables[8].size();
  //  G4cout << " Rad count = " << fMaxRadCount << G4endl;
  //  fMaxRadCount = 10  fMaxTestFluxData = 21;
  //  fMaxFluxData  = 95 fMaxFluenceData  = 102

  for (std::size_t i = 0; i < fExptRadiiTables.size(); i++){
    for (std::size_t j = 0; j < fExptRadiiTables[i].size(); j++) {
      fExptRadiiTables[i][j] *= 10.0;  // in mm now
    }
  }

}
