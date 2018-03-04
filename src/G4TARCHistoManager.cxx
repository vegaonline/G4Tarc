#include "G4TARCHistoManager.hh"

G4TARCHistoManager *  G4TARCHistoManager::fHistoManager = 0;

G4TARCHistoManager *  G4TARCHistoManager::GetPointer() {
  if ( !fHistoManager ) {
    static G4TARCHistoManager manager;
    fHistoManager = &manager;
  }
  return fHistoManager;
}

G4TARCHistoManager::G4TARCHistoManager()
:
  fHisto( 0 ),
  fHistoBooked(false),
  fPrimaryDef( 0 ),
  fProton( 0 ),
  fNHisto( 25 ),
  fNeutron( 0 ),
  fEdepMax( fMaxEVal ),
  fLength( fMaxLVal ),
  fPrimaryKineticEnergy ( fEVal0 ),
  fVerbose( 0 ),
  fNBinsE( 100),  //fMaxBin ),
  fNSlices( fMaxSlices ) {
    //fDetector = new G4TARCDetectorConstruction();
    fHisto    = new G4TARCHisto();
    fHisto->SetVerbose(fVerbose);
    fNeutron = G4Neutron::Neutron();
    fProton = G4Proton::Proton();
}

G4TARCHistoManager::~G4TARCHistoManager() {
  //delete fHisto;
  //delete G4AnalysisManager::Instance();
}

void G4TARCHistoManager::BookHisto() {
  fHistoBooked = true;
  //fHisto->Add1D("1","Energy deposition (MeV/mm/event) in the target", fNSlices, 0.0,fLength/mm,MeV/mm);
  fHisto->Add1D("1","Energy deposition (keV) in the target", fNBinsE, 0.0,1200*MeV, 1.0);
  /*
  fHisto->Add1D("2","Log10 Energy (MeV) of gammas",                   fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("3","Log10 Energy (MeV) of electrons",                fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("4","Log10 Energy (MeV) of positrons",                fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("5","Log10 Energy (MeV) of protons",                  fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("6","Log10 Energy (MeV) of neutrons",                 fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("7","Log10 Energy (MeV) of charged pions",            fNBinsE, -4.0,6.0,       1.0);
  fHisto->Add1D("8","Log10 Energy (MeV) of pi0",                      fNBinsE, -4.0,6.0,       1.0);
  fHisto->Add1D("9","Log10 Energy (MeV) of charged kaons",            fNBinsE, -4.0,6.0,       1.0);
  fHisto->Add1D("10","Log10 Energy (MeV) of neutral kaons",           fNBinsE, -4.0,6.0,       1.0);
  fHisto->Add1D("11","Log10 Energy (MeV) of deuterons and tritons",   fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("12","Log10 Energy (MeV) of He3 and alpha",           fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("13","Log10 Energy (MeV) of Generic Ions",            fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("14","Log10 Energy (MeV) of muons",                   fNBinsE, -4.0,6.0,       1.0);
  fHisto->Add1D("15","log10 Energy (MeV) of side-leaked neutrons",    fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("16","log10 Energy (MeV) of forward-leaked neutrons", fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("17","log10 Energy (MeV) of backward-leaked neutrons",fNBinsE, -5.0,5.0,       1.0);
  fHisto->Add1D("18","log10 Energy (MeV) of leaking protons",         fNBinsE, -4.0,6.0,       1.0);
  fHisto->Add1D("19","log10 Energy (MeV) of leaking charged pions",   fNBinsE, -4.0,6.0,       1.0);
  fHisto->Add1D("20","Log10 Energy (MeV) of pi+",                     fNBinsE, -4.0,6.0,       1.0);
  fHisto->Add1D("21","Log10 Energy (MeV) of pi-",                     fNBinsE, -4.0,6.0,       1.0);
  fHisto->Add1D("22","Energy deposition in the target normalized to beam energy",110,0.0,1.1,1.0);
  fHisto->Add1D("23","EM energy deposition in the target normalized to beam energy",110,0.0,1.1,1.0);
  fHisto->Add1D("24","Pion energy deposition in the target normalized to beam energy",110,0.0,1.1,1.0);
  fHisto->Add1D("25","Proton energy deposition in the target normalized to beam energy",110,0.0,1.1,1.0);
  */
}


void G4TARCHistoManager::BeginOfRun() {
  fAbsX0 = fAbsY0 = fAbsZ0 = 0.5 * fLength;
  fNevt       = 0;
  fNelec      = 0;
  fNposit     = 0;
  fNgam       = 0;
  fNstep      = 0;
  fNprot_leak = 0;
  fNmuons     = 0;
  fNions      = 0;
  fNdeut      = 0;
  fNalpha     = 0;
  fNneutron   = 0;
  fNproton    = 0;
  fNaproton   = 0;
  fNneu_forw  = 0;
  fNneu_leak  = 0;
  fNneu_back  = 0;
  fNinelastic = 0;
  fNelastic   = 0;
  fNzero      = 0;
  fEdepSum    = 0.0;
  fEdepSum2   = 0.0;
  fRange      = fMaxLVal;
  fRho        = 60.0 * mm;  //-------------------------->  Diagnose with GDML file  ********
  fLMax       = fMaxLVal;
  fLBin       = fRange / fLMax;
  fTcut       = 100.0 * CLHEP::MeV;
  fNeutronInit = fNeutronSum = fNeutronBreed = 0.0;
  fEsecond    = G4DataVector(fMaxBin, 0.0); //  fNumMax, 0.);
  fMsecond    = G4DataVector(fMaxBin, 0.0); //  fNumMax, 0.);
  // n-spectra, E, T and ET
  fTmin       = 0.001 * CLHEP::eV;
  fTmax       = fMaxEVal;   //  ***** I want to give value ~ 2.0 * incident beam energy ******* // CHECK
  // fTmax    = 1000. * CLHEP::MeV;
  fTimeMin   = 1. * nanosecond;
  fTimeMax   = 1. * millisecond;
  fNbin      = fMaxBin; // fNumMax; // 60; // 100;  // 1000; // 10; //
  fnEsecond  = G4PhysicsLogVector(fTmin,fTmax,fNbin);
  fnTsecond  = G4PhysicsLogVector(fTimeMin,fTimeMax,fNbin);
  fNSecondSum1 = G4DataVector(fNbin, 0.0);
  fNSecondSum2 = G4DataVector(fNbin, 0.0);
  fNSecondSum3 = G4DataVector(fNbin, 0.0);
  fnETsum    = G4DataVector(fNbin * fNbin, 0.);

  for( G4int ii = 0; ii < fNbin; ii++ )
  {
    std::vector<G4double> temp;
    for( G4int jj = 0; jj < fNbin; jj++ ) temp.push_back(0.0);
    fET.push_back(temp);
  }
  fGunParticleX   = G4DataVector(fLMax, 0.0);
  fGunParticleY   = G4DataVector(fLMax, 0.0);
  fGunParticleZ   = G4DataVector(fLMax, 0.0);
  fGunParticleTLX = G4DataVector(fLMax, 0.0);
  fGunParticleTLY = G4DataVector(fLMax, 0.0);
  fGunParticleTLZ = G4DataVector(fLMax, 0.0);
  fGunParticleDep = G4DataVector(fLMax, 0.0);
  fGunParticleRho = G4DataVector(fLMax, 0.0);
  fRangeVector    = G4DataVector(fLMax, 0.0);
  fRhoVector      = G4DataVector(fLMax, 0.0);
  fDeltaVector    = G4DataVector(fLMax, 0.0);
  for( G4int ii = 0; ii < fLMax; ii++ ) {
    fRangeVector[ii] = ii * fLBin - 0.5 * fRange; //
    fRhoVector[ii]   = ii * fRho / fLBin;
    fDeltaVector[ii] = ii * fTmax / fLBin;
  }
  fTkin = fEVal0;
  fEbin = fTkin / fNumMax;
  G4double EDelta = (fMaxEVal / fMaxBin);
  G4double MDelta = (fEVal0 / fMaxBin);
  for( G4int ii = 0; ii < fMaxBin; ii++ )
  {
    fEsecond[ii]  = ii * EDelta;
    // fnEsecond[ii] = ii * fEbin;
    fMsecond[ii]  = ii * MDelta;
  }

  if (!fHistoBooked) BookHisto();

  fHisto->Book();
  if( fVerbose > 0 )
    G4cout << "G4TARCHistoManager: Histograms are booked and run has been started"
           <<G4endl;
}


void G4TARCHistoManager::EndOfRun() {
  G4cout << "G4TARCHistoManager ; End of run actions are started" << G4endl;
  G4cout << "fNevt = " << fNevt << G4endl;
  G4cout << "EndOfRun(), fEdepSum = " << fEdepSum << G4endl;
  G4cout << "======================================================================" << G4endl;
  G4double x = ( G4double )fNevt, perN = 1.0;
  if (fNevt > 0){
    perN = 1.0 / x;
    x = perN;
  }
  TrackRun(x);  // track and get leaks
  NeutronRun(x); // neutron leak data and spectra
  GunParticleRun(x);  // beam distribution in target.

  //****************** This is for G4ParticleHPThermalScattering
  fNHisto = 1;

  for (G4int i = 0; i < fNHisto; i++)
    fHisto->ScaleH1(i, x);           // Normalize Histogram

  fHisto->Save();
}

void G4TARCHistoManager::BeginOfEvent() {
  fEdepEvt    = 0.0;
  fEdepEM     = 0.0;
  fEdepProton = 0.0;
  fNsecondary = 0.0;
  fRangeSum   = G4ThreeVector(0.0, 0.0, 0.0);
  fStepSum    = 0.0;
  fDeltaSum   = 0.0;
}

void G4TARCHistoManager::EndOfEvent() {
  fEdepSum  += fEdepEvt;
  fEdepSum2 += fEdepEvt * fEdepEvt;
  fNevt ++;
  if (fNsecondary > 0) fNinelastic++;
  //fHisto->Fill(21,fEdepEvt/fPrimaryKineticEnergy,1.0);
  //fHisto->Fill(22,fEdepEM/fPrimaryKineticEnergy,1.0);
  //fHisto->Fill(23,fEdepPI/fPrimaryKineticEnergy,1.0);
  //fHisto->Fill(24,fEdepProton/fPrimaryKineticEnergy,1.0);

  WriteEventRange( fRangeSum, fStepSum, fDeltaSum );
}


// Normalization for no hadron interaction in the target. EM and hadron elastic shuld be inactivated
void G4TARCHistoManager::AddNzero(const G4Track* myTrack, const G4Step* myStep ) {
  G4double cosTheta = myTrack->GetDynamicParticle()->GetMomentumDirection().x();
  if (cosTheta > 0.999) {
    fNzero++;
  }else {
    fNelastic++;
  }
}


void G4TARCHistoManager::GunParticleDistribution( const G4Track* myTrack, const G4Step* myStep) {
  G4double length = myStep->GetStepLength();
  G4double dEdx   = myStep->GetTotalEnergyDeposit();

  G4ThreeVector pGun1 = myStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector pGun2 = myStep->GetPostStepPoint()->GetPosition();

  fRangeSum += (pGun2 - pGun1);
  fStepSum += length;
  fDeltaSum += dEdx;
}

void G4TARCHistoManager::AddTargetStep(const G4Step* myStep) {
  fNstep++;

  G4double fEdep = myStep->GetTotalEnergyDeposit();
  if (fVerbose > 1) {
    G4cout << "TargetSD::ProcessHists: beta1 =" << myStep->GetPreStepPoint()->GetVelocity() / c_light
           << " beta2 = "                       << myStep->GetPostStepPoint()->GetVelocity() / c_light
           << " weight = "                      <<  myStep->GetTrack()->GetWeight()
           << G4endl;
  }
  //if (fEdep >= DBL_MIN) {
    const G4Track* myTrack = myStep->GetTrack();
    G4ThreeVector pos = (myStep->GetPreStepPoint()->GetPosition() + myStep->GetPostStepPoint()->GetPosition())*0.5;

    G4double x = pos.x() - fAbsX0;
    G4double y = pos.y() - fAbsY0;
    G4double z = pos.z() - fAbsZ0;

    //Scoring.
    fEdepEvt += fEdep;
    G4cout << " Energy for histogram-> " << fEdep/keV << " keV" << G4endl;
    fHisto->Fill(0, z, fEdep/keV);
    const G4ParticleDefinition* pd = myTrack->GetDefinition();

  if (pd == G4Gamma::Gamma() || pd == G4Electron::Electron() || pd == G4Positron::Positron()) {
    fEdepEM += fEdep;
  } else if (pd == G4Proton::Proton() || pd == G4AntiProton::AntiProton()){
    fEdepProton +=fEdep;
  } else if (pd == G4PionPlus::PionPlus() || pd == G4PionMinus::PionMinus()){
    fEdepPI +=fEdep;
  }

  if (fVerbose > 1) G4cout << "HistoManager::AddEnergy: E(keV) = " << fEdep / keV
                           << "; z(mm) = " << z / mm
                           << "; step(mm) = " << myStep->GetStepLength() / mm
                           << " by " << pd->GetParticleName()
                           << " E(MeV) = " << myTrack->GetKineticEnergy()/MeV
                           << G4endl;
  //}
}


void G4TARCHistoManager::ScoreNewTrack( const G4Track* myTrack) {
  const G4ParticleDefinition* pd = myTrack->GetDefinition();
  G4String name = pd->GetParticleName();
  G4double ke = myTrack->GetKineticEnergy();
  G4double TKin = myTrack->GetKineticEnergy();

// For Primary Track
  if (myTrack->GetTrackID() == 0) {
    fNevt++;
    fPrimaryKineticEnergy = ke;
    fPrimaryDef = pd;
    G4ThreeVector dir = myTrack->GetMomentumDirection();

    if (fVerbose > 1) G4cout << "### Primary "   << name << " KineticE(MeV) = " << ke / MeV
                             << "; Mass(MeV) = " << pd->GetPDGMass() / MeV
                             << "; Pos(mm) = "   << myTrack->GetPosition() / mm
                             << "; Dir = "       << myTrack->GetMomentumDirection()
                             << G4endl;
  } else {
// For Secondary tracks.
    if (fVerbose > 1) G4cout << "=== Secondary " << name << " KineticE(MeV) = " << ke / MeV
                             << "; Mass(MeV) = " << pd->GetPDGMass() / MeV
                             << "; Pos(mm) = "   << myTrack->GetPosition() / mm
                             << "; Dir = "       << myTrack->GetMomentumDirection()
                             << G4endl;
    ke = std::log10(ke/MeV);
    if (pd == G4Gamma::Gamma()){
      fNgam++;
      //fHisto->Fill(1, ke, 1.0);
    }else if (pd == G4Electron::Electron()){
      fNelec++;
      //fHisto->Fill(2, ke, 1.0);
    }else if (pd == G4Positron::Positron()) {
      fNposit++;
      //fHisto->Fill(3, ke, 1.0);
    } else if (pd == G4Proton::Proton()){
      fNproton++;
      //fHisto->Fill(4, ke, 1.0);
    } else if (pd == fNeutron && TKin < fTcut){  // <----- CHECK
      fNneutron++;
      //fHisto->Fill(5, ke, 1.0);
    } else if (pd == G4AntiProton::AntiProton()){
      fNaproton++;
    } else if ( pd == G4PionPlus::PionPlus() ) {
      fNcpions++;
      //fHisto->Fill(6, ke, 1.0);
      //fHisto->Fill(19, ke, 1.0);
    } else if ( pd == G4PionMinus::PionMinus()) {
      fNcpions++;
      //fHisto->Fill(6, ke, 1.0);
      //fHisto->Fill(20, ke, 1.0);
    } else if ( pd == G4PionZero::PionZero()) {
      fNpi0++;
      //fHisto->Fill(7, ke, 1.0);
    } else if ( pd == G4KaonPlus::KaonPlus() ||
                pd == G4KaonMinus::KaonMinus()) {
      fNkaons++;
      //fHisto->Fill(8, ke, 1.0);
    } else if ( pd == G4KaonZeroShort::KaonZeroShort() ||
                pd == G4KaonZeroLong::KaonZeroLong()) {
      fNkaons++;
      //fHisto->Fill(9, ke, 1.0);
    } else if (pd == G4Deuteron::Deuteron() || pd == G4Triton::Triton()){
      fNdeut++;
      //fHisto->Fill(10, ke, 1.0);
    } else if (pd == G4He3::He3() || pd == G4Alpha::Alpha()) {
      fNalpha++;
      //fHisto->Fill(11, ke, 1.0);
    } else if (pd->GetParticleType() == "nucleus"){
      fNions++;
      //fHisto->Fill(12, ke, 1.0);
    } else if (pd == G4MuonPlus::MuonPlus() || pd == G4MuonMinus::MuonMinus()){
      fNmuons++;
      //fHisto->Fill(13, ke, 1.0);
    }
  }
}

void G4TARCHistoManager::AddLeakingParticle(const G4Track* myTrack) {
  const G4ParticleDefinition* pd           = myTrack->GetDefinition();
  const G4DynamicParticle*    dynSecondary = myTrack->GetDynamicParticle();
  G4double                    Tkin         = dynSecondary->GetKineticEnergy();
/*
  if (myTrack->GetKineticEnergy()/MeV <= 0.0) {
    G4String pname  = pd->GetParticleName();
    G4double mass   = pd->GetPDGMass();
    return;
  }
  */
  G4double en       = std::log10(myTrack->GetKineticEnergy()/MeV);
  G4ThreeVector pos = myTrack->GetPosition();
  G4ThreeVector dir = myTrack->GetMomentumDirection();

  G4double xx       = pos.x();
  G4double yy       = pos.y();
  G4double zz       = pos.z();
  G4bool isLeaking  = false;

  if ( (
         (xx > -fAbsX0 && dir.x() > 0.0) // - was added later originally > +fAbsX0
    ||   (yy > -fAbsY0 && dir.y() > 0.0)
    ||   (zz > -fAbsZ0 && dir.z() > 0.0)
       )
  ) { // forward
      isLeaking = true;
      if (pd == fNeutron && Tkin < fTcut) {
        ++fNneu_forw;
        //fHisto->Fill(15, en,  1.0);
      }
  }

  if ( (
         (xx < fAbsX0 && dir.x() < 0.0)
    ||   (yy < fAbsY0 && dir.y() < 0.0)
    ||   (zz < fAbsZ0 && dir.z() < 0.0)
       )
  ) { // backward
    isLeaking = true;
    if (pd == fNeutron && Tkin < fTcut)   {
       ++fNneu_back;
       //fHisto->Fill(16, en,  1.0);
     }
  }

  if ( (
         (std::abs(xx) <= fAbsX0 && (zz * dir.z()  + yy * dir.y()) > 0.0 )
    ||   (std::abs(yy) <= fAbsY0 && (xx * dir.x()  + zz * dir.z()) > 0.0 )
    ||   (std::abs(zz) <= fAbsZ0 && (xx * dir.x()  + yy * dir.y()) > 0.0 )
       )
  ) { // side
    isLeaking = true;
    if (pd == fNeutron && Tkin < fTcut)  {
      ++fNneu_leak;
      //fHisto->Fill(14, en,  1.0);
    }
  }

  if (isLeaking) {
    if (pd == G4Proton::Proton() && myTrack->GetTrackID() == 1){
       ++fNprot_leak;
       //fHisto->Fill(17, en,  1.0);
     }
    if (pd == G4PionPlus::PionPlus() || pd == G4PionMinus::PionMinus()) {
      ++fNPionleak;
      //fHisto->Fill(18, en,  1.0);
    }
  }
}


// pA->nX low energy data neutron final state
void G4TARCHistoManager::NeutFinalState(const G4Track* myTrack, const G4Step* myStep) {
  fNsecondary++;
  G4double TKin = 0.0, Tspectra = 0.0, stepLen, time;
  G4int ix;
  size_t ii, jj;

  if (myTrack->GetDynamicParticle()->GetDefinition()->GetParticleName() == "neutron"){
    fNeutronSum+= 1.0;
    // future plan for neutron breed at breeder with X
    TKin = myTrack->GetDynamicParticle()->GetKineticEnergy();
    stepLen = myStep->GetStepLength();
    time = myTrack->GetGlobalTime();
    for (ii = 0; ii < fNbin; ii++) {
      if (TKin <= fnEsecond.GetLowEdgeEnergy(ii)){
        jj = (ii == 0) ? ii : ii - 1;
        break;
      }
    }
    jj = (ii - fNbin) ? fNbin-1 : jj;

  }
}


void G4TARCHistoManager::TargetProfile(const G4Track* myTrack, const G4Step* myStep) {
  G4double TKinE = 0.0, TSpectr = 0.0, stepLen, time, xn;
  G4int ix = 0;

  if (myTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "neutron") {
    if (
           (myTrack->GetLogicalVolumeAtVertex()->GetName() == "blockB_log")
           && (  myTrack->GetVertexPosition().x() >= -35.0 && myTrack->GetVertexPosition().x() <= 35.0
           || myTrack->GetVertexPosition().y() >= -35.0 && myTrack->GetVertexPosition().y() <= 35.0
           || myTrack->GetVertexPosition().z() >= -1500.0 && myTrack->GetVertexPosition().z() <= -299.0
           )
           && (myTrack->GetParentID() == 1)
     ) {
       fNeutronInit += 1.0;
       xn = myTrack->GetVertexPosition().x() + 0.5 * fRange;
       ix    = G4int(xn / fLBin + 0.5);
       if (ix >= 0 && ix < fLMax) fGunParticleX[ix] += 1.0;

       xn = myTrack->GetVertexPosition().y() + 0.5 * fRange;
       ix    = G4int(xn / fLBin + 0.5);
       if (ix >= 0 && ix < fLMax) fGunParticleY[ix] += 1.0;

       xn = myTrack->GetVertexPosition().z() + 0.5 * fRange;
       ix    = G4int(xn / fLBin + 0.5);
       if (ix >= 0 && ix < fLMax) fGunParticleZ[ix] += 1.0;
     }
  }
}


void G4TARCHistoManager::AddEnergyTime(const G4Track* myTrack, const G4Step* myStep) {
  G4double Tkin = 0.0, myTime, myStepLength;
  size_t ii, jj, eii, tjj;

  if (myTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "neutron") {
    Tkin = myTrack->GetDynamicParticle()->GetKineticEnergy();
    myStepLength = myStep->GetStepLength();
    myTime = myTrack->GetProperTime();

    for (ii = 0; ii < fNbin; ii++) {
      if (Tkin <= fnEsecond.GetLowEdgeEnergy(ii)) {
        eii = (ii == 0) ? ii : ii - 1;
        break;
      }
    }
    if (ii == fNbin) eii = fNbin -1 ;
    for (jj = 0; jj < fNbin; jj++) {
      if (myTime <= fnTsecond.GetLowEdgeEnergy(jj)) {
        tjj = (jj == 0) ?  jj : jj - 1;
        break;
      }
    }
    if (jj == fNbin) tjj = fNbin - 1;
    fET[tjj][eii] += 1.0;
  }
}

void G4TARCHistoManager::WriteEventRange(G4ThreeVector rsum, G4double lsum, G4double dsum) {
  G4int ix;
  G4double x, y, z, rho, range, xbin, ybin, zbin, rbin, dbin, dEdx, length;

  x = rsum.x();
  y = rsum.y();
  z = rsum.z();

  rho = std::sqrt(x * x + y * y);

  xbin = ybin = zbin = fRange / fMaxBin;   // 100 can be changed checking the binning
  rbin               = fRho / fMaxBin;
  dbin               = fTmax  / fMaxBin;

  range = x / (xbin) + 0.5;
  ix = G4int(range);
  if (ix > 0 && ix < fMaxBin) fGunParticleX[ix] += 1.0;

  range = y / (ybin) + 0.5;
  ix = G4int(range);
  if (ix > 0 && ix < fMaxBin) fGunParticleY[ix] += 1.0;

  range = z / (zbin) + 0.5;
  ix = G4int(range);
  if (ix > 0 && ix < fMaxBin) fGunParticleZ[ix] += 1.0;

  range = lsum / xbin + 0.5;
  ix = (G4int)(range);
  if (ix > 0 && ix < fMaxBin) fGunParticleTLX[ix] += 1.0;

  range = lsum / ybin + 0.5;
  ix = (G4int)(range);
  if (ix > 0 && ix < fMaxBin) fGunParticleTLY[ix] += 1.0;

  range = lsum / zbin + 0.5;
  ix = (G4int)(range);
  if (ix > 0 && ix < fMaxBin) fGunParticleTLZ[ix] += 1.0;

  range = dsum / dbin + 0.5;
  ix = (G4int)(range);
  if (ix > 0 && ix < fMaxBin) fGunParticleDep[ix] += 1.0;

  range = rho / rbin + 0.5;
  ix = (G4int)(range);
  if (ix > 0 && ix < fMaxBin) fGunParticleRho[ix] += 1.0;
}

void G4TARCHistoManager::TrackRun(G4double x) {
  std::ofstream trackout("track.dat", std::ios::out);
  G4double xStep         = x * (G4double)fNstep;
  G4double xElectron     = x * (G4double)fNelec;
  G4double xPositron     = x * (G4double)fNposit;
  G4double xGamma        = x * (G4double)fNgam;
  G4double xProton       = x * (G4double)fNproton;
  G4double xAntiProton   = x * (G4double)fNaproton;
  G4double xProtonLeak   = x * (G4double)fNprot_leak;
  G4double xNeutron      = x * (G4double)fNneutron;
  G4double xNeutronLeak  = x * (G4double)fNneu_leak;
  G4double xneuF         = x * (G4double)fNneu_forw;
  G4double xneuB         = x * (G4double)fNneu_back;
  G4double xMuons        = x * (G4double)fNmuons;
  G4double xIons         = x * (G4double)fNions;
  G4double xAlpha        = x * (G4double)fNalpha;
  G4double xDeut         = x * (G4double)fNdeut;

  trackout << " Total Eenergy Deposited (MeV) = " << fEdepSum/MeV << G4endl;
  fEdepSum  *= x;
  fEdepSum2 *= x;
  fEdepSum2 -= fEdepSum * fEdepSum;
  fEdepSum2  = (fEdepSum2 > 0.0) ? std::sqrt(fEdepSum2) : 0.0;


  trackout << " x = 1/fNevt = "               << x                                          << G4endl;
  trackout << "x * Total Energy deposited (MeV) = " << fEdepSum/MeV                         << G4endl;
  if (fPrimaryDef){
    trackout << "Beam Particle "                << fPrimaryDef->GetParticleName()           << G4endl;
    trackout << "Beam Energy   "                << fPrimaryKineticEnergy/MeV      << " MeV" << G4endl;
  }
  trackout << "Number of Events "             << fNevt                                    << G4endl;
  trackout <<                                                                                G4endl;

  trackout << " Production in Target:"                                                                                                << G4endl;
  trackout << std::setprecision(4) << "Mean Energy deposit  "         << fEdepSum/MeV << " MeV " << " +- " << fEdepSum2/MeV << " MeV" << G4endl;
  trackout << std::setprecision(4) << "Average Number of steps "      << xStep                                                        << G4endl;
  trackout << std::setprecision(4) << "Average Number of Gammas "     << xGamma                                                       << G4endl;
  trackout << std::setprecision(4) << "Average Number of Electrons "  << xElectron                                                    << G4endl;
  trackout << std::setprecision(4) << "Average Number of Positrons "  << xPositron                                                    << G4endl;
  trackout << std::setprecision(4) << "Average Number of Protons "    << xProton                                                      << G4endl;
  trackout << std::setprecision(4) << "Average Number of AntiProton " << xAntiProton                                                  << G4endl;
  trackout << std::setprecision(4) << "Average Number of Neutrons "   << xNeutron                                                     << G4endl;
  trackout << std::setprecision(4) << "Average Number of Muons "      << xMuons                                                       << G4endl;
  trackout <<                                                                                                                            G4endl;

  trackout << " leakage from the system: "                                                                                            << G4endl;
  trackout << std::setprecision(4) << "Average Number of Leaked Neutrons " << xNeutronLeak << G4endl;
  trackout << std::setprecision(4) << "Average Number of Leaked Protons "  << xProtonLeak  << G4endl;
  trackout <<                                                                                                                            G4endl;

  G4double kEffective, rho, rat, react, perN=x;
  kEffective = fNeutronSum / fNeutronInit;
  rho        = (kEffective - 1.0) / kEffective;  // reactivity :: deviation from criticality
  rat        = std::log(kEffective);
  react      = rat / (1.0 + rat);

  trackout << " IMP Parameters : "    << G4endl;
  trackout << " Neutron_Init/p = "    << fNeutronInit* perN << ",  neutron_Sum/p = " << fNeutronSum * perN << G4endl;
  trackout << " kEffective = "        << kEffective         << ", Rho = "            << rho                << G4endl;
  trackout << "Estimated reactivity " << react                                       << G4endl             << G4endl;
  trackout << "==========================================================================================" << G4endl;
  trackout << G4endl;

}


void G4TARCHistoManager::NeutronRun(G4double x) {
  G4double perN = x, ngSum = 0.0;
  std::ofstream nspec("neutronSpectra.dat", std::ios::out);

  nspec << fNSecondSum1.size() << G4endl;

  for (size_t k = 0; k < fNSecondSum1.size(); k++){
      nspec << fnEsecond.GetLowEdgeEnergy(k)/MeV
            << "  " << perN * fNSecondSum1[k] << "   " << perN * fNSecondSum2[k] << "  " << perN * fNSecondSum3[k] << G4endl;
      ngSum += fNSecondSum1[k];
  }
}


void G4TARCHistoManager::GunParticleRun(G4double x) {
  std::ofstream gunspectrum("gun.dat", std::ios::out);

  gunspectrum << fRangeVector.size() << G4endl;

  for (size_t k = 0; k < fRangeVector.size(); k++) {
    gunspectrum << fRangeVector[k] << "  " << x * fGunParticleX[k]  << "  " << x * fGunParticleY[k]   << "  "   << x * fGunParticleZ[k]
                                   << "  " << x * fGunParticleTLX[k]
                                   << "  " << x * fGunParticleTLY[k]
                                   << "  " << x * fGunParticleTLZ[k]
                                   << "  " << fDeltaVector[k]        << "  "   << x * fGunParticleDep[k]
                                   << "  " << fRhoVector[k]       << "  " << x * fGunParticleRho[k] << G4endl;
  }
}

void G4TARCHistoManager::Fill(G4int id, G4double x, G4double w) {
  fHisto->Fill(id, x, w);
}


  ////////////////////////////////////////////////////////////////////
/*
  void G4TARCHistoManager::SetVerbose(G4int val)
  {
    fVerbose = val;
    fHisto->SetVerbose(val);
  }
*/
