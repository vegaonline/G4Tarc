#include "G4TARCHisto.hh"

G4TARCHisto::~G4TARCHisto() {
  std::vector<G4double>().swap(fExptEnergyBin);
  std::vector<G4double>().swap(fInnerRadiusofShell);
  std::vector<G4double>().swap(fOuterRadiusofShell);
}


void G4TARCHisto::SetGeomParam(G4int& fShellN, G4double& fShellThk,
  std::vector<G4double>& inrad, std::vector<G4double>& outrad) {
  fInnerRadiusofShell = inrad;
  fOuterRadiusofShell = outrad;
  fShellNumber       = fShellN;
  fShellThickness    = fShellThk;
  fMinInnerRadiusofShell  = fInnerRadiusofShell[0];
}

void G4TARCHisto::SetExptDataParam( G4int& maxFl, G4int& maxFlu, G4int& maxTstFl,
  std::vector<G4double>& energyBinTab, std::vector<std::vector<G4double> >& energyTab,
  std::vector<std::vector<G4double> >& fluxTab, std::vector<std::vector<G4double> >& fluxErrTab){
  fMaxlowFluxData   = maxFl;
  fMaxFluxData      = maxFl;
  fMaxFluenceData   = maxFlu;
  fMaxTestFluxData  = maxTstFl;
  fExptEnergyBin     = energyBinTab;
  fExptEnergyTables  = energyTab;
  fExptFluxTables    = fluxTab;
  fExptFluxErrTables = fluxErrTab;
  fMaxFluenceSpectrumData = 10 * fMaxFluenceData;
}


void G4TARCHisto::InitializeParams(){
  fN_max                       = 0;
  fFractional_bin_width        = 0.2;
  fTotal_flux                  = 0.0;
  fEnergy_integral             = G4DataVector(fMaxDimEnInt, 0.0);
  fNeutFlux                    = G4DataVector(fMaxDimEnInt, 0.0);
  fEnFlux                      = G4DataVector(fMaxDimEnInt, 0.0);
  fFluence_spectrum            = G4DataVector(fMaxFluenceData, 0.0);
  fRadial_energies             = G4DataVector(fMaxRadIndex, 0.0);
  fRadii                       = G4DataVector(fMaxRadIndex, 0.0);
  fFlux_radius.resize(fMaxRadIndex, std::vector<G4double>(fMaxRadIndex));
  fRadial_fluence_step.resize(fMaxFluxData, std::vector<G4double>(fMaxRadIndex));
  fFine_energy                 = G4DataVector(fMaxFluxData, 0.0);
  fFlux                        = G4DataVector(fMaxTestFluxData, 0.0);
  fFlux_data                   = G4DataVector(fMaxTestFluxData, 0.0);
  fCos_flux                    = G4DataVector(fMaxTestFluxData, 0.0);
  fFluence                     = G4DataVector(fMaxTestFluxData, 0.0);
  fFluence_spectrum            = G4DataVector(fMaxFluenceSpectrumData, 0.0);
  fFluence_step                = G4DataVector(fMaxTestFluxData, 0.0);
  fFluence_cyl                 = G4DataVector(fMaxTestFluxData, 0.0);
  fFluence_step_cyl            = G4DataVector(fMaxTestFluxData, 0.0);
  fFluence_step_shell          = G4DataVector(fMaxTestFluxData, 0.0);
  fFluence_front_step          = G4DataVector(fMaxTestFluxData, 0.0);
  fEFlux                       = G4DataVector(fMaxTestFluxData, 0.0);
  fFine_eFlux                  = G4DataVector(2 * fMaxTestFluxData, 0.0);
  fFlux_energy                 = G4DataVector(fMaxFluxData, 0.0);
  fFlux_stat_error             = G4DataVector(fMaxFluxData, 0.0);
  fFlux_systematic_error       = G4DataVector(fMaxFluxData, 0.0);
  fLowEnergy                   = G4DataVector(fMaxlowFluxData, 0.0);
  fLowFlux                     = G4DataVector(fMaxlowFluxData, 0.0);
  fLowFluxData                 = G4DataVector(fMaxlowFluxData, 0.0);
  fLow_stat                    = G4DataVector(fMaxlowFluxData, 0.0);
  fLow_syst                    = G4DataVector(fMaxlowFluxData, 0.0);
  fCos_low_flux                = G4DataVector(fMaxlowFluxData, 0.0);
  fLow_fluence                 = G4DataVector(fMaxlowFluxData, 0.0);
  fLow_fluence_step            = G4DataVector(fMaxlowFluxData, 0.0);
  fLow_fluence_front_step      = G4DataVector(fMaxlowFluxData, 0.0);
  fLow_fluence_cyl             = G4DataVector(fMaxlowFluxData, 0.0);
  fLow_fluence_step_cyl        = G4DataVector(fMaxlowFluxData, 0.0);
  fLow_fluence_step_shell      = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_energy              = G4DataVector(fMaxlowFluxData, 0.0) ;
  fLithium_radial_energy_lower = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_radial_energy_upper = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_radial_mean         = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_radial_fluence_step = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_flux                = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_flux_data           = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_stat                = G4DataVector(fMaxlowFluxData, 0.0);
  fCos_lithium_flux            = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_fluence             = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_fluence_Step        = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_fluence_front_Step  = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_fluence_cyl         = G4DataVector(fMaxlowFluxData, 0.0);
  fLithium_Zflux               = G4DataVector(fMaxlowFluxData, 0.0);
  fGamma_flux               = 0;
  fExiting_flux             = 0;
  fExiting_check_flux       = 0;
  fNeutron_fluence          = 0;
  fNeutron_fluence_cyl      = 0;
  fIntegral_flux_p          = 0;
  fIntegral_Eflux_p         = 0.0;
  fNeutron_flux             = 0;
  fNeutron_check            = 0;
  fElectron_flux            = 0;
  fPiminus_flux             = 0;
  fPiplus_flux              = 0;
  fPizero_flux              = 0;
  fPositron_flux            = 0;
  fProton_flux              = 0;
  fMuon_flux                = 0;
  fOther_flux               = 0;
  fExiting_energy           = 0.0;
  fIntegral_scintillation   = 0;
  fIntegral_scintillation_E = 0.0;
}

void G4TARCHisto::StartProcessing(G4AnalysisManager* fAnalysisManager){
  fMaxRadIndex = fRadiiBin.size();
  InitializeParams();
  //CreateTuples(fAnalysisManager);
  //FillRadialExperimentalData(fAnalysisManager);

  for (G4int ii = 0; ii < fMaxRadIndex; ii++) {
    fRadii[ii] = fRadiiBin[ii];
    fRadial_energies[ii] = fExptEnergyBin[ii];
  }

  G4double mean_energy = 0.0, fLithium_mean_energy = 0.0;
  for (G4int ii = 0; ii < fMaxTestFluxData; ii++){
    fFlux_energy[ii]     = fExptEnergyTables[2][ii];
    fFlux_data[ii]       = fExptFluxTables[2][ii];
    mean_energy         = 0.5 * (fFlux_energy[ii] + fExptEnergyTables[2][ii + 1]);
    fFlux[ii]            = fExptFluxTables[2][ii];
    fEFlux[ii]           = fExptFluxTables[2][ii] * mean_energy;
    fIntegral_Eflux     += fExptFluxTables[2][ii] * mean_energy / 1.0e6;
    fFine_energy[ii]     = fExptFluxTables[2][ii];
    fFlux_stat_error[ii] = fExptFluxErrTables[2][ii];
  }

  G4int k_idx = 0, m_idx = 0;
  G4double start_energy   = 0.01;
  G4int radial_index = 0;
  G4double bin_width      = (std::log(1.0e-5) - std::log(start_energy)) / 100; // check Appendix A NIM paper
  fLowEnergy[0]            = start_energy;
  fLithium_energy[0]       = start_energy;
  for (G4int ii = 0; ii < fMaxlowFluxData; ii++) {
    fLowEnergy[ii + 1]      = std::exp(bin_width + std::log(fLowEnergy[ii]));
    fLithium_energy[ii + 1] = std::exp(bin_width + std::log(fLithium_energy[ii]));
    mean_energy            = std::exp(0.5 * (std::log(fLowEnergy[ii + 1]) + std::log(fLowEnergy[ii])));
    fLithium_mean_energy    = std::exp(0.5 * (std::log(fLithium_energy[ii + 1]) + std::log(fLithium_energy[ii])));
    if (   (fLithium_energy[ii] < fExptEnergyBin[radial_index])
        && (fLithium_energy[ii + 1] > fExptEnergyBin[radial_index])
      ) {
        fLithium_radial_energy_lower[radial_index] = fLithium_energy[ii];
        fLithium_radial_energy_upper[radial_index] = fLithium_energy[ii + 1];
        fLithium_radial_mean[radial_index]         = fExptEnergyBin[radial_index];
        ++radial_index;
    }
    if (std::abs(mean_energy/fLowEnergy[k_idx] -1.0) < 0.05){
      fLowFluxData[ii]    = fExptFluxTables[0][k_idx];
      fLow_stat[ii]       = fExptFluxErrTables[0][k_idx];
      ++k_idx;
    }
    if (m_idx < 95 && (std::abs(fLithium_mean_energy/fExptEnergyTables[1][m_idx] - 1.0) < 0.05)){
      fLithium_flux_data[ii] = fExptFluxTables[1][m_idx];
      fLithium_stat[ii]      = fExptFluxErrTables[1][m_idx];
      fLithium_integral     += fExptFluxTables[1][m_idx];
      ++m_idx;
    }
  }
}


// trying to access HitsMap of RUN by MultiFunctional detector name and collection name
G4THitsMap<G4double>* G4TARCHisto::GetHitsMap(const G4String& detName, const G4String& colName) {
  G4String fullName = detName + "/" + colName;
  return GetHitsMap(fullName);
}

// trying to access HitsMap of RUN by <Multifunctional Det name>/<Primitive scorer>
G4THitsMap<G4double>* G4TARCHisto::GetHitsMap(const G4String& fullName) {
  for(unsigned ii = 0; ii < fCollName.size(); ii++){
    if (fCollName[ii] == fullName)
      return fRunMap[ii];
  }
  return NULL;
}

/*
void G4TARCHisto::Merge(const G4Run* aRun) {
  const G4TARCHisto* localRun = static_cast<const G4TARCHisto *> (aRun);
  for (G4int ii = 0; ii < localRun->fCollID.size(); ii++) {
    if (localRun->fCollID[ii] >= 0)
      *fRunMap[ii]            += *localRun->fRunMap[ii];
  }
  fExiting_flux                += localRun->fExiting_flux;
  fExiting_energy              += localRun->fExiting_energy;
  fExiting_check_flux          += localRun->fExiting_check_flux;
  fGamma_flux                  += localRun->fGamma_flux;
  fProton_flux                 += localRun->fProton_flux;
  fNeutron_flux                += localRun->fNeutron_flux;
  fElectron_flux               += localRun->fElectron_flux;
  fPositron_flux               += localRun->fPositron_flux;
  fPiminus_flux                += localRun->fPiminus_flux;
  fPiplus_flux                 += localRun->fPiplus_flux;
  fPizero_flux                 += localRun->fPizero_flux;
  fMuon_flux                   += localRun->fMuon_flux;
  fOther_flux                  += localRun->fOther_flux;
  fNeutron_check               += localRun->fNeutron_check;
  fNeutron_fluence             += localRun->fNeutron_fluence;
  fIntegral_flux_p             += localRun->fIntegral_flux_p;
  fIntegral_Eflux_p            += localRun->fIntegral_Eflux_p;
  fTotal_flux                  += localRun->fTotal_flux;

  for (int ii = 0; ii < fMaxFluxData; ii++){
    fFlux[ii]                   += localRun->fFlux[ii];
    fCos_flux[ii]               += localRun->fCos_flux[ii];
    fFluence_step[ii]           += localRun->fFluence_step[ii];
    fFluence_front_step[ii]     += localRun->fFluence_front_step[ii];
    fFluence_step_cyl[ii]       += localRun->fFluence_step_cyl[ii];
    fFluence_step_shell[ii]     += localRun->fFluence_step_shell[ii];
    fEFlux[ii]                  += localRun->fEFlux[ii];
    fFine_eFlux[ii]             += localRun->fFine_eFlux[ii];
  }

     for (int ii = 0; ii < fMaxlowFluxData; ii++) {
       fLowFlux[ii]                     += localRun->fLowFlux[ii];
       fCos_low_flux[ii]                += localRun->fCos_low_flux[ii];
       fLow_fluence_step[ii]            += localRun->fLow_fluence_step[ii];
       fLow_fluence_front_step[ii]      += localRun->fLow_fluence_front_step[ii];
       fLow_fluence_step_cyl[ii]        += localRun->fLow_fluence_step_cyl[ii];
       fLow_fluence_step_shell[ii]      += localRun->fLow_fluence_step_shell[ii];
       fLithium_flux[ii]                += localRun->fLithium_flux[ii];
       fLithium_fluence[ii]             += localRun->fLithium_fluence[ii];
       fLithium_fluence_Step[ii]        += localRun->fLithium_fluence_Step[ii];
       fLithium_fluence_front_Step[ii]  += localRun->fLithium_fluence_front_Step[ii];
       fLithium_fluence_cyl[ii]         += localRun->fLithium_fluence_cyl[ii];
       fLithium_Zflux[ii]               += localRun->fLithium_Zflux[ii];
     }
     for (G4int ii = 0; ii < fMaxRadIndex; ii++) {
       for (G4int jj = 0; jj < (fMaxFluxData + fMaxDimEnInt + 1); jj++){
         fRadial_fluence_step[jj][ii] += localRun->fRadial_fluence_step[jj][ii];
          }
        }
        G4Run::Merge(aRun);
      }
*/
