#include <string.h>
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "Riostream.h"

void testRootMacro() {
  TFile*   tf1 = new TFile("res.root", "r");

  TH2* hgGamma = (TH2*) tf1->Get("Gamma");
  TH2* hgNeutEnergy = (TH2*) tf1->Get("NeutronEnergy");

  TNtuple* exitingTuple = (TNtuple*) tf1->Get("h3_N_Exiting");
  TNtuple* flux4002     = (TNtuple*) tf1->Get("h4_Flux_4002");
  TNtuple* flux4004     = (TNtuple*) tf1->Get("h5_Flux_4004");
  TNtuple* flux4005     = (TNtuple*) tf1->Get("h6_Flux_4005");
  TNtuple* tarcRad      = (TNtuple*) tf1->Get("h9_Rad_Shell_Fluence");
  TNtuple* tarcRadLi    = (TNtuple*) tf1->Get("h10_Rad_Fluence_Expt_Li_Data");
  TNtuple* tarcRadHe3   = (TNtuple*) tf1->Get("h11_Rad_Fluence_Expt_He3_Data");

  TTree* h9RadShellFluence = (TTree*) tf1->Get("h9_Rad_Shell_Fluence");

  const int xHiCnt = 22;
  const int xLoCnt = 102;
  Float_t fXbinHi[xHiCnt] = {59500, 109000, 158500, 208000, 257500, 307000, 356500, 406000,
    455500, 505000, 554500, 604000, 653500, 703000, 752500, 802000, 901000, 1000000, 1162308, 1350960, 1570232, 1825092};

  Float_t fXbinLo[xLoCnt];
  Float_t fBinWidth = (std::log(1.0e+6) - std::log(0.01)) / 100.0;
  fXbinLo[0] = 0.01;
  int index = 0;
  for(int i = 1; i < xLoCnt; i++) {
    fXbinLo[i] = std::exp(fBinWidth+std::log(fXbinLo[i-1]));
  }

  TH1F* ExitingSpec            = new TH1F("ExitSpect", "TARC Exiting Neutron Spectrum", 100, 0.01, 2.0);

  TH1F* TARCDataFluenceHi      =  new TH1F("FluenceDataHi", "TARC Fluence Data (High)", (xHiCnt - 1), fXbinHi);
  TH1F* TARCDataFluenceHe3     =  new TH1F("FluenceDataHe", "TARC Fluence Data (He3)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCDataFluenceLi      =  new TH1F("FluenceDataLi", "TARC Fluence Data (Li)", (xLoCnt - 1), fXbinLo);

  TH1F* TARCDataFluenceHiErr   =  new TH1F("FluenceDataErrHi", "TARC Fluence Data (High) Error", (xHiCnt - 1), fXbinHi);
  TH1F* TARCDataFluenceHe3Err  =  new TH1F("FluenceDataErrHe", "TARC Fluence Data (He3) Error", (xLoCnt - 1), fXbinLo);
  TH1F* TARCDataFluenceLiErr   =  new TH1F("FluenceDataErrLi", "TARC Fluence Data (Li) Error", (xLoCnt - 1), fXbinLo);

  TH1F* TARCG4FluenceHi        =  new TH1F("FluenceG4Hi", "TARC G4 Fluence (High)", (xHiCnt - 1), fXbinHi);
  TH1F* TARCG4FluenceHe3       =  new TH1F("FluenceG4He", "TARC G4 Fluence (He3)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCG4FluenceLi        =  new TH1F("FluenceG4Li", "TARC G4 Fluence (Li)", (xLoCnt - 1), fXbinLo);

  TH1F* TARCG4RatioHi          =  new TH1F("RatioG4Hi", "TARC Fluence Ratio G4 Fluence (High)", (xHiCnt - 1), fXbinHi);
  TH1F* TARCG4RatioHe3         =  new TH1F("RatioG4He", "TARC Fluence Ratio G4 Fluence (He3)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCG4RatioLi          =  new TH1F("RatioG4Li", "TARC Fluence Ratio G4 Fluence (Li)", (xLoCnt - 1), fXbinLo);


  TH2F* radialHisto            = new TH2F("Radial", "TARC radial", 10000, -1000, 1000, 1000, 1.0, 50.0e+6);

  double energy, tarcflux, g4flux, g4Err, g4perp, g4fluence, g4shell, rawflux, errstat;

  for (int irow = 0; irow < exitingTuple->GetEntries(); irow++){
    exitingTuple->SetBranchAddress("energy", &energy);
    exitingTuple->GetEntry(irow);
    ExitingSpec->Fill(energy);
  }

  for (int irow = 0; irow < flux4002->GetEntries(); irow++){
    flux4002->SetBranchAddress("energy", &energy);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("tarcflux", &tarcflux);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4flux", &g4flux);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4err", &g4Err);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4perp", &g4perp);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4fluence", &g4fluence);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("rawflux", &rawflux);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("errstat", &errstat);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4_shell", &g4shell);

    TARCDataFluenceHi->Fill(energy, tarcflux);
    TAxis* xAxis = TARCDataFluenceHi->GetXaxis();
    Int_t binX   = xAxis->FindBin(energy);
    TARCDataFluenceHi->SetBinError(binX, errstat);

    TARCDataFluenceHiErr->Fill(energy, tarcflux);
    TAxis* xAxis2 = TARCDataFluenceHiErr->GetXaxis();
    Int_t binX2   = xAxis2->FindBin(energy);
    TARCDataFluenceHiErr->SetBinError(binX2, errstat);

    double corrG4perp = g4shell;
    TARCG4FluenceHi->Fill(energy, corrG4perp);
    xAxis = TARCG4FluenceHi->GetXaxis();
    binX = xAxis->FindBin(energy);
    TARCG4FluenceHi->SetBinError(binX, g4Err);
    if (tarcflux !=0.0){
      double ratio = corrG4perp / tarcflux;
      TARCG4RatioHi->Fill(energy, ratio);
    }    
  }

  for (int irow = 0; irow < flux4004->GetEntries(); irow++){
    flux4004->SetBranchAddress("energy", &energy);
    flux4004->GetEntry(irow);
    flux4004->SetBranchAddress("tarcflux", &tarcflux);
    flux4004->GetEntry(irow);
    flux4004->SetBranchAddress("g4flux", &g4flux);
    flux4004->GetEntry(irow);
    flux4004->SetBranchAddress("g4perp", &g4perp);
    flux4004->GetEntry(irow);
    flux4004->SetBranchAddress("g4fluence", &g4fluence);
    flux4004->GetEntry(irow);
    flux4004->SetBranchAddress("rawflux", &rawflux);
    flux4004->GetEntry(irow);
    flux4004->SetBranchAddress("g4_shell", &g4shell);
    flux4004->GetEntry(irow);
    flux4004->SetBranchAddress("g4err", &g4Err);
    flux4004->GetEntry(irow);
    flux4004->SetBranchAddress("errstat", &errstat);
    flux4004->GetEntry(irow);

    TARCDataFluenceHe3->Fill(energy, tarcflux);
    TAxis* xAxis = TARCDataFluenceHe3->GetXaxis();
    Int_t binX   = xAxis->FindBin(energy);
    TARCDataFluenceHe3->SetBinError(binX, errstat);

    TARCDataFluenceHe3Err->Fill(energy, tarcflux);
    TAxis* xAxis2 = TARCDataFluenceHe3Err->GetXaxis();
    Int_t binX2   = xAxis2->FindBin(energy);
    TARCDataFluenceHe3Err->SetBinError(binX2, errstat);

    double corrG4perp = g4shell;
    TARCG4FluenceHe3->Fill(energy, corrG4perp);
    xAxis = TARCG4FluenceHe3->GetXaxis();
    binX = xAxis->FindBin(energy);
    TARCG4FluenceHe3->SetBinError(binX, g4Err);
    if (tarcflux != 0.0){
      double ratio = corrG4perp / tarcflux;
      TARCG4RatioHe3->Fill(energy, ratio);
    }
  }


  for (int irow = 0; irow < flux4005->GetEntries(); irow++){
    flux4005->SetBranchAddress("energy", &energy);
    flux4005->GetEntry(irow);
    flux4005->SetBranchAddress("tarcflux", &tarcflux);
    flux4005->GetEntry(irow);
    flux4005->SetBranchAddress("g4flux", &g4flux);
    flux4005->GetEntry(irow);
    flux4005->SetBranchAddress("g4perp", &g4perp);
    flux4005->GetEntry(irow);
    flux4005->SetBranchAddress("g4fluence", &g4fluence);
    flux4005->GetEntry(irow);
    flux4005->SetBranchAddress("rawflux", &rawflux);
    flux4005->GetEntry(irow);
    flux4005->SetBranchAddress("errstat", &errstat);
    flux4005->GetEntry(irow);
    flux4005->SetBranchAddress("g4err", &g4Err);
    flux4005->GetEntry(irow);
    flux4005->SetBranchAddress("g4_shell", &g4shell);
    flux4005->GetEntry(irow);

    TARCDataFluenceLi->Fill(energy, tarcflux);
    TAxis* xAxis = TARCDataFluenceLi->GetXaxis();
    Int_t binX   = xAxis->FindBin(energy);
    TARCDataFluenceLi->SetBinError(binX, errstat);

    TARCDataFluenceLiErr->Fill(energy, tarcflux);
    TAxis* xAxis2 = TARCDataFluenceLiErr->GetXaxis();
    Int_t binX2   = xAxis2->FindBin(energy);
    TARCDataFluenceLiErr->SetBinError(binX2, errstat);

    double corrG4perp = g4shell;
    TARCG4FluenceLi->Fill(energy, corrG4perp);
    xAxis = TARCG4FluenceLi->GetXaxis();
    binX = xAxis->FindBin(energy);
    TARCG4FluenceLi->SetBinError(binX, g4Err);
    if (tarcflux != 0.0){
      double ratio = corrG4perp / tarcflux;
      TARCG4RatioLi->Fill(energy, ratio);
    }
  }

  TARCDataFluenceHi->SetMarkerStyle(5);
  TARCDataFluenceHe3->SetMarkerStyle(6);
  TARCDataFluenceLi->SetMarkerStyle(7);
  TARCDataFluenceHiErr->SetMarkerStyle(8);
  TARCDataFluenceHe3Err->SetMarkerStyle(9);
  TARCDataFluenceLiErr->SetMarkerStyle(10);
  
  TARCDataFluenceHi->SetLineColor(kRed);
  TARCDataFluenceHi->SetLineWidth(0.85);
  TARCDataFluenceHe3->SetLineColor(kBlue);
  TARCDataFluenceHe3->SetLineWidth(0.85);
  TARCDataFluenceLi->SetLineColor(kCyan);
  TARCDataFluenceLi->SetLineWidth(0.85);
  TARCDataFluenceHiErr->SetLineColor(kGreen);
  TARCDataFluenceHiErr->SetLineWidth(0.7);
  TARCDataFluenceHe3Err->SetLineColor(kYellow);
  TARCDataFluenceHe3Err->SetLineWidth(0.7);
  TARCDataFluenceLiErr->SetLineColor(kMagenta);
  TARCDataFluenceLiErr->SetLineWidth(0.7);
  TARCDataFluenceHi->Draw("SAME E1");
  TARCDataFluenceHe3->Draw("SAME E1");
  TARCDataFluenceLi->Draw("SAME E1");
  TARCDataFluenceHiErr->Draw("SAME E");
  TARCDataFluenceHe3Err->Draw("SAME E");
  TARCDataFluenceLiErr->Draw("SAME E");

  TCanvas*  c1 = new TCanvas("c1", "TARC Study", 700, 500);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.2);
  gPad->SetLogx();
  gPad->SetLogy();

  gPad->DrawFrame(0.001, 1000, 1000000, 25000000," Energy/eV; EdF/dE n/cm^{2}/10^{9}p")->GetXaxis()->SetTitleOffset(1.2);
  TString thisTitlePart1 = "TARC Output";
  TString thisTitlePart10 = "TARC_Output";
  TString thisTitlePart2 = "1.5 GeV/c protons";
  TString thisTitlePart3 = "100 events.";
  TString plotTitle1 = thisTitlePart1 + thisTitlePart2 + thisTitlePart3;
  TLatex* tlx = new TLatex(0.23, 0.93, plotTitle1);
  tlx->SetNDC(kTRUE);
  tlx->Draw();
  
  TString savedFile1 = thisTitlePart10 + "_fluence1.png";
  c1->Print(savedFile1);
  c1->Close();

  TCanvas*  c2 = new TCanvas("c2", "TARC Ratio G4/Data Study", 700, 500);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.2);
  gPad->SetLogx();
  //  gPad->SetLogy();
  gPad->DrawFrame(0.001, 1000, 1000000, 25000000," Energy/eV; G4/Data")->GetXaxis()->SetTitleOffset(1.2);

  TARCG4RatioHi->SetMarkerStyle(15);
  TARCG4RatioHe3->SetMarkerStyle(16);
  TARCG4RatioLi->SetMarkerStyle(17);
  
  TARCG4RatioHi->SetMarkerColor(kRed);
  TARCG4RatioHi->SetMarkerStyle(21);
  TARCG4RatioHi->SetLineColor(kRed);
  TARCG4RatioHi->SetLineWidth(0.85);
  TARCG4RatioHe3->SetMarkerColor(kBlue);
  TARCG4RatioHe3->SetMarkerStyle(22);
  TARCG4RatioHe3->SetLineColor(kBlue);
  TARCG4RatioHe3->SetLineWidth(0.85);
  TARCG4RatioLi->SetMarkerColor(kCyan);
  TARCG4RatioLi->SetMarkerStyle(23);
  TARCG4RatioLi->SetLineColor(kCyan);
  TARCG4RatioLi->SetLineWidth(0.85);
  TARCG4RatioHi->Draw("SAME P");
  TARCG4RatioHe3->Draw("SAME P");
  TARCG4RatioLi->Draw("SAME P");

  TString savedFile2 = thisTitlePart10 + "_Ratio.png";
  c2->Print(savedFile2);
  c2->Close();
 

  tf1->Close();
}
