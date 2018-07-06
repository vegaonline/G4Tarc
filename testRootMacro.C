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
  TFile*   tf1 = TFile::Open("res.root");

  TH2* hgGamma                  = (TH2*) tf1->Get("Gamma");
  TH2* hgNeutEnergy           = (TH2*) tf1->Get("NeutronEnergy");
  TH2* hgElectronEnergy     = (TH2*) tf1->Get("ElectronEdep");
  TH2* hgPositEnergy          = (TH2*) tf1->Get("PositronEdep");
  TH2* hgNeutEnergyTime  = (TH2*) tf1->Get("NeutronEnergyTime");
  TH2* hgOtherEnergyTime = (TH2*) tf1->Get("OtherParticleEnergyTime");
  TH2* histo1                        = (TH2*)tf1->Get("Gamma");
  TH2* histo2                        = (TH2*)tf1->Get("NeutronEnergy");

  TNtuple* exitingTuple        = (TNtuple*) tf1->Get("h3_N_Exiting");
  TNtuple* flux4002             = (TNtuple*) tf1->Get("h4_Flux_4002");
  TNtuple* flux4004             = (TNtuple*) tf1->Get("h5_Flux_4004");
  TNtuple* flux4005             = (TNtuple*) tf1->Get("h6_Flux_4005");
  TNtuple* tarcRad               = (TNtuple*) tf1->Get("h9_Rad_Shell_Fluence");
  TNtuple* tarcRadLi           = (TNtuple*) tf1->Get("h10_Rad_Fluence_Expt_Li_Data");
  TNtuple* tarcRadHe3        = (TNtuple*) tf1->Get("h11_Rad_Fluence_Expt_He3_Data");

  TTree* h9RadShellFluence = (TTree*) tf1->Get("h9_Rad_Shell_Fluence");

  const int xHiCnt =  22;
  const int xLoCnt = 102;
  Float_t fXbinHi[xHiCnt] = {59500, 109000, 158500, 208000, 257500, 307000, 356500, 406000,
    455500, 505000, 554500, 604000, 653500, 703000, 752500, 802000, 901000, 1000000, 1162308, 1350960, 1570232, 1825092};

  Float_t fXbinLo[xLoCnt];
  Float_t fBinWidth = (std::log(1.0e+5) - std::log(0.01)) / 100.0;
  fXbinLo[0] = 0.01;
  int index = 0;
  for(int i = 1; i < xLoCnt; i++) {
    fXbinLo[i] = std::exp(fBinWidth + std::log(fXbinLo[i-1]));
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


  TH2F* radialHisto            = new TH2F("Radial", "TARC radial", 1000, -1000, 1000, 1000, 1.0, 2.5e+7);

  double energy, tarcflux, g4flux, g4Err, g4perp, g4fluence, g4shell, rawflux, g4shellerr, eflux, errstat;

  for (int irow = 0; irow < exitingTuple->GetEntries(); ++irow){
    exitingTuple->SetBranchAddress("energy", &energy);
    exitingTuple->GetEntry(irow);
    ExitingSpec->Fill(energy);
  }

  for (int irow = 0; irow < flux4002->GetEntries(); ++irow){
    flux4002->SetBranchAddress("energy", &energy);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("tarcflux", &tarcflux);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("errstat", &errstat);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4flux", &g4flux);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4perp", &g4perp);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4fluence", &g4fluence);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4err", &g4Err);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("rawflux", &rawflux);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("eflux", &eflux);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4_shell", &g4shell);
    flux4002->GetEntry(irow);
    flux4002->SetBranchAddress("g4_shell_err", &g4shellerr);
    flux4002->GetEntry(irow);

    TARCDataFluenceHi->Fill(energy, tarcflux);
    TAxis* xAxis1 = TARCDataFluenceHi->GetXaxis();
    Int_t binX1   = xAxis1->FindBin(energy);
    TARCDataFluenceHi->SetBinError(binX1, errstat);

    TARCDataFluenceHiErr->Fill(energy, tarcflux);
    TAxis* xAxis2 = TARCDataFluenceHiErr->GetXaxis();
    Int_t binX2   = xAxis2->FindBin(energy);
    TARCDataFluenceHiErr->SetBinError(binX2, errstat);

    double corrG4perp = g4shell ;  // g4perp;
    TARCG4FluenceHi->Fill(energy, corrG4perp);
    TAxis* xAxis3 = TARCG4FluenceHi->GetXaxis();
    Int_t binX3 = xAxis3->FindBin(energy);
    TARCG4FluenceHi->SetBinError(binX3, g4Err);
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
    //flux4004->SetBranchAddress("g4err", &g4Err);
    //flux4004->GetEntry(irow);
    flux4004->SetBranchAddress("errstat", &errstat);
    flux4004->GetEntry(irow);

    TARCDataFluenceHe3->Fill(energy, tarcflux);
    TAxis* xAxis1 = TARCDataFluenceHe3->GetXaxis();
    Int_t binX1   = xAxis1->FindBin(energy);
    TARCDataFluenceHe3->SetBinError(binX1, errstat);

    TARCDataFluenceHe3Err->Fill(energy, tarcflux);
    TAxis* xAxis2 = TARCDataFluenceHe3Err->GetXaxis();
    Int_t binX2   = xAxis2->FindBin(energy);
    TARCDataFluenceHe3Err->SetBinError(binX2, errstat);

    double corrG4perp =  g4shell ; // g4perp;
    TARCG4FluenceHe3->Fill(energy, corrG4perp);
    TAxis* xAxis3 = TARCG4FluenceHe3->GetXaxis();
    Int_t binX3 = xAxis3->FindBin(energy);
    TARCG4FluenceHe3->SetBinError(binX3, g4Err);
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
    //flux4005->SetBranchAddress("g4err", &g4Err);
    //flux4005->GetEntry(irow);
    flux4005->SetBranchAddress("g4_shell", &g4shell);
    flux4005->GetEntry(irow);

    TARCDataFluenceLi->Fill(energy, tarcflux);
    TAxis* xAxis1 = TARCDataFluenceLi->GetXaxis();
    Int_t binX1   = xAxis1->FindBin(energy);
    TARCDataFluenceLi->SetBinError(binX1, errstat);

    TARCDataFluenceLiErr->Fill(energy, tarcflux);
    TAxis* xAxis2 = TARCDataFluenceLiErr->GetXaxis();
    Int_t binX2   = xAxis2->FindBin(energy);
    TARCDataFluenceLiErr->SetBinError(binX2, errstat);

    double corrG4perp =  g4shell ;  // g4perp;
    TARCG4FluenceLi->Fill(energy, corrG4perp);
    TAxis* xAxis3 = TARCG4FluenceLi->GetXaxis();
    Int_t binX3 = xAxis3->FindBin(energy);
    TARCG4FluenceLi->SetBinError(binX3, g4Err);
    if (tarcflux != 0.0){
      double ratio = corrG4perp / tarcflux;
      TARCG4RatioLi->Fill(energy, ratio);
    }
  }

  TLatex* tlx;

  TCanvas* c0 = new TCanvas("c0", "TARC Summary Report", 1020, 800);
  c0->Divide(2, 2);
  c0->cd(1);  // left top gPad
  gPad->SetLogx();
  gPad->SetLogy();
  TARCG4FluenceHi->SetTitle("G4 Flux");
  TARCG4FluenceHi->GetXaxis()->SetTitle("Energy / eV");
  TARCG4FluenceHi->GetYaxis()->SetTitle("Flux ( dN/dE / source gamma)");
  TARCG4FluenceHi->GetYaxis()->SetTitleOffset(1.2);
  TARCG4FluenceHi->SetMarkerStyle(3);
  TARCG4FluenceHi->SetMarkerColor(kBlue + 3);
  TARCG4FluenceHi->Draw("SAME");

  c0->cd(2); // Right Top gPad
  gPad->SetLogx();
  gPad->SetLogy();
  TARCG4FluenceLi->SetTitle("G4 Fluence Li Data");
  TARCG4FluenceLi->GetXaxis()->SetTitle("Energy / eV");
  TARCG4FluenceLi->GetYaxis()->SetTitle("Flux ( dN/dE / source gamma)");
  TARCG4FluenceLi->GetYaxis()->SetTitleOffset(1.4);
  TARCG4FluenceLi->SetLineColor(kBlue - 2);
  TARCG4FluenceLi->Draw("SAME");

  c0->cd(3); // Left Bottom gPad
  gPad->SetLogx();
  gPad->SetLogy();
  gStyle->SetHistLineWidth(3);

  tlx=new TLatex(0.12, 0.93, "TARC fluence comparison  for 1.5 GeV/c protons ");
  tlx->SetTextSize(0.5);
  tlx->SetNDC(kTRUE);
  tlx->Draw();
  gPad->DrawFrame(1.0e-3, 1.0e3, 5.0e6, 1.0e8,"; Energy/eV; EdF/dE n/cm^{2}/10^{9}p")->GetXaxis()->SetTitleOffset(1.2);

  TARCDataFluenceHi->SetMarkerStyle(kFullCircle);
  TARCDataFluenceHi->SetMarkerSize(0.8);
  TARCDataFluenceHi->SetMarkerColor(kBlue );
  TARCDataFluenceHi->Draw("SAME E1");

  TARCDataFluenceHe3->SetMarkerStyle(kOpenSquare);
  TARCDataFluenceHe3->SetMarkerSize(0.7);
  TARCDataFluenceHe3->SetMarkerColor(kRed );
  TARCDataFluenceHe3->Draw("SAME E1");

  TARCDataFluenceLi->SetMarkerStyle(kDiamond);
  TARCDataFluenceLi->SetMarkerColor(kYellow -1);
  TARCDataFluenceLi->Draw("SAME E1");

  //TARCDataFluenceHiErr->Draw("SAME E1");
  //TARCDataFluenceHe3Err->Draw("SAME E1");
  //TARCDataFluenceLiErr->Draw("SAME E1");

  double start = 0.13, stop = 0.9;
  double xwidth  = 0.35, ywidth = 0.14;

  auto legend = new TLegend(start, stop - ywidth, start + xwidth, stop);
  legend->SetTextFont(62);
  legend->SetHeader("Fluence comparison","C"); // option "C" allows to center the header
  legend->SetTextFont(42);
  legend->AddEntry(TARCDataFluenceHi,"Histogram for G4 fluence","ep");
  legend->AddEntry(TARCDataFluenceHe3,"Histogram for He3 data fluence","ep");
  legend->AddEntry(TARCDataFluenceLi,"Histogram for Li data fluence","ep");   // ep for errors and points
  legend->Draw();

  c0->cd(4); // Bottom Right gPad
  gPad->SetLogx();
  //gPad->SetLogy();
  //hgGamma->Draw("COLZ");
  hgNeutEnergy->SetTitle("Neutron Deposition");
  hgNeutEnergy->SetLineColor(kRed);
  hgNeutEnergy->SetMarkerColor(kRed);
  hgNeutEnergy->Draw("COLZ");
  hgNeutEnergy->GetXaxis()->SetTitle("log10(Energy (eV))");
  hgNeutEnergy->GetXaxis()->SetTitleSize(0.03);
  hgNeutEnergy->GetXaxis()->SetTitleOffset(1.2);
  hgNeutEnergy->GetYaxis()->SetTitle("Neutrons");
  hgNeutEnergy->GetYaxis()->SetTitleOffset(1.3);
  //hgNeutEnergy->Draw("SAME COLZ");
  //hgGamma->Draw("SAME COLZ");

  c0->Print("TARC_Report_Summary.png");
  c0->Close();


  TCanvas*  c1 = new TCanvas("c1", "TARC Study", 900, 700);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.2);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->DrawFrame(1.0e-3, 1.0e3, 5.0e6, 1.0e8,"; log10(Energy/eV); EdF/dE n/cm^{2}/10^{9}p")->GetXaxis()->SetTitleOffset(1.2);
  TString thisTitlePart1 = "TARC Output";
  TString thisTitlePart10 = "TARC_Output";
  TString thisTitlePart2 = " 1.5 GeV/c protons";
  TString thisTitlePart3 = " for 500 events.";
  TString plotTitle1 = thisTitlePart1 + thisTitlePart2 + thisTitlePart3;

  tlx = new TLatex(0.15, 0.93, "TARC Fluence comparison for 1.5 GeV/c proton");    // plotTitle1);
  tlx->SetTextSize(0.04);
  tlx->SetNDC(kTRUE);
  tlx->Draw();

  TARCDataFluenceHi->SetMarkerStyle(kFullCircle);
  TARCDataFluenceHi->Draw("SAME E1");
  TARCDataFluenceHe3->SetMarkerStyle(kFullSquare);
  TARCDataFluenceHe3->Draw("SAME E1");
  TARCDataFluenceLi->SetMarkerStyle(kFullTriangleUp);
  TARCDataFluenceLi->Draw("SAME E1");

  TARCDataFluenceHiErr->Draw("SAME E");
  TARCDataFluenceHe3Err->Draw("SAME E");
  TARCDataFluenceLiErr->Draw("SAME E");

  start = 0.13, stop = 0.9;
  xwidth  = 0.35, ywidth = 0.14;

  legend = new TLegend(start, stop - ywidth, start + xwidth, stop);
  legend->SetTextFont(62);
  legend->SetHeader("Fluence comparison","C"); // option "C" allows to center the header
  legend->SetTextFont(42);
  legend->AddEntry(TARCDataFluenceHi,"Histogram for G4 fluence","ep");
  legend->AddEntry(TARCDataFluenceHe3,"Histogram for He3 data fluence","ep");
  legend->AddEntry(TARCDataFluenceLi,"Histogram for Li data fluence","ep");   // ep for errors and points
  legend->Draw();

  TString savedFile1 = thisTitlePart10 + "_fluence.png";
  c1->Print(savedFile1);
  c1->Close();

  TCanvas*  c2 = new TCanvas("c2", "TARC Ratio G4/Data Study", 900, 700);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.2);
  gPad->SetLogx();
  gPad->SetLogy();
  tlx = new TLatex(0.15, 0.93, "TARC G4/Data Ratio comparison for 1.5 GeV/c proton");    // plotTitle1);
  tlx->SetTextSize(0.04);
  tlx->SetNDC(kTRUE);
  tlx->Draw();

  gPad->DrawFrame(5.0e-3, 5.0e-5, 5e6, 5.0e-1,"; log10(Energy/eV); G4 Shell / Data")->GetXaxis()->SetTitleOffset(1.2);
  //gPad->DrawFrame(4e4, 1.0e-3, 5e6, 5.0e-1,"; Energy/eV; G4/Data")->GetXaxis()->SetTitleOffset(1.2);  // for the Hi Data
  TARCG4RatioHi->SetMarkerColor(kRed);
  TARCG4RatioHi->SetMarkerStyle(kFullCircle);
  TARCG4RatioHi->SetLineColor(kRed);
  TARCG4RatioHi->SetLineWidth(0.85);
  TARCG4RatioHi->SetStats(kTRUE);
  TARCG4RatioHi->Draw("SAME P");

  TARCG4RatioHe3->SetMarkerColor(kBlue);
  TARCG4RatioHe3->SetMarkerStyle(kFullSquare);
  TARCG4RatioHe3->SetLineColor(kBlue);
  TARCG4RatioHe3->SetLineWidth(0.85);
  TARCG4RatioHe3->SetStats(kTRUE);
  TARCG4RatioHe3->Draw("SAME PLC PMC");

  TARCG4RatioLi->SetMarkerColor(kGreen);
  TARCG4RatioLi->SetMarkerStyle(kFullTriangleUp);
  TARCG4RatioLi->SetLineColor(kGreen);
  TARCG4RatioLi->SetLineWidth(0.85);
  TARCG4RatioLi->SetStats(kTRUE);
  TARCG4RatioLi->Draw("SAME PLC PMC");

  start = 0.13, stop = 0.9;
  xwidth  = 0.35, ywidth = 0.14;
  legend = new TLegend(start, stop - ywidth, start + xwidth, stop);
  legend->SetTextFont(62);
  legend->SetHeader("G4/Data Ratio Comparison","C"); // option "C" allows to center the header
  legend->SetTextFont(42);
  legend->AddEntry(TARCG4RatioHi,"Shell / Data Fluence Ratio","ep");
  legend->AddEntry(TARCG4RatioHe3,"Shell / He3 data Fluence Ratio","ep");
  legend->AddEntry(TARCG4RatioLi,"Shell / Li data Fluence Ratio","ep");   // ep for errors and points
  legend->Draw();
  TString savedFile2 = thisTitlePart10 + "_Ratio.png";
  gStyle->SetTitle(thisTitlePart10+" Ratio");
  c2->Print(savedFile2);
  c2->Close();

  TCanvas*  c3 = new TCanvas("c3", "TARC Gamma Energy Deposition Study", 900, 500);
  c3->Divide(2, 1);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.2);
  hgGamma->SetMarkerStyle(7);
  hgGamma->SetMarkerColor(kBlue);
  hgGamma->SetLineColor(kBlue);
  hgGamma->GetXaxis()->SetTitle("Energy (eV)");
  hgGamma->GetYaxis()->SetTitle("Gamma");
  hgGamma->GetYaxis()->SetTitleOffset(1.4);
  hgGamma->SetTitle("Gamma Deposition");

  c3->cd(1);
  hgGamma->Draw();
  c3->cd(2);
  gPad->SetLogx();
  TH1* h1 = hgGamma->DrawCopy();
  h1->GetXaxis()->SetTitle("log10(Energy (eV))");
  c1->cd(0);
  c3->Print("TARC_Gamma_Edep.png");
  c3->Close();

  TCanvas*  c4 = new TCanvas("c4", "TARC Neutron Energy Deposition Study", 900, 500);
  c4->Divide(2, 1);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.2);
  hgNeutEnergy->SetMarkerStyle(7);
  hgNeutEnergy->SetMarkerColor(kBlue);
  hgNeutEnergy->SetLineColor(kBlue);
  hgNeutEnergy->GetXaxis()->SetTitle("Energy (eV)");
  hgNeutEnergy->GetYaxis()->SetTitle("Neutron ");
  hgNeutEnergy->GetYaxis()->SetTitleOffset(1.4);
  hgNeutEnergy->SetTitle("Neutron energy Deposition");
  c4->cd(1);
  hgNeutEnergy->Draw();
  c4->cd(2);
  gPad->SetLogx();
  h1 = hgNeutEnergy->DrawCopy();
  h1->GetXaxis()->SetTitle("log10(Energy (eV))");
  c1->cd(0);
  c4->Print("TARC_Neutron_Edep.png");
  c4->Close();

  TCanvas*  c5 = new TCanvas("c5", "TARC Electron Energy Deposition Study", 1200, 500);
  c5->Divide(2, 1);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.5);
  hgElectronEnergy->GetXaxis()->SetTitle("Energy (eV)");
  hgElectronEnergy->GetYaxis()->SetTitle("Electron");
  hgElectronEnergy->SetMarkerStyle(7);
  hgElectronEnergy->SetMarkerColor(kBlue);
  hgElectronEnergy->SetLineColor(kBlue);
  hgElectronEnergy->GetYaxis()->SetTitleOffset(1.3);
  hgElectronEnergy->SetTitle("Electron Deposition");
  c5->cd(1);
  hgElectronEnergy->Draw();
  c5->cd(2);
  gPad->SetLogx();
  TH1* h = hgElectronEnergy->DrawCopy();
  h->GetXaxis()->SetTitle("log10(Energy / eV)");
  h->GetYaxis()->SetTitleOffset(1.3);
  c1->cd(0);
  c5->Print("TARC_Electron_Edep.png");
  c5->Close();


  TCanvas*  c6 = new TCanvas("c6", "TARC Positron Energy Deposition Study", 700, 500);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.2);
  hgPositEnergy->SetTitle("Positron Deposition");
  hgPositEnergy->SetMarkerStyle(7);
  hgPositEnergy->SetMarkerColor(kBlue);
  hgPositEnergy->SetLineColor(kBlue);
  hgPositEnergy->GetXaxis()->SetTitle("Energy (eV)");
  hgPositEnergy->Draw();
  c6->Print("TARC_Positron_Edep.png");
  c6->Close();

  TCanvas*  c7 = new TCanvas("c7", "TARC Neutron Energy-Time Study", 700, 500);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.2);
  hgNeutEnergyTime->SetMarkerStyle(7);
  hgNeutEnergyTime->SetMarkerColor(kBlue);
  hgNeutEnergyTime->SetLineColor(kBlue);
  hgNeutEnergyTime->SetTitle("Neutron Energy - Time characteristics");
  hgNeutEnergyTime->GetXaxis()->SetTitle("log10(Time (#mus))");
  hgNeutEnergyTime->GetYaxis()->SetTitle("Neutron Energy (eV)");
  hgNeutEnergyTime->Draw();
  c7->Print("TARC_Neutron_Energy_Time.png");
  c7->Close();

  TCanvas*  c8 = new TCanvas("c8", "TARC Other Particles Energy-Time Study", 700, 500);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.2);
  hgOtherEnergyTime->SetMarkerStyle(7);
  hgOtherEnergyTime->SetMarkerColor(kBlue);
  hgOtherEnergyTime->SetLineColor(kBlue);
  hgOtherEnergyTime->SetTitle("Time - Energy characteristics for other particles.");
  hgOtherEnergyTime->GetXaxis()->SetTitle("Time (#mus)");
  hgOtherEnergyTime->GetYaxis()->SetTitle("Energy (eV)");
  hgOtherEnergyTime->Draw();
  c8->Print("TARC_Other_Energy_Time.png");
  c8->Close();


  TCanvas* c9 = new TCanvas("c9","TARC Exiting Neutron Spectrum", 1200, 700);
  c9->Divide(2,1);
  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(0.3);
  gStyle->SetTitleX(0.2);
  ExitingSpec->SetLineColor(kBlue);
  ExitingSpec->GetXaxis()->SetTitle("Energy (eV)");
  ExitingSpec->GetYaxis()->SetTitle("Exiting Neutron from the System");
  c9->cd(1);
  ExitingSpec->Draw();
  c9->cd(2);
  gPad->SetLogx();
  gPad->SetLogy();
  h1 = ExitingSpec->DrawCopy();
  h1->GetXaxis()->SetTitle("log10(Energy (eV))");
  //h1->GetYaxis()->SetTitle("log10(Exiting Neutron from the System)");
  c1->cd(0);
  c9->Print("TARC_Exiting_Neutron_Spectrum.png");
  c9->Close();

TCanvas* c10 = new TCanvas("c10", "TARC Radial Fluence", 900, 700);
gStyle->SetHistLineWidth(3);
gStyle->SetLineWidth(0.3);
gStyle->SetTitleX(0.2);
gPad->SetLogy();
gStyle->SetTitle("fluence");
gPad->DrawFrame(-220, 1.0e-1, 220.0, 1.0e9, "; Radial Distance / cm; dF/dE (n/cm^{2}/eV/10^{9} p)")->GetXaxis()->SetTitleOffset(1.2);
tarcRad->SetMarkerStyle(21);
tarcRad->SetMarkerColor(kRed);
tarcRad->SetMarkerSize(0.8);
tarcRad->Scan();
tarcRad->Draw("fluence/energy  : radius / 10.0", "", "SAME"); // to convert to cm^{2}
tarcRadLi->SetMarkerStyle(28);
tarcRadLi->SetMarkerSize(0.8);
tarcRadLi->SetMarkerColor(kBlue);
tarcRadLi->Draw("data : radius / 10.0", "", "SAME");   // changing to cm
tarcRadHe3->SetMarkerStyle(28);
tarcRadHe3->SetMarkerSize(0.8);
tarcRadHe3->SetMarkerColor(kGreen);
tarcRadHe3->Draw("data : radius / 10.0", "", "SAME");   // changing to cm
start = 0.11, stop = 0.9;
xwidth  = 0.28, ywidth = 0.12;
legend = new TLegend(start, stop - ywidth, start + xwidth, stop);
legend->SetTextFont(62);
legend->SetHeader("Radial Fluence distribution","C"); // option "C" allows to center the header
legend->SetTextFont(42);
legend->AddEntry(tarcRad,"Distribution for TARC simulation","ep");
legend->AddEntry(tarcRadHe3,"Distribution for He3 Data","ep");
legend->AddEntry(tarcRadLi,"Distribution for Li data","ep");   // ep for errors and points
legend->Draw();
c10->Print("TARC_Output_radial_X.png");
c10->Close();

tf1->Close();

}
