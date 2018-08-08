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

void RootDataPlotting(){
  // ROOT Data file loaded
  TFile* tf1 = TFile::Open("res.root");

  // Declaring some parameters could be used for scaling/binning
  const int xHiCnt = 22;
  const int xLoCnt = 103;
  Float_t fXbinHi[xHiCnt] = {59500, 109000, 158500, 208000, 257500, 307000, 356500, 406000, 455500, 505000, 554500,
                             604000, 653500, 703000, 752500, 802000, 901000, 1000000, 1162308, 1350960, 1570232, 1825092};
  Float_t fXbinLo[xLoCnt];
  fXbinLo[0] = 0.01;
  Float_t fBinWidth = (std::log(1.0e+5) - std::log(0.01)) / 100.0;
  for (int i = 1; i < xLoCnt; i++) fXbinLo[i] = std::exp(fBinWidth + std::log(fXbinLo[i - 1]));

  // NTuples read from the histograms in the ROOT file
  TH2* GammaED             = (TH2*) tf1->Get("Gamma");
  TH2* NeutEnergy          = (TH2*) tf1->Get("NeutronEnergy");
  TH2* ElecED              = (TH2*) tf1->Get("ElectronEdep");
  TH2* PositED             = (TH2*) tf1->Get("PositronEdep");
  TH2* OtherED             = (TH2*) tf1->Get("OtherEdep");
  TH2* PStack              = (TH2*) tf1->Get("ParticleStack");
  TH2* NeutPerEvent        = (TH2*) tf1->Get("NeutronPerEvent");
  TH2* ProtPerEvent        = (TH2*) tf1->Get("ProtonPerEvent");
  TH2* NeutronET           = (TH2*) tf1->Get("NeutronET");
  TH2* OtherET             = (TH2*) tf1->Get("OtherPartET");

  // NTuples read from Tree  in the ROOT file and kept as tuples for further computation
  TNtuple* h1SecTuple      = (TNtuple*) tf1->Get("h1_Secondary");
  TNtuple* h2NET           = (TNtuple*) tf1->Get("h2_N_ET");
  TNtuple* h3NExiting      = (TNtuple*) tf1->Get("h3_N_Exiting");
  TNtuple* Flux4002        = (TNtuple*) tf1->Get("h4_Flux_4002");
  TNtuple* Flux4004        = (TNtuple*) tf1->Get("h5_Flux_4004");
  TNtuple* Flux4005        = (TNtuple*) tf1->Get("h6_Flux_4005"); 
  TNtuple* h8SpallN        = (TNtuple*) tf1->Get("h8_Created_N");
  TNtuple* RadShellFluence = (TNtuple*) tf1->Get("h9_Rad_Shell_Fluence");  
  TNtuple* ExptLiData      = (TNtuple*) tf1->Get("h10_Rad_Fluence_Expt_Li_Data");
  TNtuple* ExptHeData      = (TNtuple*) tf1->Get("h11_Rad_Fluence_Expt_He3_Data");
  TNtuple* h15OtherET      = (TNtuple*) tf1->Get("h15_Other_ET");


  // Declare Histograms to produce  
  TH1F* ExitingSpec            = new TH1F("ExitSpect", "TARC Exiting Neutron Spectrum", 300, 0.01, 1.5);
  TH1F* TARCDataFluenceHi      = new TH1F("FluenceDataHi", "TARC Fluence Data (High)", (xHiCnt - 1), fXbinHi);
  TH1F* TARCDataFluenceHe3     = new TH1F("FluenceDataHe", "TARC Fluence Data (He3)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCDataFluenceLi      = new TH1F("FluenceDataLi", "TARC Fluence Data (Li)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCDataFluenceHiErr   = new TH1F("FluenceDataErrHi", "TARC Fluence Data (High) Error", (xHiCnt - 1), fXbinHi);
  TH1F* TARCDataFluenceHe3Err  = new TH1F("FluenceDataErrHe", "TARC Fluence Data (He3) Error", (xLoCnt - 1), fXbinLo);
  TH1F* TARCDataFluenceLiErr   = new TH1F("FluenceDataErrLi", "TARC Fluence Data (Li) Error", (xLoCnt - 1), fXbinLo);
  TH1F* TARCG4FluenceHi        = new TH1F("FluenceG4Hi", "TARC G4 Fluence (High)", (xHiCnt - 1), fXbinHi);
  TH1F* TARCG4FluenceHe3       = new TH1F("FluenceG4He", "TARC G4 Fluence (He3)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCG4FluenceLi        = new TH1F("FluenceG4Li", "TARC G4 Fluence (Li)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCG4RatioHi          = new TH1F("RatioG4Hi", "TARC Fluence Ratio G4 Fluence (High)", (xHiCnt - 1), fXbinHi);
  TH1F* TARCG4RatRatHi         = new TH1F("Ratio/(1-Ratio)", "TARC Fluence Ratio G4 Fluence (High)", (xHiCnt - 1), fXbinHi);
  TH1F* TARCG4RatioHe3         = new TH1F("RatioG4He", "TARC Fluence Ratio G4 Fluence (He3)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCG4RatRatHe3        = new TH1F("Ratio/(1-Ratio)", "TARC Fluence Ratio G4 Fluence (He3)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCG4RatioLi          = new TH1F("RatioG4Li", "TARC Fluence Ratio G4 Fluence (Li)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCG4RatRatLi         = new TH1F("Ratio/(1-Ratio)", "TARC Fluence Ratio G4 Fluence (Li)", (xLoCnt - 1), fXbinLo);
  TH2F* radialHisto            = new TH2F("Radial", "TARC radial", 1000, -1000, 1000, 1000, 1.0, 2.5e+7);

  // Variables for histograms
  double energy, tarcflux, errstat, g4flux, g4perp, g4fluence, g4error, rawflux, eflux, tarcmeanflux, abseflux;
  double radius, g4shell, g4shellerr, tarcmeanstep, tarccyl, tarcfrontflux;
  TLatex* tlx;

  // Fill the histograms
  for (int irow = 0; irow < h3NExiting->GetEntries(); ++irow){
    h3NExiting->SetBranchAddress("energy", &energy);
    h3NExiting->GetEntry(irow);
    ExitingSpec->Fill(energy);
  }

  for (int irow = 0; irow < Flux4002->GetEntries(); ++irow) {
     Flux4002->SetBranchAddress("energy",     &energy);        Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("tarcflux",   &tarcflux);      Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("errstat",    &errstat);       Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("g4flux",     &g4flux);        Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("g4perp",     &g4perp);        Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("g4fluence",  &g4fluence);     Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("g4err",      &g4error);       Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("rawflux",    &rawflux);       Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("eflux",      &eflux);         Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("g4eflux",    &tarcmeanflux);  Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("gstep",      &tarcmeanstep);  Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("gfl_cyl",    &tarccyl);       Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("g4front",    &tarcfrontflux); Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("g4shell",    &g4shell);       Flux4002->GetEntry(irow);
     Flux4002->SetBranchAddress("g4shellerr", &g4shellerr);    Flux4002->GetEntry(irow);
     
     TARCDataFluenceHi->Fill(energy, tarcflux);
     TAxis* xAxis1 = TARCDataFluenceHi->GetXaxis();
     Int_t binX1 = xAxis1->FindBin(energy);
     TARCDataFluenceHi->SetBinError(binX1, errstat);
     
     TARCDataFluenceHiErr->Fill(energy, tarcflux);
     TAxis* xAxis2 = TARCDataFluenceHiErr->GetXaxis();
     Int_t binX2 = xAxis2->FindBin(energy);
     TARCDataFluenceHiErr->SetBinError(binX2, errstat);

     double corr = g4shell;
     TARCG4FluenceHi->Fill(energy, corr);
     TAxis* xAxis3 = TARCG4FluenceHi->GetXaxis();
     Int_t binX3 = xAxis3->FindBin(energy);
     TARCG4FluenceHi->SetBinError(binX3, g4error);
     if (tarcflux != 0.0) {
       double ratio = corr / tarcflux;
       double ratrat = ratio/(1.0-ratio);
       // std::cout << energy << "    " << ratio << std::endl;
       TARCG4RatRatHi->Fill(energy, ratrat);
       TARCG4RatioHi->Fill(energy, ratio);
     }
  }

  for (int irow = 0; irow < Flux4004->GetEntries(); ++irow) {
    Flux4004->SetBranchAddress("energy",     &energy);        Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("tarcflux",   &tarcflux);      Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("errstat",    &errstat);       Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("g4flux",     &g4flux);        Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("g4perp",     &g4perp);        Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("g4fluence",  &g4fluence);     Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("g4err",      &g4error);       Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("rawflux",    &rawflux);       Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("gstep",      &tarcmeanstep);  Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("gfl_cyl",    &tarccyl);       Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("g4front",    &tarcfrontflux); Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("g4shell",    &g4shell);       Flux4004->GetEntry(irow);
    Flux4004->SetBranchAddress("g4shellerr", &g4shellerr);    Flux4004->GetEntry(irow);
   
    TARCDataFluenceHe3->Fill(energy, tarcflux);
    TAxis* xAxis1 = TARCDataFluenceHe3->GetXaxis();
    Int_t binX1 = xAxis1->FindBin(energy);
    TARCDataFluenceHe3->SetBinError(binX1, errstat);

    TARCDataFluenceHe3Err->Fill(energy, tarcflux);
    TAxis* xAxis2 = TARCDataFluenceHe3Err->GetXaxis();
    Int_t binX2 = xAxis2->FindBin(energy);
    TARCDataFluenceHe3Err->SetBinError(binX2, errstat);
  
    double corr = g4shell;
    TARCG4FluenceHe3->Fill(energy, corr);
    TAxis* xAxis3 = TARCG4FluenceHe3->GetXaxis();
    Int_t binX3 = xAxis3->FindBin(energy);
    TARCG4FluenceHe3->SetBinError(binX3, g4error);
    if (tarcflux != 0.0) {
      double ratio = corr / tarcflux;
      double ratrat = ratio/(1.0-ratio);
      TARCG4RatRatHe3->Fill(energy, ratrat);
      TARCG4RatioHe3->Fill(energy, ratio);
    }
  }

  for (int irow = 0; irow < Flux4005->GetEntries(); ++irow) {
    Flux4005->SetBranchAddress("energy",  &energy);          Flux4005->GetEntry(irow);
    Flux4005->SetBranchAddress("tarcflux",  &tarcflux);      Flux4005->GetEntry(irow);
    Flux4005->SetBranchAddress("errstat",   &errstat);       Flux4005->GetEntry(irow);
    Flux4005->SetBranchAddress("g4flux",    &g4flux);        Flux4005->GetEntry(irow);
    Flux4005->SetBranchAddress("g4perp",    &g4perp);        Flux4005->GetEntry(irow);
    Flux4005->SetBranchAddress("g4fluence", &g4fluence);     Flux4005->GetEntry(irow);
    Flux4005->SetBranchAddress("g4err",     &g4error);       Flux4005->GetEntry(irow);
    Flux4005->SetBranchAddress("rawflux",   &rawflux);       Flux4005->GetEntry(irow);
    Flux4005->SetBranchAddress("gstep",     &tarcmeanstep);  Flux4005->GetEntry(irow);
    Flux4005->SetBranchAddress("gfl_cyl",   &tarccyl);       Flux4005->GetEntry(irow);
    Flux4005->SetBranchAddress("g4front",   &tarcfrontflux); Flux4002->GetEntry(irow);
    Flux4005->SetBranchAddress("g4shell",   &g4shell);       Flux4005->GetEntry(irow);
     
     
    TARCDataFluenceLi->Fill(energy, tarcflux);
    TAxis* xAxis1 = TARCDataFluenceLi->GetXaxis();
    Int_t binX1 = xAxis1->FindBin(energy);
    TARCDataFluenceLi->SetBinError(binX1, errstat);
     
    TARCDataFluenceLiErr->Fill(energy, tarcflux);
    TAxis* xAxis2 = TARCDataFluenceLiErr->GetXaxis();
    Int_t binX2 = xAxis2->FindBin(energy);
    TARCDataFluenceLiErr->SetBinError(binX2, errstat);
 
    double corr = g4shell;
    TARCG4FluenceLi->Fill(energy, corr);
    TAxis* xAxis3 = TARCG4FluenceLi->GetXaxis();
    Int_t binX3 = xAxis3->FindBin(energy);
    TARCG4FluenceLi->SetBinError(binX3, g4error);
    if (tarcflux != 0.0) {
      double ratio = corr / tarcflux;
      double ratrat = ratio/(1.0-ratio);
      TARCG4RatRatLi->Fill(energy, ratrat);
      TARCG4RatioLi->Fill(energy, ratio);
    }
  }
     

  // Draw on the Canvas
  //******************************** c0 *******************************************
  TCanvas* c0 = new TCanvas("c0","", 1020, 800);
  c0->Divide(2, 2);
  c0->cd(1);
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
  TARCG4FluenceLi->SetTitle("G4 Flux Li Data");
  TARCG4FluenceLi->SetMarkerStyle(3);
  TARCG4FluenceLi->GetXaxis()->SetTitle("Energy / eV");
  TARCG4FluenceLi->GetYaxis()->SetTitle("Flux ( dN/dE / source gamma)");
  TARCG4FluenceLi->GetYaxis()->SetTitleOffset(1.4);
  TARCG4FluenceLi->SetLineColor(kBlue - 2);
  TARCG4FluenceLi->Draw("SAME E1");
      
  c0->cd(3); // Left Bottom gPad
  gPad->SetLogx();
  gPad->SetLogy();
  gStyle->SetHistLineWidth(3);
   
  tlx=new TLatex(0.12, 0.93, "TARC fluence comparison  for 1.5 GeV/c protons ");
  tlx->SetTextSize(0.5);
  tlx->SetNDC(kTRUE);
  tlx->Draw();
  gPad->DrawFrame(1.0e-3, 1.0e5, 5.0e7, 1.0e10,"; Energy/eV; EdF/dE n/cm^{2}/10^{9}p")->GetXaxis()->SetTitleOffset(1.2);
      
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
  gPad->SetLogy(); 
  NeutEnergy->SetTitle("Neutron Deposition");
  NeutEnergy->SetLineColor(kRed);
  NeutEnergy->SetMarkerColor(kRed);
  NeutEnergy->Draw("SAME COLZ");
  NeutEnergy->GetXaxis()->SetTitle("log10(Energy (eV))");
  NeutEnergy->GetYaxis()->SetTitle("log10(Flux(dN/dE))");
  NeutEnergy->GetXaxis()->SetTitleSize(0.03);
  NeutEnergy->GetXaxis()->SetTitleOffset(1.2);
  //NeutEnergy->GetYaxis()->SetTitle("Neutrons");
  NeutEnergy->GetYaxis()->SetTitleOffset(1.3);
     
  c0->Print("TARC_Report_Summary.png");
  c0->Close();

  //******************************** c1 *******************************************
  TCanvas*  c1 = new TCanvas("c1", "TARC Study", 900, 700);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleX(0.2);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->DrawFrame(1.0e-3, 5.0e5, 5.0e6, 1.0e10,"; log10(Energy/eV); EdF/dE n/cm^{2}/10^{9}p")->GetXaxis()->SetTitleOffset(1.2);
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
   
  TString savedFileN = thisTitlePart10 + "_fluence.png";
  c1->Print(savedFileN);
  c1->Close();

  //******************************** c2 *******************************************
  TCanvas*  c2 = new TCanvas("c2", "TARC Ratio G4/Data Study", 900, 700);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleX(0.2);
  gPad->SetLogx();
  gPad->SetLogy();
  tlx = new TLatex(0.15, 0.93, "TARC G4/Data Ratio comparison for 1.5 GeV/c proton");    // plotTitle1);
  tlx->SetTextSize(0.04);
  tlx->SetNDC(kTRUE);
  tlx->Draw();
   
  gPad->DrawFrame(5.0e-3, 1.0e-17, 5e6, 5.0e-11,"; log10(Energy/eV); G4 Shell / Data")->GetXaxis()->SetTitleOffset(1.2);
  //gPad->DrawFrame(4e4, 1.0e-3, 5e6, 5.0e-1,"; Energy/eV; G4/Data")->GetXaxis()->SetTitleOffset(1.2);  // for the Hi Data
  TARCG4RatioHi->SetMarkerColor(kRed);
  TARCG4RatioHi->SetMarkerStyle(kFullCircle);
  TARCG4RatioHi->SetLineColor(kRed);
  TARCG4RatioHi->SetStats(kTRUE);
  TARCG4RatioHi->Draw("SAME P");
   
  TARCG4RatioHe3->SetMarkerColor(kBlue);
  TARCG4RatioHe3->SetMarkerStyle(kFullSquare);
  TARCG4RatioHe3->SetLineColor(kBlue);
  TARCG4RatioHe3->SetStats(kTRUE);
  TARCG4RatioHe3->Draw(" SAME PLC");// PMC");
   
  TARCG4RatioLi->SetMarkerColor(kGreen);
  TARCG4RatioLi->SetMarkerStyle(kFullTriangleUp);
  TARCG4RatioLi->SetLineColor(kGreen);
  TARCG4RatioLi->SetStats(kTRUE);
  TARCG4RatioLi->Draw(" SAME PLC");// PMC");
   
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

  //******************************** c3 *******************************************  
  TCanvas*  c3 = new TCanvas("c3", "TARC Gamma Energy Deposition Study", 900, 500);
  c3->Divide(2, 1);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleX(0.2);
  GammaED->SetMarkerStyle(7);
  GammaED->SetMarkerColor(kBlue);
  GammaED->SetLineColor(kBlue);
  GammaED->GetXaxis()->SetTitle("Energy (eV)");
  GammaED->GetYaxis()->SetTitle("Gamma");
  GammaED->GetYaxis()->SetTitleOffset(1.4);
  GammaED->SetTitle("Gamma Deposition");
   
  c3->cd(1);
  GammaED->Draw();
  c3->cd(2);
  gPad->SetLogx();
  TH1* h1 = GammaED->DrawCopy();
  h1->GetXaxis()->SetTitle("log10(Energy (eV))");
  c1->cd(0);
  savedFileN = thisTitlePart10 + "_Gamma_Edep.png";
  c3->Print(savedFileN);
  c3->Close();
  

}