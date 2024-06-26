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
  TFile* tf1 = TFile::Open("G4TARC_output.root");

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
  TNtuple* h1SecTuple      = (TNtuple*) tf1->Get("h0_Secondary");
  TNtuple* h2NET           = (TNtuple*) tf1->Get("h1_N_ET");
  TNtuple* h3NExiting      = (TNtuple*) tf1->Get("h2_N_Exiting");
  TNtuple* Flux4002        = (TNtuple*) tf1->Get("h3_Flux_4002");
  TNtuple* Flux4004        = (TNtuple*) tf1->Get("h4_Flux_4004");
  TNtuple* Flux4005        = (TNtuple*) tf1->Get("h5_Flux_4005"); 
  TNtuple* h8SpallN        = (TNtuple*) tf1->Get("h7_Created_N");
  TNtuple* RadShellFluence = (TNtuple*) tf1->Get("h8_Rad_Shell_Fluence");  
  TNtuple* ExptLiData      = (TNtuple*) tf1->Get("h9_Rad_Fluence_Expt_Li_Data");
  TNtuple* ExptHe3Data     = (TNtuple*) tf1->Get("h10_Rad_Fluence_Expt_He3_Data");
  TNtuple* h15OtherET      = (TNtuple*) tf1->Get("h12_Other_ET");


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
  //  TH1F* TARCG4RatRatHi         = new TH1F("HiRatioBYRatio", "TARC Fluence Ratio G4 Fluence (High)", (xHiCnt - 1), fXbinHi);
  TH1F* TARCG4RatioHe3         = new TH1F("RatioG4He", "TARC Fluence Ratio G4 Fluence (He3)", (xLoCnt - 1), fXbinLo);
  //  TH1F* TARCG4RatRatHe3        = new TH1F("He3RatioBYRatio", "TARC Fluence Ratio G4 Fluence (He3)", (xLoCnt - 1), fXbinLo);
  TH1F* TARCG4RatioLi          = new TH1F("RatioG4Li", "TARC Fluence Ratio G4 Fluence (Li)", (xLoCnt - 1), fXbinLo);
  //  TH1F* TARCG4RatRatLi         = new TH1F("LiRatioBYRatio", "TARC Fluence Ratio G4 Fluence (Li)", (xLoCnt - 1), fXbinLo);
  TH2F* radialHisto            = new TH2F("Radial", "TARC radial", 1000, -1000, 1000, 1000, 1.0, 2.5e+7);
  TH2F* TARCNeutCorr           = new TH2F("Correlation", "TARC Correlation function", 5000, 0.52, 50, 5000,0.00057, 1.8e9);

  // Variables for histograms
  double energy, time, tarcflux, errstat, g4flux, g4perp, g4fluence, g4error, rawflux, eflux, tarcmeanflux, abseflux;
  double radius, g4shell, g4shellerr, tarcmeanstep, tarccyl, tarcfrontflux;
  double t0Corr = 0.37; //us
  double KCorr = 173.3; // keV X us^2
  TLatex* tlx;
  double maxNewX = -99999.99, minNewX = -1.0 * maxNewX;
  double maxY = maxNewX, minY = minNewX;
  
  
  // Fill the histograms
  for (int irow = 0; irow < h3NExiting->GetEntries(); ++irow){
    h3NExiting->SetBranchAddress("energy", &energy);
    h3NExiting->GetEntry(irow);
    ExitingSpec->Fill(energy);
  }

  //std::ofstream ETout ("testET.dat", std::ios::out);
  /*
  for (int irow = 0; irow < h2NET->GetEntries(); irow++) {
    h2NET->SetBranchAddress("energy", &energy);
    h2NET->SetBranchAddress("time", &time); // this time is already in microsecond
    h2NET->GetEntry(irow);
    double newX = (time + t0Corr) * sqrt(energy/1000.0);
    // maxNewX = std::max(maxNewX, newX);
    // minNewX = std::min(minNewX, newX);
    // maxY = std::max(maxY, energy);
    // minY = std::min(minY, energy);
    // ETout << newX << "   " << energy << std::endl;
    TARCNeutCorr->Fill(newX, energy);
  }
  */
  //cout << "   maxX: " << maxNewX << "  minX: " << minNewX << "  maxY: " << maxY << "  minY: " << minY << endl;
  //ETout.close();

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
       //double ratrat = ratio/(1.0-ratio);
       //std::cout << energy << "    " << ratio << std::endl;
       //TARCG4RatRatHi->Fill(energy, ratrat);
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
      //double ratrat = ratio/(1.0-ratio);
      //TARCG4RatRatHe3->Fill(energy, ratrat);
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
      //double ratrat = ratio/(1.0-ratio);
      //TARCG4RatRatLi->Fill(energy, ratrat);
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
  TARCG4FluenceHi->SetTitle("G4 based Neutron Flux");
  TARCG4FluenceHi->GetXaxis()->SetTitle("Energy / eV");
  TARCG4FluenceHi->GetYaxis()->SetTitle("Neutron Flux ( dN/dE / source gamma)");
  TARCG4FluenceHi->GetYaxis()->SetTitleOffset(1.2);
  TARCG4FluenceHi->SetMarkerStyle(3);
  TARCG4FluenceHi->SetMarkerColor(kBlue + 3);
  TARCG4FluenceHi->Draw("SAME");
        
  c0->cd(2); // Right Top gPad
  gPad->SetLogx();
  gPad->SetLogy();
  TARCG4FluenceLi->SetTitle("Neutron Flux at Li detector Data");
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
  //tlx->Draw();
  // gPad->DrawFrame(1.0e-3, 1.0e5, 1.0e7, 1.0e10,"; Energy/eV; EdF/dE n/cm^{2}/10^{9}p")->GetXaxis()->SetTitleOffset(1.2);
      
  TARCDataFluenceHi->SetMarkerStyle(kFullCircle);
  TARCDataFluenceHi->SetMarkerSize(0.8);
  TARCDataFluenceHi->SetMarkerColor(kBlue );
  TARCDataFluenceHi->Draw("SAME E1");
      
  TARCDataFluenceHe3->SetMarkerStyle(kOpenSquare);
  TARCDataFluenceHe3->SetMarkerSize(0.7);
  TARCDataFluenceHe3->SetMarkerColor(kRed );
  //  TARCDataFluenceHe3->Draw("SAME E1");
      
  TARCDataFluenceLi->SetMarkerStyle(kDiamond);
  TARCDataFluenceLi->SetMarkerColor(kYellow -1);
  //TARCDataFluenceLi->Draw("SAME E1");
   
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
  NeutEnergy->SetTitle("Neutron Deposition");
  NeutEnergy->SetLineColor(kRed);
  NeutEnergy->SetMarkerColor(kRed);
  NeutEnergy->Draw("SAME COLZ");
  NeutEnergy->GetXaxis()->SetTitle("log10(Energy (eV))");
  NeutEnergy->GetYaxis()->SetTitle("Flux(dN/dE)");
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
   
  gPad->DrawFrame(5.0e-3, 1.0e-9, 5e6, 5.0e-4,"; log10(Energy/eV); G4 Shell / Data")->GetXaxis()->SetTitleOffset(1.2);
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
  GammaED->GetXaxis()->SetTitle("log10(Energy (eV))");
  GammaED->GetYaxis()->SetTitle("Gamma");
  GammaED->GetYaxis()->SetTitleOffset(1.4);
  GammaED->SetTitle("Gamma Deposition");
  //std::cout << GammaED->GetMinimum() << "   " << GammaED->GetMaximum() << std::endl;
   
  c3->cd(1);
  gPad->SetLogx();
  GammaED->Draw();
  c3->cd(2);
  gPad->SetLogx();
  gPad->SetLogy();
  TH1* h1 = GammaED->DrawCopy();
  h1->GetXaxis()->SetTitle("log10(Energy (eV))");
  c1->cd(0);
  savedFileN = thisTitlePart10 + "_Gamma_Edep.png";
  c3->Print(savedFileN);
  c3->Close();
  
  //******************************** c4 *******************************************
  TCanvas*  c4 = new TCanvas("c4", "TARC Neutron Energy Deposition Study", 900, 500);
  c4->Divide(2, 1);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleX(0.2);
  NeutEnergy->SetMarkerStyle(7);
  NeutEnergy->SetMarkerColor(kBlue);
  NeutEnergy->SetLineColor(kBlue);
  NeutEnergy->GetXaxis()->SetTitle("log10(Energy (eV))");
  NeutEnergy->GetYaxis()->SetTitle("Neutron ");
  NeutEnergy->GetXaxis()->SetTitleOffset(1.6);
  NeutEnergy->GetYaxis()->SetTitleOffset(1.4);
  NeutEnergy->SetTitle("Neutron energy Deposition");
  c4->cd(1);
  gPad->SetLogx();
  NeutEnergy->Draw();
  c4->cd(2);
  gPad->SetLogx();
  gPad->SetLogy();
  NeutEnergy->DrawCopy();
  h1->GetXaxis()->SetTitle("log10(Energy (eV))");
  c1->cd(0);
  savedFileN = thisTitlePart10 + "_Neutron_Edep.png";
  c4->Print(savedFileN);
  c4->Close();

  //******************************** c5 *******************************************  
  TCanvas*  c5 = new TCanvas("c5", "TARC Electron Energy Deposition Study", 1200, 500);
  c5->Divide(2, 1);
  gStyle->SetHistLineWidth(3);  
  gStyle->SetTitleX(0.5);
  ElecED->GetXaxis()->SetTitle("Energy (eV)");
  ElecED->GetYaxis()->SetTitle("Electron");
  ElecED->SetMarkerStyle(7);
  ElecED->SetMarkerColor(kBlue);
  ElecED->SetLineColor(kBlue);
  ElecED->GetYaxis()->SetTitleOffset(1.3);
  ElecED->SetTitle("Electron Deposition");
  //std::cout << ElecED->GetMinimum() << "  " << ElecED->GetMaximum() << std::endl;
  c5->cd(1);
  ElecED->Draw();
  c5->cd(2);
  gPad->SetLogx();
  TH1* h = ElecED->DrawCopy();
  h->GetXaxis()->SetTitle("log10(Energy / eV)");
  h->GetYaxis()->SetTitleOffset(1.3);
  c1->cd(0);
  savedFileN = thisTitlePart10 + "_Electron_Edep.png";
  c5->Print(savedFileN);
  c5->Close();

  //******************************** c6 *******************************************  
  TCanvas*  c6 = new TCanvas("c6", "TARC Positron Energy Deposition Study", 700, 500);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleX(0.2);
  //gPad->SetLogx();
  gPad->SetLogy();
  PositED->SetTitle("Positron Deposition");
  PositED->SetMarkerStyle(7);
  PositED->SetMarkerColor(kBlue);
  PositED->SetLineColor(kBlue);
  PositED->GetXaxis()->SetTitle("Energy (eV)");
  PositED->Draw();
  savedFileN = thisTitlePart10 + "_Positron_Edep.png";
  c6->Print(savedFileN);
  c6->Close();

  //******************************** c7 *******************************************  
  TCanvas*  c7 = new TCanvas("c7", "TARC Neutron Energy-Time Study", 900, 600);
  c7->Divide(2, 1);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleX(0.2);
  NeutronET->SetMarkerStyle(7);
  NeutronET->SetMarkerColor(kBlue);
  NeutronET->SetLineColor(kBlue);
  NeutronET->SetTitle("Neutron Energy - Time characteristics");
  NeutronET->GetXaxis()->SetTitle("log10(Time (#mus))");
  NeutronET->GetYaxis()->SetTitle("log10(Neutron Energy (eV))");
  c7->cd(1);
  NeutronET->Draw();
  c7->cd(2);
  NeutronET->Draw("CONT1Z");
  c7->cd(0);
  savedFileN = thisTitlePart10 + "_Neutron_Energy_Time.png";
  c7->Print(savedFileN);
  c7->Close();

  //******************************** c8 *******************************************
  TCanvas*  c8 = new TCanvas("c8", "TARC Other Particles Energy-Time Study", 1200, 500);
  c8->Divide(2, 1);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleX(0.2);
  OtherET->SetMarkerStyle(7);
  OtherET->SetMarkerColor(kBlue);
  OtherET->SetLineColor(kBlue);
  OtherET->SetTitle("Time - Energy characteristics for other particles.");
  OtherET->GetXaxis()->SetTitle("log10(Time (#mus))");
  OtherET->GetYaxis()->SetTitle("log10(Energy (eV))");
  c8->cd(1);
  OtherET->Draw();
  c8->cd(2);
  OtherET->Draw("CONT1Z");
  c8->cd(0);
  savedFileN = thisTitlePart10 + "_Other_Energy_Time.png";
  c8->Print(savedFileN);
  c8->Close();
 
  //******************************** c9 *******************************************
  TCanvas* c9 = new TCanvas("c9","TARC Exiting Neutron Spectrum", 1200, 700);
  c9->Divide(2,1);
  gStyle->SetHistLineWidth(3);
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
  savedFileN = thisTitlePart10 + "_Exiting_Neutron_Spectrum.png";
  c9->Print(savedFileN);
  c9->Close();


  //******************************** c10 *******************************************
  TCanvas* c10 = new TCanvas("c10", "TARC Radial Fluence", 900, 700);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleX(0.2);
  gPad->SetLogy();
  gStyle->SetTitle("fluence");
  // gPad->DrawFrame(-220, 1e-3, 220.0, 5.0e11, "; Radial Distance / cm; dF/dE (n/cm^{2}/eV/10^{9} p)")->GetXaxis()->SetTitleOffset(1.2);
  gPad->DrawFrame(-10.0, 1e-3, 220.0, 5.0e11, "; Radial Distance / cm; dF/dE (n/cm^{2}/eV/10^{9} p)")->GetXaxis()->SetTitleOffset(1.2);
  RadShellFluence->SetMarkerStyle(28); //7);
  RadShellFluence->SetMarkerColor(kRed);
  RadShellFluence->SetMarkerSize(0.8);
  //RadShellFluence->Scan();
  //RadShellFluence->Draw("radialFluenceTrueMean/radialEnergyTrueMean  : radius / 10.0", "", "SAME"); // to convert to cm^{2}
  //  RadShellFluence->Draw("fluence/energy  : radius / 10.0", "", "SAME"); // to convert to cm^{2}
  RadShellFluence->Draw("truef/truee * 100.0  : radius / 10.0", "", "SAME"); // to convert to cm^{2}
  ExptLiData->SetMarkerStyle(28);
  ExptLiData->SetMarkerSize(0.8);
  ExptLiData->SetMarkerColor(kBlue);
  //ExptLiData->Print();
  ExptLiData->Draw("data/energy : radius / 10.0", "", "SAME");   // changing to cm
  ExptHe3Data->SetMarkerStyle(28);
  ExptHe3Data->SetMarkerSize(0.8);
  ExptHe3Data->SetMarkerColor(kGreen);
  //ExptHe3Data->Print();
  ExptHe3Data->Draw("data/energy : radius / 10.0", "", "SAME");   // changing to cm
  start = 0.11, stop = 0.9;
  xwidth  = 0.28, ywidth = 0.12;
  legend = new TLegend(start, stop - ywidth, start + xwidth, stop);
  legend->SetTextFont(62);
  legend->SetHeader("Radial Fluence distribution","C"); // option "C" allows to center the header
  legend->SetTextFont(42);
  legend->AddEntry(RadShellFluence,"Distribution for TARC simulation","ep");
  legend->AddEntry(ExptHe3Data,"Distribution for He3 Data","ep");
  legend->AddEntry(ExptLiData,"Distribution for Li data","ep");   // ep for errors and points
  legend->Draw();
  savedFileN = thisTitlePart10 + "_radial.png";
  c10->Print(savedFileN);
  c10->Close();


  //******************************** c11 *******************************************
  /*
  TCanvas* c11 = new TCanvas("c11", "TARC Neutron ET correlation", 900, 700);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleX(0.2);
  //gPad->DrawFrame(0.0, 1.0, 60.0, 1.0e10, "; Radial Distance / cm; dF/dE (n/cm^{2}/eV/10^{9} p)")->GetXaxis()->SetTitleOffset(1.2);
  gPad->SetLogy();
  //gPad->SetLogx();
  TARCNeutCorr->Print();
  gStyle->SetTitle("Neutron energy-time correlation");
  TARCNeutCorr->Draw();
  savedFileN = thisTitlePart10 + "_Neutron_ET_Correlation.png";
  c11->Print(savedFileN);
  c11->Close();
  */
 //******************************** C12 *******************************************
 TCanvas* c12 = new TCanvas("c12","", 1020, 800);
 gPad->SetLogx();
 gPad->SetLogy();
 TARCG4FluenceHi->SetTitle("G4 based Neutron Flux");
 TARCG4FluenceHi->GetXaxis()->SetTitle("Energy / eV");
 TARCG4FluenceHi->GetYaxis()->SetTitle("Neutron Flux ( dF/dE (n/cm^{2}/eV/10^{9} proton))");
 TARCG4FluenceHi->GetYaxis()->SetTitleOffset(1.2);
 TARCG4FluenceHi->SetMarkerStyle(3);
 TARCG4FluenceHi->SetMarkerColor(kBlue + 3);
 TARCG4FluenceHi->Draw(""); //"SAME");
 savedFileN = thisTitlePart10 + "_G4_Neutron_Flux.png";
 c12->Print(savedFileN);
 c12->Close();
 //******************************** END ***********************************************
															     
  tf1->Close();
}
