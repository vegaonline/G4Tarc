#include <string.h>
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TPad.h"
#include "Riostream.h"

void neutron() {
	TString dir = "./";
	dir.ReplaceAll("testPlotneutron.C","");
	dir.ReplaceAll("/./","/");
	printf("%s", dir.Data());
	ifstream in;
	in.open(Form("%snspectra.dat", dir.Data()));

	Float_t a, b, c, d;
	TFile *f = new TFile("nspect.root", "RECREATE");
	TCanvas *c1 = new TCanvas("test", "h1");
	TH1F *h1 = new TH1F("h1", "a dist.", 50, 0, 100);
	TTree *T = new TTree("nData", "data from ascii file");
	Long64_t nlines = T->ReadFile(Form("%sneutronSpectra.dat", dir.Data()), "a:b:c:d");

	printf(" found %lld points\n", nlines);
	T->Draw("a","b>=0");
	h1->Draw();
	c1->Update();
	in.close();
	T->Write();

}
