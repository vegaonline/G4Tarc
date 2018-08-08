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

	TFile *f = new TFile("neutSpect.root", "RECREATE");
	TH1F *h1 = new TH1F("h1", "a dist.", 200, -5, 8);
	TTree *T = new TTree("nData", "data from ascii file");
	Long64_t nlines = T->ReadFile(Form("%sneutSpec.dat", dir.Data()), "a:b");

	printf(" found %lld points\n", nlines);
	T->Draw("a","b");
	T->Write();

}
