#include <string.h>
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TPad.h"
#include "Riostream.h"

void simple(){
	TFile  *tf1 = new TFile("G4TARC_output.root");
	TCanvas *c1 = new TCanvas("c1", "");
	TPad  *npad = new TPad("npad", "", 0.01, 0.01, 0.98, 0.98);
	npad->Draw();
	tf1->ls();
}


