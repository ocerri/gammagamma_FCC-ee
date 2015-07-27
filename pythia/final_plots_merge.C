#include <Riostream.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "TPythia6.h"
#include "TMCParticle.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TError.h"
#include "TMath.h"
//#include "TTree.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TDatime.h"
#include "TLegend.h"
#include "TColor.h"


//Create histos-----------------------------------------------------
TH1F* h_dsdeta10 = 0;
TH1F* h_dsdeta30 = 0;

TH1F* h_dsdpt10 = 0;
TH1F* h_dsdpt30 = 0;

///Start merging procedure------------------------------------------
void final_plots_merge(){

 TFile* file10_dsdpt = new TFile("_root/MSTP14-10_final_out.root");
 h_dsdpt10 = (TH1F*)file10_dsdpt->Get("hdsigmadpT");
 h_dsdpt10->SetMarkerColor(2);
 h_dsdpt10->SetLineColor(2);
 h_dsdpt10->SetTitle("PYTHIA simulation of d#sigma/dp_{t} for e^{+}e^{-} #rightarrow h^{#pm} at 160 GeV");
 
 
 TFile* file30_dsdpt = new TFile("_root/MSTP14-30_final_out.root");
 h_dsdpt30 = (TH1F*)file30_dsdpt->Get("hdsigmadpT");
 h_dsdpt30->SetLineColor(4);
 //h_dsdpt30->SetLineStyle(4);
 h_dsdpt30->SetTitle("PYTHIA simulation of d#sigma/dp_{t} for {e}^{+} {e}^{-} -> {h}^{#pm} at 160 GeV");


 //Draw on Canvas---------------------------------------------------
 TCanvas* c_dsdpt = new TCanvas("c_dsdpt","dsigmadipt",700,700);
 //gStyle->SetOptStat(0);
 c_dsdpt->cd();

 h_dsdpt10->Draw();
 h_dsdpt30->Draw("same");
 c_dsdpt->SetLogy();
 c_dsdpt->SetGridx();
 c_dsdpt->SetGridy();


 TLegend* leg = c_dsdpt->BuildLegend();
 // TLegend* leg = new TLegend(0.7,0.7,0.48,0.9);
 // leg->AddEntry("h_dsdpt10" ,"MSTP(14)=10");
 // leg->AddEntry("h_dsdpt30" ,"MSTP(14)=30");
 leg->Draw();
 



 
}
