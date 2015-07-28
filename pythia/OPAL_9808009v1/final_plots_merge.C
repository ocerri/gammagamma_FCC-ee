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
#include "TPaveText.h"
#include "TPaveStats.h"



//Create histos-----------------------------------------------------
TH1F* h_dsdeta10 = 0;
TH1F* h_dsdeta30 = 0;
TH1F* h_dsdeta_exp = 0;

TH1F* h_dsdpt10 = 0;
TH1F* h_dsdpt30 = 0;

///Start merging procedure------------------------------------------
void final_plots_merge(){

 TFile* file10_dsdpt = new TFile("_root/MSTP14-10_final_out.root");
 h_dsdpt10 = (TH1F*)file10_dsdpt->Get("hdsigmadpT");
 h_dsdpt10->SetMarkerColor(2);
 h_dsdpt10->SetLineColor(2);
 h_dsdpt10->SetAxisRange(0.01,50000,"Y");
 h_dsdpt10->SetTitle("PYTHIA simulation of d#sigma/dp_{t} for e^{+}e^{-} #rightarrow h^{#pm} at 160 GeV");
 
 
 TFile* file30_dsdpt = new TFile("_root/MSTP14-30_final_out.root");
 h_dsdpt30 = (TH1F*)file30_dsdpt->Get("hdsigmadpT");
 h_dsdpt30->SetLineColor(4);
 h_dsdpt30->SetMarkerColor(4);
 h_dsdpt10->SetAxisRange(0.01,50000,"Y");
 h_dsdpt30->SetTitle("PYTHIA simulation of d#sigma/dp_{t} for {e}^{+} {e}^{-} -> {h}^{#pm} at 160 GeV");

 // TFile* file10_dsdeta = new TFile("_root/MSTP14-10_final.root");
 // h_dsdeta10 = (TH1F*)file10_dsdeta->Get("hdsigmadeta");
 // h_dsdeta10->SetMarkerColor(2);
 // h_dsdeta10->SetLineColor(2);
 // h_dsdeta10->SetAxisRange(8.,24.,"Y");
 // h_dsdeta10->SetTitle("PYTHIA simulation of d#sigma/deta for e^{+}e^{-} #rightarrow h^{#pm} at 160 GeV");
 
 
 // TFile* file30_dsdeta = new TFile("_root/MSTP14-30_final.root");
 // h_dsdeta30 = (TH1F*)file30_dsdeta->Get("hdsigmadeta");
 // h_dsdeta30->SetMarkerColor(4);
 // //h_dsdeta30->SetLineColor(4);
 // h_dsdeta30->SetAxisRange(8.,24.,"Y");
 // h_dsdeta30->SetTitle("PYTHIA simulation of d#sigma/deta for {e}^{+} {e}^{-} -> {h}^{#pm} at 160 GeV");

 //Plotting experimental data from txt-------------------------------------------

 std::ifstream myfile;
 myfile.open("_txt/dsdpt_opal.txt");

 const Int_t n_data = 21; //number of data (rows) in .txt

 Double_t x[n_data];
 Double_t y[n_data];
 Double_t ex[n_data];
 Double_t ey[n_data];

 int j=0;

 //Read data from file .txt
 while(true) {
   myfile >> x[j];
   myfile >> ex[j];
   myfile >> y[j];
   myfile >> ey[j];
   if( myfile.eof() ) break;
   ++j;
 }

 myfile.close();

 TGraphErrors* gr = new TGraphErrors(n_data,x,y,ex,ey); //Create TGraph
 
 gr->SetName("gr");
 //gr->SetLineWidth(2);
 gr->SetMarkerStyle(5);
 gr->SetMarkerSize(1);
 gr->SetMarkerColor(1);

 

 //Draw on Canvas---------------------------------------------------
   
 TCanvas* c_dsdpt = new TCanvas("c_dsdpt","dsigmadipt",700,700);
 gStyle->SetOptStat(0);
 c_dsdpt->cd();

 h_dsdpt10->Draw();
 h_dsdpt30->Draw("same");
 gr->Draw("P");
 c_dsdpt->SetLogy();
 c_dsdpt->SetGridx();
 c_dsdpt->SetGridy();

 // TCanvas* c_dsdeta = new TCanvas("c_dsdeta","dsigmadeta",700,700);
 // //gStyle->SetOptStat(0);
 // c_dsdeta->cd();

 // h_dsdeta10->Draw();
 // h_dsdeta30->Draw("same");
 // gr->Draw("P");
 // c_dsdeta->SetGridx();
 // c_dsdeta->SetGridy();

 //Building Legend and adding comments-----------------------------------------------------

 TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
 //leg->SetHeader("Legend");
 leg->AddEntry(h_dsdpt10,"Histogram MSTP(14)=10","lep");
 leg->AddEntry(h_dsdpt30,"Histogram MSTP(14)=30","lep");
 leg->AddEntry("gr","Experimental data from OPAL","lep");
 leg->Draw();

 TLegend* leg2 =new TLegend(0.1,0.7,0.5,0.9,"#splitline{#cbar#eta#cbar< 1.5}{10<W<125 GeV}");
 leg2->Draw();
 
}
