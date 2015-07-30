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
#include "TGraphAsymmErrors.h"



//Create histos-----------------------------------------------------

TH1F* h_dsdeta90 = 0;
TH1F* h_dsdeta160 = 0;
TH1F* h_dsdeta240 = 0;
TH1F* h_dsdeta350 = 0;

///Start merging procedure------------------------------------------
void final_plots_merge_dNdeta(){
  // for _tot analysis-----------------------------
  Double_t xsection[]={8294.51,13027.4,17500.8,21962.6};
  Double_t ntrials[]={307970/25000,351324/25000,380337/25000,414176/25000};
  Double_t total_hadrons[4];

 TFile* file90_dsdeta = new TFile("90GeV/_root/pythia6_gammagamma_hadrons_90GeV_seed479395016_Nevts25000.root");
 h_dsdeta90 = (TH1F*)file90_dsdeta->Get("hdsigmadeta");
 h_dsdeta90->SetMarkerStyle(24);
 h_dsdeta90->SetMarkerColor(1);
 h_dsdeta90->SetLineColor(1);
 h_dsdeta90->SetAxisRange(0.001,1200,"Y");
 h_dsdeta90->SetXTitle("#eta");
 h_dsdeta90->SetYTitle("dN/d#eta");
 h_dsdeta90->SetTitle("#splitline{PYTHIA simulation of dN/d#eta for e^{+}e^{-}#rightarrow e^{+}e^{-}/#gamma*#gamma* #rightarrow h}{Total event simulated: 25k}");
 h_dsdeta90->Scale(ntrials[0]/xsection[0]);
 total_hadrons[0]=h_dsdeta90->Integral();

 TFile* file160_dsdeta = new TFile("160GeV/_root/pythia6_gammagamma_hadrons_160GeV_seed479395646_Nevts25000.root");
 h_dsdeta160 = (TH1F*)file160_dsdeta->Get("hdsigmadeta");
 h_dsdeta160->SetMarkerStyle(25);
 h_dsdeta160->SetMarkerColor(2);
 h_dsdeta160->SetLineColor(2);
 h_dsdeta160->SetAxisRange(0.001,1200,"Y");
 h_dsdeta160->SetTitle("PYTHIA simulation of d#N/#eta for e^{+}e^{-}#rightarrow e^{+}e^{-}/#gamma*#gamma* #rightarrow h^{#pm}");
 h_dsdeta160->Scale(ntrials[1]/xsection[1]);
 total_hadrons[1]=h_dsdeta160->Integral();

 TFile* file240_dsdeta = new TFile("240GeV/_root/pythia6_gammagamma_hadrons_240GeV_seed479395489_Nevts25000.root");
 h_dsdeta240 = (TH1F*)file240_dsdeta->Get("hdsigmadeta");
 h_dsdeta240->SetMarkerStyle(26);
 h_dsdeta240->SetMarkerColor(3);
 h_dsdeta240->SetLineColor(3);
 h_dsdeta240->SetAxisRange(0.001,1200,"Y");
 h_dsdeta240->SetTitle("PYTHIA simulation of d#sigma/#eta for e^{+}e^{-}#rightarrow e^{+}e^{-}/#gamma*#gamma* #rightarrow h^{#pm}");
 h_dsdeta240->Scale(ntrials[2]/xsection[2]);
 total_hadrons[2]=h_dsdeta240->Integral();

 TFile* file350_dsdeta = new TFile("350GeV/_root/pythia6_gammagamma_hadrons_350GeV_seed479395808_Nevts25000.root");
 h_dsdeta350 = (TH1F*)file350_dsdeta->Get("hdsigmadeta");
 h_dsdeta350->SetMarkerStyle(27);
 h_dsdeta350->SetMarkerColor(4);
 h_dsdeta350->SetLineColor(4);
 h_dsdeta350->SetAxisRange(0.001,1200,"Y");
 h_dsdeta350->SetXTitle("#eta");
 h_dsdeta350->SetYTitle("dN/d#eta");
 h_dsdeta350->SetTitle("PYTHIA simulation of dN/d#eta for e^{+}e^{-}#rightarrow e^{+}e^{-}/#gamma*#gamma* #rightarrow h^{#pm}");
 h_dsdeta350->Scale(ntrials[3]/xsection[3]);
 total_hadrons[3]=h_dsdeta350->Integral();

 

 //Draw on Canvas---------------------------------------------------
   
 TCanvas* c_dsdeta = new TCanvas("c_dsdeta","dsdeta",700,700);
 //gStyle->SetOptStat(0);
 c_dsdeta->cd();

 h_dsdeta350->Draw();
 h_dsdeta90->Draw("same");
 h_dsdeta160->Draw("same");
 h_dsdeta240->Draw("same");
 
 c_dsdeta->SetGridx();
 c_dsdeta->SetGridy();

 //Building Legend-----------------------------------------------------

 TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
 //leg->SetHeader("Legend");
 leg->AddEntry(h_dsdeta90,Form("#sqrt{s}=90,  %.2f total hadrons",total_hadrons[0]),"lep");
 leg->AddEntry(h_dsdeta160,Form("#sqrt{s}=160,  %.2f total hadrons",total_hadrons[1]),"lep");
 leg->AddEntry(h_dsdeta240,Form("#sqrt{s}=240,  %.2f total hadrons",total_hadrons[2]),"lep");
 leg->AddEntry(h_dsdeta350,Form("#sqrt{s}=350,  %.2f total hadrons",total_hadrons[3]),"lep");
 leg->Draw();
cout <<"ciao" << endl;

//TLegend *leg2 =new TLegend(0.5945559,0.7688889,0.8968481,0.8859259,"#splitline{|#eta|<1.5}{5 < W < 35 GeV}");
//leg2->Draw();


 // the following line is needed to avoid that the automatic redrawing of stats
 h_dsdeta90->SetStats(0);
 h_dsdeta160->SetStats(0);
 h_dsdeta240->SetStats(0);
 h_dsdeta350->SetStats(0);

 c_dsdeta->Modified();


}
