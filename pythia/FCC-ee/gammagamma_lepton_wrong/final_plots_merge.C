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

TH1F* h_dNdeta90 = 0;
TH1F* h_dNdeta160 = 0;
TH1F* h_dNdeta240 = 0;
TH1F* h_dNdeta350 = 0;

///Start merging procedure------------------------------------------
void final_plots_merge(){
  Double_t total_leptons[4];

 TFile* file90_dNdeta = new TFile("90GeV/_root/pythia6_gammagamma_leptons_90GeV_seed479533328_Nevts10000.root");
 h_dNdeta90 = (TH1F*)file90_dNdeta->Get("hdNdeta");
 h_dNdeta90->SetMarkerStyle(24);
 h_dNdeta90->SetMarkerColor(1);
 h_dNdeta90->SetLineColor(1);
 h_dNdeta90->SetAxisRange(0.001,0.15,"Y");
 h_dNdeta90->SetXTitle("#eta");
 h_dNdeta90->SetYTitle("dN/d#eta (per event)");
 h_dNdeta90->SetTitle("PYTHIA simulation of dN/d#eta for e^{+}e^{-}#rightarrow e^{+}e^{-}/#gamma*#gamma* #rightarrow l^{#pm}, total generated event: 10000");
 total_leptons[0]=h_dNdeta90->Integral();

 TFile* file160_dNdeta = new TFile("160GeV/_root/pythia6_gammagamma_leptons_160GeV_seed479533187_Nevts10000.root");
 h_dNdeta160 = (TH1F*)file160_dNdeta->Get("hdNdeta");
 h_dNdeta160->SetMarkerStyle(25);
 h_dNdeta160->SetMarkerColor(2);
 h_dNdeta160->SetLineColor(2);
 //h_dNdeta160->SetAxisRange(0.001,45,"Y");
 h_dNdeta160->SetTitle("PYTHIA simulation of dN/d#eta for e^{+}e^{-}#rightarrow e^{+}e^{-}/#gamma*#gamma* #rightarrow l^{#pm}");
 total_leptons[1]=h_dNdeta160->Integral();

 TFile* file240_dNdeta = new TFile("240GeV/_root/pythia6_gammagamma_leptons_240GeV_seed479533073_Nevts10000.root");
 h_dNdeta240 = (TH1F*)file240_dNdeta->Get("hdNdeta");
 h_dNdeta240->SetMarkerStyle(26);
 h_dNdeta240->SetMarkerColor(3);
 h_dNdeta240->SetLineColor(3);
 //h_dNdeta240->SetAxisRange(0.001,45,"Y");
 h_dNdeta240->SetTitle("PYTHIA simulation of dN/d#eta for e^{+}e^{-}#rightarrow e^{+}e^{-}/#gamma*#gamma* #rightarrow l^{#pm}");
 total_leptons[2]=h_dNdeta240->Integral();

 TFile* file350_dNdeta = new TFile("350GeV/_root/pythia6_gammagamma_leptons_350GeV_seed479533325_Nevts10000.root");
 h_dNdeta350 = (TH1F*)file350_dNdeta->Get("hdNdeta");
 h_dNdeta350->SetMarkerStyle(27);
 h_dNdeta350->SetMarkerColor(4);
 h_dNdeta350->SetLineColor(4);
 //h_dNdeta350->SetAxisRange(0.001,45,"Y");
 h_dNdeta350->SetTitle("PYTHIA simulation of dN/d#eta for e^{+}e^{-}#rightarrow e^{+}e^{-}/#gamma*#gamma* #rightarrow l^{#pm}");
 total_leptons[3]=h_dNdeta350->Integral();

 

 //Draw on Canvas---------------------------------------------------
   
 TCanvas* c_dNdeta = new TCanvas("c_dNdeta","dNdeta",700,700);
 //gStyle->SetOptStat(0);
 c_dNdeta->cd();

 h_dNdeta90->Draw();
 h_dNdeta160->Draw("same");
 h_dNdeta240->Draw("same");
 h_dNdeta350->Draw("same");
 
 c_dNdeta->SetGridx();
 c_dNdeta->SetGridy();

 //Building Legend-----------------------------------------------------

 TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
 //leg->SetHeader("Legend");
 leg->AddEntry(h_dNdeta90,Form("#sqrt{s}=90,  %.2f charged leptons",total_leptons[0]),"lep");
 leg->AddEntry(h_dNdeta160,Form("#sqrt{s}=160,  %.2f charged leptons",total_leptons[1]),"lep");
 leg->AddEntry(h_dNdeta240,Form("#sqrt{s}=240,  %.2f charged leptons",total_leptons[2]),"lep");
 leg->AddEntry(h_dNdeta350,Form("#sqrt{s}=350,  %.2f charged leptons",total_leptons[3]),"lep");
 leg->Draw();
cout <<"ciao" << endl;

//TLegend *leg2 =new TLegend(0.5945559,0.7688889,0.8968481,0.8859259,"#splitline{|#eta|<1.5}{5 < W < 35 GeV}");
//leg2->Draw();


 // the following line is needed to avoid that the automatic redrawing of stats
 h_dNdeta90->SetStats(0);
 h_dNdeta160->SetStats(0);
 h_dNdeta240->SetStats(0);
 h_dNdeta350->SetStats(0);

 c_dNdeta->Modified();


}
