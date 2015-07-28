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
 h_dsdpt10->SetAxisRange(0.01,100000,"Y");
 h_dsdpt10->SetTitle("PYTHIA simulation of d#sigma/dp_{t} for e^{+}e^{-} #rightarrow h^{#pm} at 195.5 GeV");
 
 
 TFile* file30_dsdpt = new TFile("_root/MSTP14-30_final_out.root");
 h_dsdpt30 = (TH1F*)file30_dsdpt->Get("hdsigmadpT");
 h_dsdpt30->SetLineColor(4);
 h_dsdpt30->SetMarkerColor(4);
 h_dsdpt10->SetAxisRange(0.01,100000,"Y");
 h_dsdpt30->SetTitle("PYTHIA simulation of d#sigma/dp_{t} for {e}^{+} {e}^{-} -> {h}^{#pm} at 195.5 GeV");

 // TFile* file10_dsdeta = new TFile("_root/MSTP14-10_final.root");
 // h_dsdeta10 = (TH1F*)file10_dsdeta->Get("hdsigmadeta");
 // h_dsdeta10->SetMarkerColor(2);
 // h_dsdeta10->SetLineColor(2);
 // h_dsdeta10->SetAxisRange(8.,24.,"Y");
 // h_dsdeta10->SetTitle("PYTHIA simulation of d#sigma/deta for e^{+}e^{-} #rightarrow h^{#pm} at 195.5 GeV");
 
 
 // TFile* file30_dsdeta = new TFile("_root/MSTP14-30_final.root");
 // h_dsdeta30 = (TH1F*)file30_dsdeta->Get("hdsigmadeta");
 // h_dsdeta30->SetMarkerColor(4);
 // //h_dsdeta30->SetLineColor(4);
 // h_dsdeta30->SetAxisRange(8.,24.,"Y");
 // h_dsdeta30->SetTitle("PYTHIA simulation of d#sigma/deta for {e}^{+} {e}^{-} -> {h}^{#pm} at 195.5 GeV");

 //Plotting experimental data from txt-------------------------------------------

 // std::ifstream myfile;
 // myfile.open("_txt/dsdpt_opal.txt");

 // const Int_t n_data = 21; //number of data (rows) in .txt

 // Double_t x[n_data];
 // Double_t y[n_data];
 // Double_t ex[n_data];
 // Double_t ey[n_data];

 // int j=0;

 // //Read data from file .txt
 // while(true) {
 //   myfile >> x[j];
 //   myfile >> ex[j];
 //   myfile >> y[j];
 //   myfile >> ey[j];
 //   if( myfile.eof() ) break;
 //   ++j;
 // }

 // myfile.close();

 // TGraphErrors* gr = new TGraphErrors(n_data,x,y,ex,ey); //Create TGraph
 
 // gr->SetName("gr");

 //Exp data from internet-----------------------------------------------------------
  double p7661_d1x1y2_xval[] = { 1.76, 2.17, 2.58, 2.98, 3.38, 3.78, 4.18, 4.59, 4.99, 
    5.39, 5.79, 6.19, 6.59, 6.98, 7.38, 7.78, 8.44, 9.47, 10.87, 
    13.53 };
  double p7661_d1x1y2_xerrminus[] = { 0.15999999999999992, 0.16999999999999993, 0.18000000000000016, 0.18000000000000016, 0.17999999999999972, 0.17999999999999972, 0.17999999999999972, 0.1899999999999995, 0.1900000000000004, 
    0.1899999999999995, 0.1900000000000004, 0.1900000000000004, 0.1899999999999995, 0.1800000000000006, 0.17999999999999972, 0.1800000000000006, 0.4399999999999995, 0.47000000000000064, 0.8699999999999992, 
    1.5299999999999994 };
  double p7661_d1x1y2_xerrplus[] = { 0.24, 0.22999999999999998, 0.21999999999999975, 0.2200000000000002, 0.2200000000000002, 0.2200000000000002, 0.22000000000000064, 0.20999999999999996, 0.20999999999999996, 
    0.20999999999999996, 0.20999999999999996, 0.20999999999999996, 0.20999999999999996, 0.21999999999999975, 0.21999999999999975, 0.21999999999999975, 0.5600000000000005, 0.5299999999999994, 1.1300000000000008, 
    2.4700000000000006 };
  double p7661_d1x1y2_yval[] = { 300.0, 115.0, 52.3, 26.6, 16.1, 9.41, 5.54, 3.89, 2.78, 
    1.65, 1.7, 1.16, 0.834, 0.665, 0.543, 0.367, 0.265, 0.171, 0.068, 
    0.023 };
  double p7661_d1x1y2_yerrminus[] = { 60.03332407921454, 17.029386365926403, 8.238931969618392, 3.8470768123342687, 2.5495097567963922, 1.725862103413827, 0.6037383539249432, 0.5277309920783504, 0.223606797749979, 
    0.17088007490635063, 0.30610455730027936, 0.2051828452868319, 0.2987507322166759, 0.31087618113969423, 0.11140017953306897, 0.16919810873647495, 0.260823695242591, 0.1558717421471897, 0.05643580423808985, 
    0.01565247584249853 };
  double p7661_d1x1y2_yerrplus[] = { 42.04759208325728, 9.055385138137417, 2.816025568065745, 3.059411708155671, 0.7071067811865476, 1.087841900277793, 0.8918520056601319, 0.48373546489791297, 0.35227829907617075, 
    0.24839484696748443, 0.22472205054244232, 0.21213203435596426, 0.1438923208513922, 0.2757680184502909, 0.11253888216967503, 0.09953893710503443, 0.22970633426181353, 0.1651575005865613, 0.13983204210766573, 
    0.043566041821583926 };
  double p7661_d1x1y2_ystatminus[] = { 2.0, 1.0, 0.8, 0.6, 0.5, 0.35, 0.27, 0.24, 0.2, 
    0.16, 0.19, 0.15, 0.136, 0.112, 0.109, 0.092, 0.123, 0.086, 0.028, 
    0.007 };
  double p7661_d1x1y2_ystatplus[] = { 2.0, 1.0, 0.8, 0.6, 0.5, 0.35, 0.27, 0.24, 0.2, 
    0.16, 0.19, 0.15, 0.136, 0.112, 0.109, 0.092, 0.123, 0.086, 0.028, 
    0.007 };
  int p7661_d1x1y2_numpoints = 20;
  TGraphAsymmErrors *p7661_d1x1y2 = new TGraphAsymmErrors(p7661_d1x1y2_numpoints, p7661_d1x1y2_xval, p7661_d1x1y2_yval, p7661_d1x1y2_xerrminus, p7661_d1x1y2_xerrplus, p7661_d1x1y2_yerrminus, p7661_d1x1y2_yerrplus);
  p7661_d1x1y2->SetName("Delphi data");
  p7661_d1x1y2->SetTitle("Delphi data");
  p7661_d1x1y2->SetLineWidth(2);
  p7661_d1x1y2->SetMarkerStyle(5);
  p7661_d1x1y2->SetMarkerSize(1);
  p7661_d1x1y2->SetMarkerColor(1);

 

 //Draw on Canvas---------------------------------------------------
   
 TCanvas* c_dsdpt = new TCanvas("c_dsdpt","dsigmadipt",700,700);
 //gStyle->SetOptStat(0);
 c_dsdpt->cd();

 h_dsdpt10->Draw();
 h_dsdpt30->Draw("same");
 p7661_d1x1y2->Draw("AP");
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

 //Building Legend-----------------------------------------------------

 TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
 //leg->SetHeader("Legend");
 leg->AddEntry(h_dsdpt10,"Histogram MSTP(14)=10","lep");
 leg->AddEntry(h_dsdpt30,"Histogram MSTP(14)=30","lep");
 leg->AddEntry("p7661_d1x1y2","Experimental data from DELPHI","lep");
 leg->Draw();

 //Adding comments in the statistics Box-------------------------------------------


 // TPaveText *pt = new TPaveText(.05,.1,.95,.8);
 // pt->AddText("A TPaveText can contain severals line of text.");
 // pt->AddText("They are added to the pave using the AddText method.");
 // pt->AddLine(.0,.5,1.,.5);
 // pt->AddText("Even complex TLatex formulas can be added:");
 // pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
 // pt->Draw();

 TPaveStats *ps = (TPaveStats*)c_dsdpt->GetPrimitive("stats");
 ps->SetName("mystats");
 TList *list = ps->GetListOfLines();

 // Remove the statistics lines
 TText *tconst = ps->GetLineWith("RMS");
 list->Remove(tconst);
 tconst = ps->GetLineWith("hdsigmadpT");
 list->Remove(tconst);
 tconst = ps->GetLineWith("Entries");
 list->Remove(tconst); 
 tconst = ps->GetLineWith("Mean");
 list->Remove(tconst); 

 // Add a new line in the stat box.
 TLatex *myt = new TLatex(0,0,"#cbar#eta#cbar< 1.5");
 myt ->SetTextFont(42);
 myt ->SetTextSize(0.01);
 list->Add(myt);
 myt = new TLatex(0,0,"5<W<35 GeV");
 myt ->SetTextFont(42);
 myt ->SetTextSize(0.01);
 list->Add(myt);

 // the following line is needed to avoid that the automatic redrawing of stats
 h_dsdpt10->SetStats(0);
 h_dsdpt30->SetStats(0);

 c_dsdpt->Modified();


}
