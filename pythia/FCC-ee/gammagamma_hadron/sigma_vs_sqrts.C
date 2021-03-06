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

TGraphErrors* g_sigma_tot = 0;
TGraphErrors* g_sigma_c = 0;

void sigma_vs_sqrts(){

Double_t energy[]={90,160,240,350};//Gev
Double_t de[]={0,0,0,0};
Double_t sigma_tot[]={102.2, 183.0, 266.2, 363.9};//nb
Double_t dsigma_tot[]={0.003*102.2, 0.003*183.0, 0.003*266.2, 0.003*363.9};//nb
Double_t sigma_c[]={51.4,90.9,131.6,181.4}; //nb
Double_t dsigma_c[]={0.003*51.4,0.003*90.9,0.003*131.6,0.003*181.4}; //nb

g_sigma_tot = new TGraphErrors(4,energy,sigma_tot,de,dsigma_tot);
g_sigma_tot->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
g_sigma_tot->GetYaxis()->SetTitle("#sigma_{tot}");

TAxis *axis = g_sigma_tot->GetYaxis();
axis->SetLimits(0.,400.);    

g_sigma_c = new TGraphErrors(4,energy,sigma_c,de,dsigma_c);
g_sigma_c->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
g_sigma_c->GetYaxis()->SetTitle("#sigma_{tot}");
//g_sigma_c->SetMarkerStyle(24);
g_sigma_c->SetMarkerColor(2);
g_sigma_c->SetLineColor(2);


TCanvas *c_sigma = new TCanvas("c_sigma","c_sigma",700,700);
c_sigma->cd();
g_sigma_tot->Draw("AP");
g_sigma_c->Draw("P");
c_sigma->SetGridx();
c_sigma->SetGridy();

}
