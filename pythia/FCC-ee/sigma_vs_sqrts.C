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
Double_t sigma_c[]={51.4,90.9,131.6,181.4}; //nb

g_sigma_c = new TGraphErrors(4,energy,sigma_c,de,dsigma_c);
g_sigma_c->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
g_sigma_c->GetYaxis()->SetTitle("#sigma");

g_sigma_tot = new TGraphErrors(4,energy,sigma_tot,de,dsigma_tot);
g_sigma_tot->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
g_sigma_tot->GetYaxis()->SetTitle("#sigma");

TCanvas *c_sigma

}
