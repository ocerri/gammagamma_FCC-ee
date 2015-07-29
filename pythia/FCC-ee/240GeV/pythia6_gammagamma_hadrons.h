
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

//bool debug = true;
bool debug = false;

using namespace std;

TDatabasePDG *PDG = 0;

TH1F* hdsigmadeta = 0;
TH1F* hdsigmadpT = 0;
TH1F* h_had_per_ev = 0;
TH1F* hW = 0;
TH1F* h_id_part = 0;

//TH1F* hEdsigmadpT = 0;
//TNtuple* ntFFdsigmadeta = 0;

//****creating a pdf for dsigmadpT fitting

Double_t pdf(Double_t *x, Double_t*par){
  Double_t _alpha = par[2];
  Double_t _C = par[1];
  Double_t _m = par[0];

  Double_t value = _C*TMath::Power(x[0],_alpha) + _m;
  return value;
}

//____________________________________________________________________
//
// Transforms a TMCparticle into a TParticle
//

TParticle *TMCParticle2TParticle( TMCParticle *part )
{

  if (!part) return 0;

  int pdg = -999, status = -999, mother1 = -999, mother2 = -999, daughter1 = -999, daughter2 = -999; 
  float px = -999., py = -999., pz = -999., etot = -999., vx=-999., vy=-999.,vz=-999., lifetime = -999.;
  
  pdg = part->GetKF(); status = part->GetKS();
  mother1 = part->GetParent();  mother2 = part->GetParent();
  daughter1 = part->GetFirstChild(); daughter2 = part->GetLastChild();
  px = part->GetPx(); py = part->GetPy(); pz = part->GetPz();
  etot = part->GetEnergy(); vx = part->GetVx(); vy = part->GetVy(); vz = part->GetVz(); 
  lifetime= part->GetTime();
  
  TParticle *partic = new TParticle(pdg, status, mother1,mother2, daughter1,daughter2, px,py,pz,etot, vx, vy, vz, lifetime);
  return partic;
  
}

//____________________________________________________________________
//
// Print out all the PYTHIA event-record for debugging purposes ...
//

void dumpEvent( TPythia6* pythia )
{

  TClonesArray* particleArray = (TClonesArray*)pythia->GetListOfParticles();
  int N = particleArray->GetEntriesFast();
  cout << "<I> PYTHIA particles in the event: " << N << " " << endl;

  for (int ip=0; ip<N; ip++){

    TMCParticle *part = (TMCParticle*)particleArray->At(ip);
    int status = part->GetKS();
    //if (status==1) continue; // Exclude final particles
    int pdg = part->GetKF();
    //if ( abs(pdg)>5 && abs(pdg)!=21) continue; // Exclude everything but quarks (udscb) and gluons
    double ptPart  = sqrt( part->GetPx()*part->GetPx()+ part->GetPy()*part->GetPy() );
    if (debug) cout << "<I> Particle at position " << ip << " with PDG=" << pdg // << ", Energy=" << part->GetEnergy() 
		    << ", Status=" << status << " and pT=" << ptPart << endl;
  }
}

//____________________________________________________________________
//
// Apply MB trigger cuts
//

bool passEvtSelection( TPythia6* pythia )
{

  //return true;  // ACCEPT ALL EVENTS

  TClonesArray* particleArray = (TClonesArray*)pythia->GetListOfParticles();
  int N = particleArray->GetEntriesFast();

  // loop on event particles to check if it satisfies trigger cuts
  Double_t EN_tot=0;
  Double_t pz_tot=0;
  Double_t px_tot=0;
  Double_t py_tot=0;

  for (int ip=0; ip<N; ip++){

    TMCParticle *part = (TMCParticle*)particleArray->At(ip);

    int status = part->GetKS();
    if (status!=1) continue; // Only consider final-state particles

    int pdg = part->GetKF();
    double charge = PDG->GetParticle(pdg)->Charge();

    TParticle *partic = (TParticle*)TMCParticle2TParticle(part);

    double etaPart = partic->Eta();
    double enPart = partic->Energy();

    if ( abs(pdg)<11 || abs(pdg)>16 ){
      EN_tot+= enPart;
      pz_tot += partic->Pz();
      px_tot += partic->Px();
      py_tot += partic->Py();
    }

  }

  Double_t W_gen = TMath::Sqrt(EN_tot*EN_tot - px_tot*px_tot - py_tot*py_tot - pz_tot*pz_tot );

  hW->Fill(W_gen); //Fill the W histo

  //  if (W_gen<35 && W_gen>5) return true; //Cut on W !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  return true;

  //return false;

}

//____________________________________________________________________
//
// Book histos
// 

Double_t pt_max = 16.;
Double_t pt_min = 0.1;
Int_t n_bin_pt = 40;

void book_histos()
{

  char title[300];  
  sprintf(title, "hdsigmadeta");
  hdsigmadeta = new TH1F(title,title,20,-10.5,10.5);
  hdsigmadeta->SetXTitle("|#eta|");
  hdsigmadeta->SetYTitle("d#sigma_{ch}/d|#eta| (nb)");
  hdsigmadeta->SetMinimum(0.0001);
  hdsigmadeta->Sumw2();

  sprintf(title, "h_hadrons_per_event");
  h_had_per_ev = new TH1F(title,title,20,0,20);
  h_had_per_ev->SetXTitle("Number of charged hadrons per event");
  h_had_per_ev->SetYTitle("Number of events");
  h_had_per_ev->SetMinimum(0.000001);
  h_had_per_ev->Sumw2();

  sprintf(title, "h_id_particle");
  h_id_part = new TH1F(title,title,5000,0,5000);
  h_id_part->SetXTitle("Particle id");
  h_id_part->SetYTitle("Number particle detected");
  h_id_part->SetMarkerSize(1.4);
  h_id_part->SetMarkerStyle(33);
  h_id_part->Sumw2();
 
  sprintf(title, "hdsigmadpT");
  hdsigmadpT = new TH1F(title,title,n_bin_pt, pt_min , pt_max);
  hdsigmadpT->SetXTitle("p_{T} (GeV/c)");
  hdsigmadpT->SetYTitle("d#sigma/dp_{T} (pb/GeV)");
  hdsigmadpT->SetMinimum(0.001);
  hdsigmadpT->Sumw2();

  sprintf(title, "hW");
  hW = new TH1F(title,title,40,0,196);
  hW->SetXTitle("W(GeV)");
  hW->SetYTitle("Nuymber of events");
  hW->SetMinimum(0.0001);
  hW->Sumw2();

  return;

}

