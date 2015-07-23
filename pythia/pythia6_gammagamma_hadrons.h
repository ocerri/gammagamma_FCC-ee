
#include <Riostream.h>
#include <cstdlib>

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
//TH1F* hEdsigmadpT = 0;
//TNtuple* ntFFdsigmadeta = 0;

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

  

  bool cntTrack = false;

  bool HFplus = false;
  bool HFminus = false;

  bool V0plus = false;
  bool V0minus = false;

  bool fwdTrackplus1 = false;
  bool fwdTrackminus1 = false;
  bool fwdTrackplus5 = false;
  bool fwdTrackminus5 = false;

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

    if ( charge!=0 && TMath::Abs(etaPart)<2. ) cntTrack = true;


    if ( charge!=0 && etaPart>2.  && etaPart<5.6 ) fwdTrackplus5 = true;
    if ( charge!=0 && etaPart<-2. && etaPart>-5.6 ) fwdTrackminus5 = true;

    if ( charge!=0 && etaPart>1.5  && etaPart<5. ) fwdTrackplus1 = true;
    if ( charge!=0 && etaPart<-1.5 && etaPart>-5. ) fwdTrackminus1 = true;

    if ( etaPart> 2.9 && etaPart< 5.2 && enPart > 3. ) HFplus=true;
    if ( etaPart<-2.9 && etaPart>-5.2 && enPart > 3. ) HFminus=true;

    if ( etaPart>2.8 && etaPart<5.1 ) V0plus=true;
    if ( etaPart<-1.7 && etaPart>-3.7 ) V0minus=true;

  }

  Double_t W_gen = TMath::Sqrt(EN_tot*EN_tot - px_tot*px_tot - py_tot*py_tot - pz_tot*pz_tot );


  if (W_gen<125 && W_gen>10) return true; //Cut on W !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
//   // The minbias trigger was a two-arm trigger that required at least two
//   // charged particles in opposite rapidity hemispheres in the range 1.5 < eta <5.5.
//   if ( trigger.Contains("ua1") && (fwdTrackplus1==true && fwdTrackminus1==true ) ) return true;

//   if ( trigger.Contains("ua5") && (fwdTrackplus5==true && fwdTrackminus5==true) ) return true;

//   // One track within |eta|<2
//   if ( trigger.Contains("aliceMB0") && (cntTrack==true) ) return true;
//   // One particle in 2.8 < eta < 5.1 and one in -1.7 < eta < -3.7
//   else if ( trigger.Contains("aliceMB") && (V0plus==true && V0minus==true)  ) return true;

//   if ( trigger.Contains("cms") && (HFplus==true && HFminus==true) ) return true;

  return false;

}

//____________________________________________________________________
//
// Book histos
// 

void book_histos()
{

  char title[300];  
  sprintf(title, "hdsigmadeta");
  //hdsigmadeta = new TH1F(title,title,21,-10.5,10.5);
  hdsigmadeta = new TH1F(title,title,50,0,10);
  hdsigmadeta->SetXTitle("|#eta|");
  hdsigmadeta->SetYTitle("d#sigma_{ch}/d|#eta| (mb)");
  hdsigmadeta->SetMinimum(0.);
  hdsigmadeta->Sumw2();

  //   hdsigmadetaTruth = new TH1F("hdsigmadetaTruth","hdsigmadetaTruth",2000,-10,10);
  //   hdsigmadetaTruth->SetXTitle("#eta_{MC truth}");
  //   hdsigmadetaTruth->SetYTitle("d#sigma/d#eta [mb]");
  //   hdsigmadetaTruth->SetMinimum(1.e-5);
  //   hdsigmadetaTruth->Sumw2();
 
  sprintf(title, "hdsigmadpT");
  hdsigmadpT = new TH1F(title,title,1000,0.,10.);
  hdsigmadpT->SetXTitle("p_{T} (GeV/c)");
  hdsigmadpT->SetYTitle("d#sigma/dp_{T} (mb/GeV)");
  hdsigmadpT->SetMinimum(0.0001);
  hdsigmadpT->Sumw2();

//   sprintf(title, "hEdsigmadpT");
//   hEdsigmadpT = new TH1F(title,title,1000,0.,10.);
//   hEdsigmadpT->SetXTitle("p_{T} (GeV/c)");
//   hEdsigmadpT->SetYTitle("1/(2#pi p_{T})dsigma/dp_{T} (GeV/c)^{-2}");
//   hEdsigmadpT->SetMinimum(0.0001);
//   hEdsigmadpT->Sumw2();

//   sprintf(title, "hmulteta1");
//   hmulteta1 = new TH1F(title,title,1000,-0.5,999.5);
//   hmulteta1->SetXTitle("N_{ch}");
//   hmulteta1->SetYTitle("dsigma/dsigma_{ch}");
//   hmulteta1->SetMinimum(0.00001);
//   hmulteta1->Sumw2();

  // Ntuple to store the FF=FF(z,ptPartic,pdgParton)
  //sprintf(title, "ntFFdsigmadeta_z_ptpartic_pdgparton_subproc_%s");
  //ntFFdsigmadeta = new TNtuple(title,title,"zT:z:pTpartic:Epartic:pdgparton");

  return;

}

