//____________________________________________________________________
/*

 NOTE 1: To run this code, you must have a version of ROOT compiled with Pythia6 enabled 
         and installed in $HOME/progs/pythia/pythia6/libPythia6

 NOTE 2: Compile/Run the code via ACLIC:


  gSystem->Load("libEG"); 
  gSystem->Load("$HOME/progs/pythia/pythia6/libPythia6"); 
  gSystem->Load("libEGPythia6");  

  .L /home/enterria/wrk/gammagamma_hadrons_fcc_ee/pythia6_gammagamma_hadrons.C++

  //Or once compiled: gSystem->Load("pythia6_gammagamma_hadrons_C.so");  

  pythia6_gammagamma_hadrons(10000, 160. );

*/
//____________________________________________________________________

#include "pythia6_gammagamma_hadrons.h"

//____________________________________________________________________
//

void pythia6_gammagamma_hadrons( int Nevts = 5000, double sqrts = 160.) 
{

  // Instance of the Pythia event generator

  TPythia6* pythia = new TPythia6();

  PDG = TDatabasePDG::Instance(); // TDataBasePDG contains info on all particle properties

  // Random seeding
  int seed = (int)(time(NULL)/3);
  if( (seed>=0) && (seed<=900000000) ) {
    pythia->SetMRPY(1,seed); // set seed
    pythia->SetMRPY(2,0); // use new seed
    cout << endl << "<I> Random Seed : " << seed << endl; 
  } 
  else{
    cout << endl << "<I> Default random seed ... "<< seed << endl; 
  } 

  //____________________________________________________________________
  //
  // PYTHIA GENERATION SETTINGS
  //____________________________________________________________________
  

  // *******************************************************************
  // KINEMATICS: acceptance CUTS. 

  //   double etawidth = 10.;
  //   double etamin = -etawidth;
  //   double etamax = +etawidth;
  
  //   // range of allowed pseudorapidities for the particle with largest pseudorapidity in a 2 --> 2 process
  //   pythia->SetCKIN(13,etamin);     //! etamin',
  //   pythia->SetCKIN(14,etamax);    //! etamax',
  //   // range of allowed pseudorapidities for the particle with smallest pseudorapidity in a 2 --> 2 process
  //   pythia->SetCKIN(15,-etamax);   //! -etamax',
  //   pythia->SetCKIN(16,-etamin);    //! -etamin' 
  
  //   cout << endl << "<I> Selected kinematics: eta = [ " << etamin << ", " << etamax << endl;
  
  // *******************************************************************
  // PHYSICS PROCESS SELECTION
  
  pythia->SetMSEL(0); // we choose the process by hand below 
  //pythia->SetMSTP(61, 1); // switch on/off ISR
  //pythia->SetMSTP(71, 1); // switch on/off FSR
  pythia->SetMSTP(81, 1); // Multiple parton interactions for resolved gamma-gamma colls.


  // *******************************************************************
  // Final-state:

  //pythia->Pyedit(2); // Only keep final-state particles in data-record (neutrinos removed too, compared to Pyedit(1))

  pythia->SetMSTJ(22,2);      //! Decay those unstable particles',
  pythia->SetPARJ(71,10.);    //! for which ctau < 10 mm',

  //pythia->SetMSTP(41,0); // decays turned off for all resonances

  // We switch off all relevant meson decays into dsigmadeta
  // KC = 1 -- 500 compressed code for particles 
  //   int pi0_KC=pythia->Pycomp(111); // get KC code from PDG-code
  //   pythia->SetMDCY(pi0_KC,1,0);  // No decay of pi0

  //  All decays switched-off
  //for(int i=-3000;i<3000;i++) pythia->SetMDCY(i,1,0);  // No decay of anything


  // *******************************************************************
  // Detailed process selection 

  /* See PY6 example: https://pythia6.hepforge.org/examples/main69.f
     Set minimal and maximal W of gamma*gamma* system.
     CKIN(77)=50D0
     CKIN(78)=200D0
     Set minimum photon virtuality.
     CKIN(65)=0D0
     CKIN(67)=0D0    
     Set maximum photon virtuality.
     CKIN(66)=4D0
     CKIN(68)=4D0  

     //Modify contribution of resolved longitudinal photons.
     //MSTP(17)=4
     //PARP(165)=1D0 

     //Simplify simulation.
     //MSTP(61)=0 
     //MSTP(71)=0 
     //MSTP(81)=0 
     //MSTP(91)=0 
     //MSTP(111)=0 

     //Use minimal W to set minimum photon energy fractions. 
     CKIN(61)=(CKIN(77)/ECM)**2  
     CKIN(63)=CKIN(61)
     CKIN(73)=CKIN(61)
     CKIN(75)=CKIN(61)
     //Check that not too low pTmin.
     IF(MOPT.EQ.1) THEN
     PTMINS=PARP(81)*(ECM/PARP(89))**PARP(90) 
     CKIN(3)=MAX(CKIN(3),1.01D0*PTMINS)
     ENDIF

     //
     CALL PYINIT('CMS','GAMMA/E+','GAMMA/E-',ECM)

     // Book histograms
     CALL PYBOOK(1,'x = energy fraction of photons in lepton',&100,0D0,1D0)
     CALL PYBOOK(2,'W distribution',100,0D0,200D0)
     CALL PYBOOK(3,'log10(Q2)',100,-15D0,10D0)
     CALL PYBOOK(4,'x = energy fraction of partons in photon',&100,0D0,1D0)
     CALL PYBOOK(5,'pT of hard scattering',100,0D0,20D0)
     CALL PYBOOK(6,'number of hard interactions in event',&20,-0.5D0,19.5D0)

     // LOOP OVER EVTS
     ISUB=MSTI(1)
     X1G=PARI(103)
     X2G=PARI(104)
     W=PARI(11)
     Q21=PARI(105)
     Q22=PARI(106)
     // Shift a direct photon away from x_p = 1 for better histogramming.
     X1P=MIN(0.999D0,PARI(33))
     X2P=MIN(0.999D0,PARI(34))
     // pT=0 for low-pT, but nonvanishing for elastic/diffractive events.
     PT=PARI(17)
     NMUL=MSTI(31)
     
     //Fill histograms. End event loop.
     CALL PYFILL(1,X1G,1D0)
     IF(ISEL.EQ.2) CALL PYFILL(1,X2G,1D0)
     CALL PYFILL(2,W,1D0)
     CALL PYFILL(3,LOG10(Q21),1D0)
     IF(ISEL.EQ.2) CALL PYFILL(3,LOG10(Q22),1D0)
     CALL PYFILL(4,X1P,1D0)
     IF(ISEL.EQ.2) CALL PYFILL(4,X2P,1D0)
     CALL PYFILL(5,PT,1D0)
     CALL PYFILL(6,DBLE(NMUL),1D0)

  */

  //pythia->SetMSTP(14,10);  // gamma gamma -> hadrons
  pythia->SetMSTP(14,30);  // gamma gamma -> hadrons (FULL)
  pythia->SetMSEL(2); // min-bias QCD

//   cout << "****************** DEFAULT PYTHIA minbias QCD settings: *****************" << endl;
//   cout << " *     MSTP(81) =    " << pythia->GetMSTP(81) << "UE model                                     *" << endl;
//   cout << " *     MSTP(82) =    " << pythia->GetMSTP(82) << "MPI model                                     *" << endl;
//   cout << " *     PARP(82) =    " << pythia->GetPARP(82) << "   UE IR cutoff at reference ecm                *" << endl;
//   cout << " *     PARP(89) =    " << pythia->GetPARP(89) << "   UE IR cutoff reference ecm                   *" << endl;
//   cout << " *     PARP(90) =    " << pythia->GetPARP(90) << "   UE IR cutoff ecm scaling power               *" << endl; 
//   cout << " *     MSTP(51) =    " << pythia->GetMSTP(51) << "   PDF set                                 *" << endl;
//   cout << " *     MSTP(52) =    " << pythia->GetMSTP(52) << "   PDF pdflib(2/internal(1)                   *" << endl;


  // ******************************************************************


  book_histos();

  // ******************************************************************
  //  Initialise it to run e+/e- --> gamma gamma --> X

  pythia->Initialize("cms", "gamma/e+", "gamma/e-", sqrts);
  //pythia->Initialize("cms", "e+", "e-", sqrts);
  cout << "SYSTEM: e+e- at sqrt(s)=" << sqrts << " GeV" << endl; 

  TClonesArray* particlesArray = 0;

  // ******************************************************************
  // EVENT LOOP
  
  int exclEvts = 0;

  for (int ievt = 0; ievt < Nevts; ievt++) {

    if (ievt ==0) cout << "Generating event #: " << flush;
    if (ievt % 100 == 0) cout << ievt << " " << flush;


    // ******************************************************************
    // Generate event

    pythia->GenerateEvent();
    
    // ***********************************************************************************************
    // Check if event passes trigger selection
    // 
    //  if (debug) cout << "############ EVENT " << ievt << " #############" << endl;
    //	dumpEvent( pythia );

    if ( passEvtSelection(pythia)==false ) { ///Cut on W!!!!!!!!!!!!!!!
      exclEvts++;
      continue;
    }

    // Loop over particles and fill histos

    particlesArray = (TClonesArray*)pythia->GetListOfParticles();
    int N = particlesArray->GetEntriesFast();
    //cout << "<I> PYTHIA particles in event: " << N << " " << flush;

    //int Nch = 0;
    TMCParticle *part;
    TParticle *partic;

    for (int ip = 0; ip < N; ip++ ) {

      part = (TMCParticle*)particlesArray->At(ip);
      int status = part->GetKS();
      if (status!=1) continue; // Only consider final-state particles

      int pdg = part->GetKF();
      double charge = PDG->GetParticle(pdg)->Charge();
      if(charge==0) continue; // only charged particles

      //if( abs(pdg)!=211 && abs(pdg)!=321 && abs(pdg)!=2212 && abs(pdg)!=3122 && // charged pions, kaons, protons, lambdas
      //	  abs(pdg)!=11 && abs(pdg)!=13 && abs(pdg)!=15 ) continue; // leptons

      if ( abs(pdg)==11 && abs(pdg)==13 && abs(pdg)==15) continue;// no leptons
      if (abs(pdg)==24) continue; //no W ....even if they should't
      //fly for 10 mm
	

      partic = (TParticle*)TMCParticle2TParticle(part);
      double ptPartic = partic->Pt();
      double etaPartic = partic->Eta();
      //phiPartic = partic->Phi();

      if ( ptPartic<0.12 ) continue; //cut on Pt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (abs(etaPartic)>1.5) continue; //Fill only if |eta|<1.5!!!!!!!!!!!!!!!!!!!!!   
      // Histo filling

      hdsigmadeta->Fill(abs(etaPartic));
      
      // if (TMath::Abs(etaPartic)<etaRange) {
      hdsigmadpT->Fill(ptPartic);
      //hEdsigmadpT->Fill(ptPartic,1/(2*TMath::Pi()*ptPartic));
      //}
      delete partic;
      // if (TMath::Abs(etaPartic)<1.) Nch++;
    } // End loop over all particles in event



    
    //hmulteta1->Fill(Nch);

    //     // Checkpointing histo every 10000 evts.
    //     if (ievt % 10000 == 0) {
    //       //hpTdsigmadeta->Write("",TObject::kOverwrite);
    //       hdsigmadeta->Write("",TObject::kOverwrite);
    //       //ntFFdsigmadeta->Write("",TObject::kOverwrite);
    //     }
  } // END EVENT GENERATION LOOP

  // **********************************************************************************  
  // FINAL GENERATION MESSAGES:
  
  pythia->Pystat(1);
  
  double xsection = pythia->GetPARI(1) * 1e-9;//conversion form mb to pb
  int    ntrials  = pythia->GetMSTI(5);
  double sigmaweight = xsection/ntrials;

  cout << endl << "#######################" << endl 
       << endl << "<I> Event summary: " << Nevts << " generated." << endl ;
  
  //double evts = hdsigmadeta->Integral(); 
  double triggEff = (Nevts-exclEvts)/(float)Nevts;
  cout << "<I> Num. of total valid events = " << Nevts-exclEvts << " (Effic="<< triggEff*100. << "%)" << endl;

  cout << "<I> Pythia cross section: " << xsection << " mb || number of trials: " << ntrials << " || sigmaweight = "<< sigmaweight <<endl;

  // double dNchdeta = hdsigmadeta->GetBinContent(11)/ntrials; 
  // cout << "<I> dN_ch/deta|_{eta=0} = " << dNchdeta << " in e+e- --> gamma gamma --> X at sqrt(s) = " << sqrts << " GeV " << endl ;

  // **********************************************************************************  
  // Normalize histos by weighted cross-section, pseudorapidity range, and the pT bin size
       
  double etabinsize = 1.5/5; // eta binning: 20 within -10<eta<10
  hdsigmadeta->Scale(sigmaweight/etabinsize);
  double ptbinsize = (pt_max-pt_min)/n_bin_pt;
  hdsigmadpT->Scale(sigmaweight/ptbinsize);
    
  //hdsigmadetaTruth->Scale(sigmaweight*0.01); // eta binning: 2000 in -10<eta<10
  //hEdsigmadpT->Scale(ptbinsize/ntrials);
  //ntFFdsigmadeta->SetWeight(sigmaweight);


  // **********************************************************************************  
  // Plot distributions for cross-check


//   //hEdsigmadpT->Rebin(4);
  char title[300];
    sprintf(title, "cinvdsigmadpT_%iGeV",(int)sqrts);
    TCanvas *cinvdsigmadpT = new TCanvas(title,title,700,600);
    cinvdsigmadpT->SetLogy();
    cinvdsigmadpT->SetLogx();
    cinvdsigmadpT->cd();
    hdsigmadpT->Draw();
    cinvdsigmadpT->SaveAs("hdsigmadpT.png");    


    sprintf(title, "cinvdsigmadeta_%iGeV",(int)sqrts);
    TCanvas *cinvdsigmadeta = new TCanvas(title,title,700,600);
    cinvdsigmadeta->cd();
    hdsigmadeta->Draw();
    cinvdsigmadeta->SaveAs("hdsigmadeta.png");


    sprintf(title, "cinvW_%iGeV",(int)sqrts);
    TCanvas *cinvW = new TCanvas(title,title,700,600);
    cinvW->cd();
    hW->Draw();
    cinvW->SaveAs("hW.png");

    
  //*********************FITTING hdsigmadpT using a power law*******************

  TF1 *f_fit = new TF1("f_fit", pdf, pt_min, pt_max, 3);
  f_fit->SetParNames("m","C","alpha");
  
  cout << "Starting fitting to dsigmadpT..." << endl << "Give starting parameters for pdf(x) =C*x^alpha +m" << endl;

  Double_t m = 5e4;
  cout << "m :"; cin >> m;

  Double_t C = 1e4;
  cout << "C: "; cin >> C;

  Double_t alpha = -8;
  cout << "alpha: "; cin >> alpha;
  
  f_fit->SetParameter(0, m);
  f_fit->SetParameter(1, C);
  f_fit->SetParameter(2, alpha);

  hdsigmadpT->Fit(f_fit,"","",pt_min,pt_max);

  Double_t fit_chi2 = f_fit->GetChisquare();
  Int_t    fit_ndof = f_fit->GetNDF();
  Double_t fit_prob = f_fit->GetProb();
  cout << endl;
  cout << "chi2 = " << fit_chi2 << endl;
  cout << "ndof = " << fit_ndof << endl;
  cout << "prob = " << fit_prob << endl;
  cout << endl;

  // **********************************************************************************  
  // Open  output file and Close file
    
  char filename[200];
  sprintf(filename, "pythia6_gammagamma_hadrons_%iGeV_seed%d_Nevts%d.root",(int)TMath::Ceil(sqrts),seed,Nevts);

  TFile* file = TFile::Open(filename, "RECREATE");
  if (!file || !file->IsOpen()) {
    Error("pythia6_gammagamma_hadrons", "Couldn;t open file %s", filename);
    return;
  }

  file->cd();
  hdsigmadeta->Write();
  hdsigmadpT->Write();
  hW->Write();
  //file->Write("",TObject::kOverwrite);
  file->Close();
  cout << endl << "#######<I> File " << filename << " created. Take a look ... ##############" << endl << endl;

  file = TFile::Open(filename);
  file->ls();
  file->Close();
}


//____________________________________________________________________
//
