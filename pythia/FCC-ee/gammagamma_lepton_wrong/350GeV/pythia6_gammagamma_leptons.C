#include "pythia6_gammagamma_leptons.h"

//____________________________________________________________________
//

void pythia6_gammagamma_leptons( int Nevts = 10000, double sqrts = 350, int MSTP14_val=10) 
{

  //Luminosity and bunches definitions according to FCC-ee project specifications
  Double_t Luminosity=1;
  Double_t Luminosity_c=1; //with opt crab waist
  Double_t FCC_Circumference = 100000; // meters
  Double_t speed_of_light = 299792458; // m/s
  Int_t N_bunch = 1; //number of bunches per beam
  Int_t N_bunch_c= 1;//with opt Crab waist

  if (sqrts==90){
    N_bunch = 16700;
    Luminosity = 0.28; //pb^-1 s^-1
    N_bunch_c = 59581;
    Luminosity_c = 2.15;
  }
  else if (sqrts==160){
    N_bunch = 4490;
    Luminosity = 0.12; //pb^-1 s^-1
    N_bunch_c =3143;
    Luminosity_c =0.38;
  }
  else if (sqrts==240){
    N_bunch = 1360;
    Luminosity =  0.06; //pb^-1 s^-1
    N_bunch_c = 625;
    Luminosity_c = 0.087;
  }
  else if (sqrts==350){
    N_bunch = 98;
    Luminosity = 0.0018; //pb^-1 s^-1
    N_bunch_c = 68;
    Luminosity_c = 0.021;
  }
  
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
  // PHYSICS PROCESS SELECTION
  
  pythia->SetMSEL(0); // we choose the process by hand below 
  pythia->SetMSTP(81, 1); // Multiple parton interactions for resolved gamma-gamma colls.


  // *******************************************************************
  // Final-state:

  pythia->SetMSTJ(22,2);      //! Decay those unstable particles',
  pythia->SetPARJ(71,10.);    //! for which ctau < 10 mm',


  // *******************************************************************
  // Detailed process selection 


  pythia->SetMSTP(14,MSTP14_val);  // gamma gamma -> hadrons (30=FULL)
  //pythia->SetMSEL(2); // min-bias QCD
  pythia->SetMSEL(1); //less accuracy in hadronic simulation but faster(usefull if intereted in leptonic)


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

  std::ofstream myfile;
  myfile.open("_txt/generation_time.txt");
  int execution_time[Nevts/1000];

  Int_t n_lepton_total=0; //counter for the total number of charged leptons (sum on the events)

  for (int ievt = 0; ievt < Nevts; ievt++) {

    if (ievt ==0) cout << "Generating event #: " << flush;
    if (ievt % 1000 == 0) {
      cout << ievt << " " << flush;
      execution_time[ievt/1000]=time(NULL);
      if (ievt !=0)
	myfile << execution_time[ievt/1000] << "\t" << ievt << "\t" << 1000/(float)(execution_time[ievt/1000]-execution_time[ievt/1000 -1])<< endl;
    }

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

    TMCParticle *part;
    TParticle *partic;

    Int_t n_leptons = 0;

    for (int ip = 0; ip < N; ip++ ) {

      part = (TMCParticle*)particlesArray->At(ip);
      int status = part->GetKS();
      if (status!=1) continue; // Only consider final-state particles

      int pdg = part->GetKF();
      double charge = PDG->GetParticle(pdg)->Charge();
      if(charge==0) continue; // only charged particles

      if ( abs(pdg)!=11 && abs(pdg)!=13 && abs(pdg)!=15) continue;// only leptons
      

      partic = (TParticle*)TMCParticle2TParticle(part);
      double ptPartic = partic->Pt();
      double etaPartic = partic->Eta();
      //phiPartic = partic->Phi();

      //if ( ptPartic<0.15 ) continue; //cut on Pt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //if (abs(etaPartic)>1.5) continue; //Fill only if |eta|<1.5!!!!!!!!!!!!!!!!!!!!!   


      // Histo filling
      h_id_part->Fill(abs(pdg));
      hdsigmadeta->Fill(etaPartic);
      
      hdsigmadpT->Fill(ptPartic);

      n_leptons++; //Count the total number of charged hadrons in this event
      delete partic;
      // if (TMath::Abs(etaPartic)<1.) Nch++;
    } // End loop over all particles in event

    h_lep_per_ev->Fill(n_leptons);

    n_lepton_total += n_leptons;
    
    
  } // END EVENT GENERATION LOOP

  myfile.close();

  
  
  // **********************************************************************************  
  // FINAL GENERATION MESSAGES:
  
  pythia->Pystat(1);
  
  double xsection = pythia->GetPARI(1) * 1e9;//conversion form mb to pb
  int    ntrials  = pythia->GetMSTI(5);
  double sigmaweight = xsection/ntrials;
  double total_lepton_xsection = sigmaweight*n_lepton_total;

    
  cout << endl << "#######################" << endl 
       << endl << "<I> Event summary: " << Nevts << " generated." << endl ;

  
  //double evts = hdsigmadeta->Integral(); 
  double triggEff = (Nevts-exclEvts)/(float)Nevts;
  cout << "<I> Num. of total valid events = " << Nevts-exclEvts << " (Effic="<< triggEff*100. << "%)" << endl;

  cout << "<I> Pythia cross section: " << xsection << " mb || number of trials: " << ntrials << " || sigmaweight = "<< sigmaweight <<endl;

  cout << endl << "Hadrons per events: " << endl << "Mean = " << h_lep_per_ev->GetMean() << endl << "Max = " << h_lep_per_ev->GetMaximumBin() <<endl;

  cout << endl << endl << endl <<  "Total number of leptons = " << n_lepton_total << endl;
  cout << "Cross section = " << total_lepton_xsection << "pb  in e+e- --> gamma gamma --> X at sqrt(s) = " << sqrts << " GeV " << endl << endl << endl;

  Double_t N_rate = total_lepton_xsection*Luminosity;
  Double_t Pileup = (N_rate*FCC_Circumference)/(N_bunch*speed_of_light);
  Double_t N_rate_c = total_lepton_xsection*Luminosity_c;
  Double_t Pileup_c = (N_rate_c*FCC_Circumference)/(N_bunch_c*speed_of_light);

  ofstream file_out;
  char file_out_name[150];
  sprintf(file_out_name, "_txt/pythia6_gammagamma_leptons_%iGeV_seed%d_Nevts%d_output.txt",(int)TMath::Ceil(sqrts),seed,Nevts);
  file_out.open(file_out_name);

  file_out << Nevts-exclEvts << " //   Total event generated" << endl;
  file_out << xsection << " //   [pb] Pythia6 total cross section" << endl;
  file_out << total_lepton_xsection << " //   [pb] Cross section for e+e- --> gamma gamma --> l"<< endl;
  file_out << sqrts << " //   [GeV] sqrt(s)" << endl;
  file_out << n_lepton_total <<" //   Total Lepton produced" << endl;
  file_out << N_rate << " //   [Hz] Production Rate (Sigma*Lum)" << endl;
  file_out << N_rate_c << " //   [Hz] Production Rate (Sigma*Lum) (Crab waist)" << endl << endl << endl;
  
  file_out << Pileup << " //   Interactions e/gamma e/gamma -> l per bunch" << endl;
  file_out << Pileup_c << " //   Interactions e/gamma e/gamma -> l per bunch (Crab Waist)" << endl;

  
  file_out.close();



  // double dNchdeta = hdsigmadeta->GetBinContent(11)/ntrials; 
  // cout << "<I> dN_ch/deta|_{eta=0} = " << dNchdeta << " in e+e- --> gamma gamma --> X at sqrt(s) = " << sqrts << " GeV " << endl ;

  // **********************************************************************************  
  // Normalize histos by weighted cross-section, pseudorapidity range, and the pT bin size

  hdNdeta=(TH1F*)hdsigmadeta->Clone("hdNdeta");
  hdNdeta->Scale(1./(float)Nevts);
  hdNdeta->SetTitle("hdNdeta");
  hdNdeta->SetXTitle("|#eta|");
  hdNdeta->SetYTitle("dN/d|#eta| (nb)");

  double etabinsize = 20/20; // eta binning: 20 within -10<eta<10
  hdsigmadeta->Scale(1e-3*sigmaweight/etabinsize);
  //hdsigmadeta->Scale(1/Nevts);
  double ptbinsize = (pt_max-pt_min)/n_bin_pt;
  hdsigmadpT->Scale(sigmaweight/ptbinsize);
    
  //hdsigmadetaTruth->Scale(sigmaweight*0.01); // eta binning: 2000 in -10<eta<10
  //hEdsigmadpT->Scale(ptbinsize/ntrials);
  //ntFFdsigmadeta->SetWeight(sigmaweight);


  // **********************************************************************************  
  // Plot distributions for cross-check


  char title[300];
    sprintf(title, "cinvdsigmadpT_%iGeV",(int)sqrts);
    TCanvas *cinvdsigmadpT = new TCanvas(title,title,700,600);
    cinvdsigmadpT->SetLogy();
    cinvdsigmadpT->SetLogx();
    cinvdsigmadpT->cd();
    hdsigmadpT->Draw();
    cinvdsigmadpT->SaveAs("_png/histo_hdsigmadpT.png");   
    

    sprintf(title, "cinvdsigmadeta_%iGeV",(int)sqrts);
    TCanvas *cinvdsigmadeta = new TCanvas(title,title,700,600);
    cinvdsigmadeta->cd();
    hdsigmadeta->Draw();
    cinvdsigmadeta->SaveAs("_png/histo_hdsigmadeta.png");

    sprintf(title, "cinvdNdeta_%iGeV",(int)sqrts);
    TCanvas *cinvdNdeta = new TCanvas(title,title,700,600);
    cinvdNdeta->cd();
    hdNdeta->Draw();
    cinvdNdeta->SaveAs("_png/histo_hdNdeta.png");


    // sprintf(title, "cinvW_%iGeV",(int)sqrts);
    // TCanvas *cinvW = new TCanvas(title,title,700,600);
    // cinvW->cd();
    // hW->Draw();
    // cinvW->SaveAs("_png/histo_hW.png");

    sprintf(title, "cinvn_had_%iGeV",(int)sqrts);
    TCanvas *cinvn_had = new TCanvas(title,title,700,600);
    cinvn_had->cd();
    h_lep_per_ev->Draw();
    cinvn_had->SaveAs("_png/histo_charged_leptons_per_ev.png");

    // sprintf(title, "cinvp_id_%iGeV",(int)sqrts);
    // TCanvas *cinvp_id = new TCanvas(title,title,700,600);
    // cinvp_id->cd();
    // cinvp_id->SetLogy();
    // h_id_part->Draw();
    // cinvp_id->SaveAs("_png/histo_particle_id.png");

    

  // **********************************************************************************  
  // Open  output file and Close file
    
  char filename[200];
  sprintf(filename, "_root/pythia6_gammagamma_leptons_%iGeV_seed%d_Nevts%d.root",(int)TMath::Ceil(sqrts),seed,Nevts);

  TFile* file = TFile::Open(filename, "RECREATE");
  if (!file || !file->IsOpen()) {
    Error("pythia6_gammagamma_leptons", "Couldn;t open file %s", filename);
    return;
  }

  file->cd();
  hdsigmadeta->Write();
  hdsigmadpT->Write();
  hdNdeta->Write();
  h_lep_per_ev->Write();
  hW->Write();
  h_id_part->Write();
  file->Close();
  cout << endl << "#######<I> File " << filename << " created. Take a look ... ##############" << endl << endl;

  file = TFile::Open(filename);
  file->ls();
  file->Close();
}


//____________________________________________________________________
//
