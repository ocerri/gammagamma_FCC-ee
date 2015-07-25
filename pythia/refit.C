#include "pythia6_gammagamma_hadrons.h"

void refit(){
  Float_t sqrts=160.0; //energia nel Cm in GeV
  TFile *file = new TFile("_root/MSTP14-10_final.root");
  TH1F *histo = (TH1F*)file->Get("hdsigmadpT");
  //file->Close();

  TF1 *f_fit = new TF1("f_fit", pdf, pt_min, pt_max, 3);
  f_fit->SetParNames("m","C","alpha");

  cout << "Starting fitting to dsigmadpT..." << endl << "Give starting parameters for pdf(x) =C*x^alpha +m" << endl;

  Double_t m = 0;
  //cout << "m :"; cin >> m;

  Double_t C = 2500;
  //cout << "C: "; cin >> C;

  Double_t alpha = -3;
  //cout << "alpha: "; cin >> alpha;
  
  f_fit->SetParameter(0, m);
  f_fit->FixParameter(0,0);
  f_fit->SetParameter(1, C);
  f_fit->SetParameter(2, alpha);

  histo->Fit(f_fit,"","",0.6,7);
  Double_t fit_chi2 = f_fit->GetChisquare();
  Int_t    fit_ndof = f_fit->GetNDF();
  Double_t fit_prob = f_fit->GetProb();
  cout << endl;
  cout << "chi2 = " << fit_chi2 << endl;
  cout << "ndof = " << fit_ndof << endl;
  cout << "prob = " << fit_prob << endl;
  cout << endl;

  std::ofstream myfile;
  myfile.open("_txt/final_fit_results_MSRP14-10.txt");
  myfile << "fitting with C*x^alpha" << endl;
  myfile << "C = " << f_fit->GetParameter(1) << " +/- " << f_fit->GetParError(1) << endl;
  myfile << "alpha = " << f_fit->GetParameter(2) << " +/- " << f_fit->GetParError(2) << endl;
  myfile << "chi2 = " << fit_chi2 << endl;
  myfile << "ndof = " << fit_ndof << endl;
  myfile << "prob = " << fit_prob << endl;
  myfile.close();

  char title[300];
  sprintf(title, "cinvdsigmadpT_%iGeV",(int)sqrts);
  TCanvas *cinvdsigmadpT = new TCanvas(title,title,700,600);
  cinvdsigmadpT->SetLogy();
  cinvdsigmadpT->SetLogx();
  cinvdsigmadpT->cd();
  histo->Draw();
  //cinvdsigmadpT->SaveAs("_png/histo_hdsigmadpT_final_MSTP14-10.png");


  //file = new TFile("_root/MSTP14-10_final.root");
  // file->cd();
  // histo->Write();
  // file->Close();

  
  TFile* file_out = TFile::Open("_root/MSTP14-10_final_out.root", "RECREATE");
  if (!file || !file->IsOpen()) {
    Error("pythia6_gammagamma_hadrons_final_out", "Couldn;t open file %s", "_root/MSTP14-10_final_out.root");
    return;
  }

  file_out->cd();
  histo->Write();
  file->Close();
  cout << endl << "#######<I> File " << "_root/MSTP14-10_final_out.root" << " created. Take a look ... ##############" << endl << endl;
    
}    
