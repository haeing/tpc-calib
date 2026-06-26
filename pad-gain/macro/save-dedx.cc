const int runnumber[] = {2599, 2601, 2602, 2603, 2604, 2606, 2607};
const int nrun = sizeof(runnumber) / sizeof(runnumber[0]);

void save_dedx(){

  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0);
  
  TH2D *dedx_total = nullptr;
  for(int i=0;i<nrun;i++){
    TFile *file = new TFile(Form("~/data/JPARC2025Nov_root/gain_calib/run0%d_DstTPCHelixTracking.root",runnumber[i]));
    TH2D *dedx = (TH2D*)file->Get("PID_dEdx_vs_SignedMom");
    if(i==0)dedx_total = (TH2D*)dedx->Clone("PID_dEdx_vs_SignedMom");
    else{
      dedx_total->Add(dedx);
    }
  }

  dedx_total->GetYaxis()->SetTitle("TPC #LT#it{dE/dx}#GT (a.u.)");
  dedx_total->GetXaxis()->SetTitle("#it{p/z} (GeV/#it{c})");
  
  TFile *f = new TFile("save-dedx.root","RECREATE");
  dedx_total->Write();
  f->Close();
  
}
