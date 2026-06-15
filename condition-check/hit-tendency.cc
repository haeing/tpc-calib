#include "../TPCPadHelper_260416.hh"

const int start_run = 2263;
const int end_run = 3849;


void hit_tendency(){
  gROOT->SetBatch(kTRUE);

  TGraph *g_tendency[NumOfAsadTPC];
  
  for(int i=0;i<NumOfAsadTPC;i++){
    g_tendency[i] = new TGraph();
    g_tendency[i]->SetName(Form("g_asad%d",i));
    g_tendency[i]->SetTitle(Form("AsAd%d",i));
  }
  
  string outpdf = "result/hit-tendency.pdf";
  TFile *f = new TFile("result/hit-tendency.root","RECREATE");
  
  for(int run = start_run;run<=end_run;run++){
    TString fname = Form("result/tpc-condition-run0%d.root",run);
    if(gSystem->AccessPathName(fname)){
      cout << "No file : " << fname << endl;
      continue;
    }
    

    TFile *file = TFile::Open(fname);
    if(!file || file->IsZombie()){
      cout << "Cannot open : " << fname << endl;
      if(file) file->Close();
      continue;
    }
    for(int i=0;i<NumOfAsadTPC;i++){
      TH1D *h = (TH1D*)file->Get(Form("hist_hitnum%d",i));
      TF1 *fit = new TF1("fit", "[0]", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
      
      fit->SetParameter(0, h->GetMean());

      h->Fit(fit, "RQ");
      double value = fit->GetParameter(0);
      int n = g_tendency[i]->GetN();
      g_tendency[i]->SetPoint(n, run, value);
      delete fit;
    }
  }

  auto c1 = new TCanvas("c1","c1");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("hit-tendency.cc");
  p->AddText(Form("run0%d - 0%d",start_run,end_run));
  p->AddText("TPC AsAd condition check");
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());
  for(int i=0;i<NumOfAsadTPC;i++){
    c1->Clear();
    g_tendency[i]->SetMarkerStyle(4);
    g_tendency[i]->Draw("AP");
    f->cd();
    g_tendency[i]->Write();
    if(i==NumOfAsadTPC-1){c1->Print((outpdf + ")").c_str());}
    else{c1->Print(outpdf.c_str());}
    
  }

  f->Close();
  
  
}
