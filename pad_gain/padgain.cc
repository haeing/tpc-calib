#include "../TPCPadHelper_260416.hh"

void padgain(){
  TFile *f = new TFile("padgain.root");

  TH1D *hist_de[NumOfPadTPC];
  TH1D *hist_de_cut[NumOfPadTPC];


  TCanvas c("c","c",900,700);
  gStyle->SetOptStat(0);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("padgain.cc");
  p->AddText("TPC Pad Gain Calibration");
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c.Print("padgain_hist.pdf(");
  c.Clear();
  
  
  for(int i=0;i<NumOfPadTPC;i++){
    hist_de[i] = (TH1D*)f->Get(Form("hist_de%d",i));
    hist_de_cut[i] = (TH1D*)f->Get(Form("hist_de_cut%d",i));
    if(i%30==0){
      if(i!=0)c.Print("padgain_hist.pdf");
      c.Clear();
      c.Divide(6,5);
    }
    c.cd(i%30+1);
    hist_de[i]->Draw();
    hist_de_cut[i]->SetLineColor(kRed);
    hist_de_cut[i]->Draw("same");
  }
  c.Print("padgain_hist.pdf)");
}
