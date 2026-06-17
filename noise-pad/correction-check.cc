#include "../TPCPadHelper_260416.hh"
const int total_hist = 50;

void correction_check(){
  gROOT->SetBatch(kTRUE);
  string dir = "/home/had/haein/data/JPARC2025Nov_root/noise-check-hist";
  string outpdf = "/home/had/haein/work/files/tpc-calib/noise-pad/correction-check-run02343.pdf";

  auto c1 = new TCanvas("c1","c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("correction-check.cc");
  p->AddText("TPC Pad Baseline After Correction");
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());
  c1->Clear();
  
  int pads_per_hist = (NumOfPadTPC + total_hist -1) / total_hist;
  
  for(int i=0;i<50;i++){
    TFile *file = new TFile(Form("%s/run02343_TPCHit_%d.root",dir.c_str(),i));
    int start_pad = i*pads_per_hist;
    int end_pad = std::min((i+1)*pads_per_hist, NumOfPadTPC);
    for(int j=start_pad;j<end_pad;j++){
      if(j%30==0){
	if(j!=0)c1->Print(outpdf.c_str());
	c1->Clear();
	c1->Divide(6,5);
      }
      c1->cd(j%30+1);
      auto h = (TH2D*)file->Get(Form("TPC_FADC_After_Pad%d",j));
      if(h){
	h->SetName(Form("h_fadc_cor%d",j));
	h->SetTitle(Form("h_fadc_cor%d",j));
	h->GetYaxis()->SetRangeUser(-200,1000);
	gPad->SetTopMargin(0.08);
	h->Draw("colz");
      }
    }
    
  }
  c1->Print((outpdf + ")").c_str());
  
}
