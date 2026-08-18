const int runnumber = 2682;

void residual()
{
  gROOT->SetBatch(kTRUE);
  TFile *file = TFile::Open(Form("~/data/JPARC2025Nov_root/physics-755/run0%d_DstTPCK18HelixTracking.root",runnumber),"read");
  string outpdf = Form("result/residual-run0%d_yoffset0.pdf", runnumber);
  TCanvas *c1 = new TCanvas("c1","c1");
  
  gStyle->SetOptStat(0);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("residual.cc");
  p->AddText(Form("run0%d",runnumber));
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());
  c1->Clear();

  TString hname = "TPCTrk_ResY_vs_Layer_Trk";
  TH2D *h2 = (TH2D*)file->Get(hname);
  h2->GetXaxis()->SetTitle("Layer ID");
  h2->GetYaxis()->SetTitle("Residual Y");

  h2->Draw("colz");
  c1->Print((outpdf + ")").c_str());
  
}
