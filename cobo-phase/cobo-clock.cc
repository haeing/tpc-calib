const int NCobo = 8;
const int runnumber = 2448;

Double_t ClockShiftFunc(Double_t *x, Double_t *par)
{
  return par[0] * TMath::Freq((x[0] - par[1]) / par[2]);
}


void cobo_clock(){

  TFile *file = TFile::Open("~/data/JPARC2025Nov_root/cobo-phase/run02448_DstTPCHitBcOutTracking.root","read");
  TString histFmt = "TPCHit_ResY_vs_ClockTime_CoBo%d_RawClock";

  string outpdf = Form("result/cobo-clock%d.pdf", runnumber);
  TFile *fout = new TFile("result/cobo-clock.root","RECREATE");
  
  TCanvas *c1 = new TCanvas("c1","c1");
  gStyle->SetOptStat(0);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("cobo-clock.cc");
  p->AddText("CoBo Clocktime Phaseshift Correction");
  p->AddText(Form("run0%d",runnumber));
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());
  c1->Clear();
  
  TF1 *fshift[NCobo];
  double p1[NCobo] = {19.5,18.8,18.5,-30.7,-11.8,18.0,17.0,18.6};
  for (int icobo = 0; icobo < NCobo; ++icobo) {
    cout<<icobo<<endl;
    
    TString hname = Form(histFmt, icobo);
    TH2D *h2 = (TH2D*)file->Get(hname);
    h2->RebinY(2);

    int bin1 = h2->GetXaxis()->FindBin(-40.0);
    int bin2 = h2->GetXaxis()->FindBin(-35.0);
    int bin3 = h2->GetXaxis()->FindBin(30.0);
    int bin4 = h2->GetXaxis()->FindBin(40.0);
    
    TH1D *hleft = h2->ProjectionY(Form("hleft%d",icobo),bin1,bin2);
    TH1D *hright = h2->ProjectionY(Form("hright%d",icobo),bin3,bin4);

    int maxbin_left = hleft->GetMaximumBin();
    double max_left = hleft->GetBinCenter(maxbin_left);

    int maxbin_right = hright->GetMaximumBin();
    double max_right = hright->GetBinCenter(maxbin_right);
    
    TF1 *fleft = new TF1(Form("fleft%d",icobo),"gaus",max_left-5,max_left+5);
    TF1 *fright = new TF1(Form("fright%d",icobo),"gaus",max_right-5,max_right+5);

    hleft->Fit(fleft,"R");
    hright->Fit(fright,"R");
    double yleft = fleft->GetParameter(1);
    double yright = fright->GetParameter(1);
    double p0 = yright-yleft;
    cout<<yright<<","<<yleft<<endl;
    hleft->Write();
    hright->Write();
    fshift[icobo] = new TF1(Form("fshift%d",icobo),ClockShiftFunc,-60,60,3);

    double p2 = 0.01;
    fshift[icobo]->SetParameters(p0,p1[icobo],p2);
    fshift[icobo]->SetParLimits(0,p0-0.2,p0+0.2);
    fshift[icobo]->SetParLimits(1,p1[icobo]-1,p1[icobo]+1);
    fshift[icobo]->FixParameter(2,p2);
    
    h2->Fit(fshift[icobo],"R");
    
    h2->Draw("colz");
    h2->Write();
    fshift[icobo]->Write();
    if(icobo==NCobo-1)c1->Print((outpdf + ")").c_str());
    else{c1->Print(outpdf.c_str());}
    c1->Clear();
    
  }
  
  fout->Close();
  
}
