const int NCobo = 8;
const int runnumber = 3772;
const bool param_update = true;
Double_t ClockShiftFunc(Double_t *x, Double_t *par)
{
  return par[0] * TMath::Freq((x[0] - par[1]) / par[2])+par[3];
}

TH2D* CorrectHistY(TH2D *h, TF1 *fshift, const char *name)
{
  TH2D *hcorr = new TH2D(name, h->GetTitle(),
                         h->GetNbinsX(),
                         h->GetXaxis()->GetXmin(),
                         h->GetXaxis()->GetXmax(),
                         h->GetNbinsY(),
                         h->GetYaxis()->GetXmin(),
                         h->GetYaxis()->GetXmax());

  for(int ix = 1; ix <= h->GetNbinsX(); ix++){

    double x = h->GetXaxis()->GetBinCenter(ix);
    double shift = fshift->Eval(x);

    for(int iy = 1; iy <= h->GetNbinsY(); iy++){

      double y = h->GetYaxis()->GetBinCenter(iy);
      double cont = h->GetBinContent(ix, iy);

      if(cont <= 0) continue;

      double ycorr = y - shift;

      hcorr->Fill(x, ycorr, cont);
    }
  }

  return hcorr;
}

void UpdateCoboParameter(const char* infile,
                         const char* outfile,
			 const double p0fit[],
                         const double p1fit[],
                         const int NCobo)
{
  std::ifstream fin(infile);
  std::ofstream fout(outfile);

  std::string line;
  bool inCobo = false;
  int idx = 0;

  while(std::getline(fin, line)){

    if(line.find("Coboparameter") != std::string::npos){
      inCobo = true;
      idx = 0;
      fout << line << "\n";
      continue;
    }

    if(inCobo && idx >= NCobo){
      inCobo = false;
    }

    if(inCobo){
      std::stringstream ss(line);

      int cobo, dummy, aty;
      double oldp0, oldp1, oldp2;

      if(ss >> cobo >> dummy >> aty >> oldp0 >> oldp1 >> oldp2){

        fout << std::setw(8)  << cobo
             << std::setw(8)  << dummy
             << std::setw(8)  << aty
             << std::setw(16) << std::fixed << std::setprecision(8) << -80
             << std::setw(16) << std::fixed << std::setprecision(8) << p1fit[idx]
             << std::setw(16) << std::fixed << std::setprecision(8) << oldp2
             << "\n";

        std::cout << "update cobo idx " << idx
                  << " : " << oldp0 << " " << oldp1
                  << " -> " << p0fit[idx] << " " << p1fit[idx]
                  << std::endl;

        idx++;
        continue;
      }
    }

    fout << line << "\n";
  }

  std::cout << "Total updated Cobo lines = " << idx << std::endl;
}


void cobo_clock(){
  gROOT->SetBatch(kTRUE);
  TFile *file = TFile::Open(Form("~/data/JPARC2025Nov_root/cobo-phase/run0%d_DstTPCHitBcOutTracking.root",runnumber),"read");
  TString histFmt = "TPCHit_ResY_vs_ClockTime_CoBo%d_RawClock";
  TString histFmt_cor = "TPCHit_ResY_vs_ClockTime_CoBo%d";

  string outpdf = Form("result/cobo-clock-run0%d.pdf", runnumber);
  TFile *fout = new TFile(Form("result/cobo-clock-run0%d.root",runnumber),"RECREATE");
  
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
  c1->Divide(4,2);
  
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(4,2);
  TF1 *fshift[NCobo];
  TH2D *h2[NCobo];
  TH2D *h2_cor[NCobo];
  double p1[NCobo]={0.};
  if(runnumber == 2448 || runnumber == 2447){
    p1[0] = 19.5;
    p1[1] = 18.8;
    p1[2] = 18.5;
    p1[3] = -30.7;
    p1[4] = -11.8;
    p1[5] = 18.0;
    p1[6] = 17.0;
    p1[7] = 18.6;
  }
  else if(runnumber == 2570){
    p1[0] = 11.;
    p1[1] = 19.5;
    p1[2] = 19.2;
    p1[3] = -19.5;
    p1[4] = -20.2;
    p1[5] = 19.5;
    p1[6] = 18.9;
    p1[7] = 19.3;
  }
  else if(runnumber == 3772){
    p1[0] = 23.;
    p1[1] = 21.4;
    p1[2] = 21.2;
    p1[3] = 22.9;
    p1[4] = 20.2;
    p1[5] = 21.0;
    p1[6] = 20.5;
    p1[7] = 21.3;
  }
  double p0_fit[NCobo];
  double p1_fit[NCobo];
  for (int icobo = 0; icobo < NCobo; ++icobo) {
    cout<<icobo<<endl;
    
    TString hname = Form(histFmt, icobo);
    h2[icobo] = (TH2D*)file->Get(hname);
    h2[icobo]->RebinY(2);

    int bin1 = h2[icobo]->GetXaxis()->FindBin(-40.0);
    int bin2 = h2[icobo]->GetXaxis()->FindBin(-35.0);
    int bin3 = h2[icobo]->GetXaxis()->FindBin(30.0);
    int bin4 = h2[icobo]->GetXaxis()->FindBin(40.0);
    
    TH1D *hleft = h2[icobo]->ProjectionY(Form("hleft%d",icobo),bin1,bin2);
    TH1D *hright = h2[icobo]->ProjectionY(Form("hright%d",icobo),bin3,bin4);

    int maxbin_left = hleft->GetMaximumBin();
    double max_left = hleft->GetBinCenter(maxbin_left);

    int maxbin_right = hright->GetMaximumBin();
    double max_right = hright->GetBinCenter(maxbin_right);
    
    TF1 *fleft = new TF1(Form("fleft%d",icobo),"gaus",max_left-5,max_left+5);
    TF1 *fright = new TF1(Form("fright%d",icobo),"gaus",max_right-5,max_right+5);

    hleft->Fit(fleft,"RQ0");
    hright->Fit(fright,"RQ0");
    double yleft = fleft->GetParameter(1);
    double yright = fright->GetParameter(1);
    double p0 = yright - yleft;

    hleft->Write();
    hright->Write();
    fshift[icobo] = new TF1(Form("fshift%d",icobo),ClockShiftFunc,-60,60,4);

    double p2 = 0.01;
    double p3 = yleft;
    fshift[icobo]->SetParameters(p0,p1[icobo],p2,p3);
    
    fshift[icobo]->SetParLimits(0,p0-0.1,p0+0.1);
    fshift[icobo]->SetParLimits(1,p1[icobo]-0.5,p1[icobo]+0.5);
    //fshift[icobo]->SetParLimits(2,0,p2+0.02);
    //fshift[icobo]->FixParameter(0,p0);
    fshift[icobo]->FixParameter(2,p2);
    fshift[icobo]->SetParLimits(3,p3-0.2,p3+0.2);
    
    c1->cd(icobo+1);
    h2[icobo]->Fit(fshift[icobo],"R");
    h2[icobo]->SetTitle(Form("CoBo%d",icobo));
    h2[icobo]->Draw("colz");
    
    h2[icobo]->Write();
    fshift[icobo]->Write();
    p0_fit[icobo] = fshift[icobo]->GetParameter(0) * -1 / 0.055;
    p1_fit[icobo] = fshift[icobo]->GetParameter(1);
    
  }
  //c1->Modified();
  //c1->Update();
  c1->Print(outpdf.c_str());
  c1->Clear();
  c1->Divide(4,2);
  for(int icobo=0;icobo<NCobo;icobo++){
    
    TH2D *hcorr = CorrectHistY(h2[icobo], fshift[icobo], Form("hcorr_cobo%d",icobo));
    c1->cd(icobo+1);
    hcorr->Draw("colz");
    hcorr->Write();

  }

  
  

  //c1->Modified();
  //c1->Update();
  c1->Print(outpdf.c_str());

  c1->Clear();
  c1->Divide(4,2);
  for(int icobo=0;icobo<NCobo;icobo++){

    TString hname = Form(histFmt_cor, icobo);
    h2_cor[icobo] = (TH2D*)file->Get(hname);
    h2_cor[icobo]->RebinY(2);
    c1->cd(icobo+1);
    h2_cor[icobo]->Draw("colz");
    h2_cor[icobo]->Write();
    
  }
  //c1->Modified();
  //c1->Update();
  c1->Print((outpdf + ")").c_str());


  
  fout->Close();

  if(param_update){
    UpdateCoboParameter("param_history/TPCParam_e72_20260616",
			Form("param_history/TPCParam_e72_run0%d",runnumber),
			p0_fit, p1_fit, NCobo);
  }
}
