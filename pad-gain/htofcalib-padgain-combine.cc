#include "../TPCPadHelper_260416.hh"

const int runnumber[8] = {2489,2599, 2601, 2602, 2603, 2604, 2606, 2607};
struct FitResult {
  TF1  *f = nullptr;
  TH1D *h = nullptr;
  //int rebin = 1;
};

bool IsGoodLandauFit(TH1D* h, TF1* f)
{
  if(!h || !f) return false;
  if(h->GetEntries() < 80) return false;
  if(f->GetNDF() <= 0) return false;
  
  double chi2ndf = f->GetChisquare() / f->GetNDF();
  //if(chi2ndf > 5.0) return false;
  
  double mpv   = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  if(mpv < 50 || mpv > 300) return false;
  //if(sigma > 100) return false;

  return true;
}

FitResult FitLandauWithRetry(TH1D *hraw){
  std::vector<int> rebinList = {1, 2, 4};
  for(auto rebin : rebinList){
    TH1D *h = (TH1D*)hraw->Clone(Form("%s_tmp_rebin%d", hraw->GetName(), rebin));
    h->SetDirectory(0);
    h->Rebin(rebin);
    TF1 *f = new TF1(Form("f_%s", hraw->GetName()),"landau",100,700);
    f->SetParameter(1,200);
    f->SetParameter(2,33);
    f->SetLineColor(kRed);
    int status = h->Fit(f,"RQ0N");
    if(IsGoodLandauFit(h, f)){
      //h->GetListOfFunctions()->Add(f);
      //return {f, h, rebin};
      return {f, h};
    }
    delete f;
    delete h;
  }
  return {};
}

void htofcalib_padgain_combine(){
  gROOT->SetBatch(kTRUE);
  string outpdf = "result/260629/htofcalib-padgain-combine-260629.pdf";
  
  TH1D *hist_de[NumOfPadTPC];

  TH2Poly *TPC_count = new TH2Poly("TPC_count","TPC_count;Z;X",MinZ,MaxZ,MinX,MaxX);
  auto TPC_gain = new TH2Poly("TPC_gain","TPC_gain;Z;X",MinZ,MaxZ,MinX,MaxX);
  TGraph *graph_gain = new TGraph();
  TGraph *graph_chi = new TGraph();
  TGraph *graph_sigma = new TGraph();
  graph_gain->SetName("graph_gain");
  graph_chi->SetName("graph_chi");
  graph_sigma->SetName("graph_sigma");
  TH1D *hist_chi = new TH1D("hist_chi","hist_chi;#chi^2/NDF;Counts",500,0,10);
  TH1D *hist_sigma = new TH1D("hist_sigma","hist_sigma;#sigma;Counts",100,0,200);

  double l = (586./2.)/(1+sqrt(2.));
  Double_t px[9] = {-l*(1+sqrt(2.)),-l,l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.))};
  Double_t py[9] = {l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.)),-l,l};
  
  auto pLine = new TPolyLine(9,px,py);

  Double_t X[5];
  Double_t Y[5];
  for (Int_t l=0; l<NumOfLayersTPC; ++l) {
    Double_t pLength = tpc::padParameter[l][5];
    Double_t st      = (180.-(360./tpc::padParameter[l][3]) *
                        tpc::padParameter[l][1]/2.);
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./tpc::padParameter[l][3])/180.*TMath::Pi();
    Double_t cRad    = tpc::padParameter[l][2];
    Int_t    nPad    = tpc::padParameter[l][1];
    for (Int_t j=0; j<nPad; ++j) {
      X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[0] = X[4];
      Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[0] = Y[4];
      for (Int_t k=0; k<5; ++k) X[k] += ZTarget;
      TPC_count->AddBin(5, X, Y);
      TPC_gain->AddBin(5, X, Y);
    }
  }
  
  for(int i=0;i<8;i++){
    TFile *file = new TFile(Form("result/260629/htofcalib-padgain-run0%d-260629.root",runnumber[i]));
    TH2Poly *TPC_tr_cluster = (TH2Poly*)file->Get("TPC_tr_cluster");
    
    for(int n=0;n<TPC_tr_cluster->GetNumberOfBins();n++){
      double run_cnt = TPC_tr_cluster->GetBinContent(n+1);

      double total_cnt = TPC_count->GetBinContent(n+1);
      TPC_count->SetBinContent(n+1,total_cnt + run_cnt);

      auto hist = (TH1D*)file->Get(Form("hist_de%d",n));
      if(i==0){
        hist_de[n] = (TH1D*)hist->Clone(Form("hist_de%d",n));
	hist_de[n] ->SetDirectory(0);
      }
      else{
        hist_de[n]->Add(hist);
      }
    }
  }


  TFile *f = new TFile("result/260629/htofcalib-padgain-combine-260629.root","RECREATE");
  
  auto c1 = new TCanvas("c1","c1");
  gStyle->SetOptStat(0);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("htofcalib-padgain-combine.cc");
  p->AddText("TPC Pad Gain Calibration");
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());




  for(int ipad = 0;ipad <NumOfPadTPC;ipad++){

    if(ipad%30 == 0){
      if(ipad !=0)c1->Print(outpdf.c_str());
      c1->Clear();
      c1->Divide(6,5);
    }
    c1->cd(ipad%30+1);
    FitResult res = FitLandauWithRetry(hist_de[ipad]);
    if(res.f && res.h){
      auto fit_param = res.f->GetParameters();
      double chi2ndf = res.f->GetChisquare() / res.f->GetNDF();
      hist_chi->Fill(chi2ndf);
      hist_sigma->Fill(fit_param[2]);
      graph_gain->AddPoint(ipad, fit_param[1]);
      graph_chi->AddPoint(ipad, chi2ndf);
      graph_sigma->AddPoint(ipad, fit_param[2]);

      TPC_gain->SetBinContent(ipad+1, fit_param[1]);
      res.h->Draw();
      res.f->Draw("same");
      //res.h->Write(hist_de[ipad]->GetName());
      res.h->Write(Form("hist_de%d",ipad));
      //res.f->Write(Form("%s_fitfunc", hist_de[ipad]->GetName()));

    }
    
      else{
	hist_de[ipad]->Draw();
	hist_de[ipad]->Write();
      }
    
  }

  c1->Print(outpdf.c_str());
  
  c1->Clear();
  
  TPC_gain->SetMinimum(0);
  TPC_gain->SetMaximum(400);
  
  TPC_gain->Draw("colz");
  TPC_gain->Write();
  c1->Print(outpdf.c_str());

  c1->Clear();
  auto mg = new TMultiGraph();
  graph_gain->SetMarkerStyle(4);
  mg->Add(graph_gain);
  graph_gain->Write();

  graph_chi->SetMarkerStyle(4);
  graph_chi->SetMarkerColor(kRed);
  mg->Add(graph_chi);
  graph_chi->Write();

  graph_sigma->SetMarkerStyle(4);
  graph_sigma->SetMarkerColor(kBlue);
  mg->Add(graph_sigma);
  graph_sigma->Write();
  mg->Draw("AP");
  c1->Print(outpdf.c_str());

  
  
  c1->Clear();
  gPad->SetLogz();
  TPC_count->Draw("colz");
  TPC_count->Write();
  c1->Print((outpdf + ")").c_str());
  hist_chi->Write();
  hist_sigma->Write();
  f->Close();
}
