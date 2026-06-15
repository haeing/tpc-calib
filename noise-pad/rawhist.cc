#include "../TPCPadHelper_260416.hh"

void rawhist(){
  const int runnumber = 2343;
  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/noise-check";
  TFile *file = new TFile(Form("%s/run0%d_TPCHit.root",dir.c_str(),runnumber));

  TH2D *hist_fadc[NumOfPadTPC];
  TH1D *hist_rms[NumOfPadTPC];

  for(int i=0; i<NumOfPadTPC; i++){
    hist_fadc[i] = (TH2D*)file->Get(Form("TPC_FADC_Pad%d",i));
    hist_rms[i]  = (TH1D*)file->Get(Form("TPC_RMS_Pad%d",i));
  }

  auto TPC_rms = new TH2Poly("TPC_rms","TPC_rms;Z;X",MinZ,MaxZ,MinX,MaxX);
  
  
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
      TPC_rms->AddBin(5, X, Y);
    }
  }

  for(int i=0; i<NumOfPadTPC; i++){
    if(!hist_rms[i]) continue;
    TF1 *fit = new TF1(Form("fit_rms_%d", i), "gaus",
                       hist_rms[i]->GetMean() - 3*hist_rms[i]->GetRMS(),
                       hist_rms[i]->GetMean() + 3*hist_rms[i]->GetRMS());
    hist_rms[i]->Fit(fit, "RQ0");
    double mean  = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    TPC_rms->SetBinContent(i+1, sigma);
  }

  TString pdfname =Form("run0%d_rawhist.pdf",runnumber);
  TCanvas *copen = new TCanvas();
  copen->Print(pdfname+"[");

  TCanvas *cmap = new TCanvas("cmap","cmap",900,800);
  TPC_rms->Draw("colz");
  cmap->Print(pdfname);

  TCanvas *c1 = new TCanvas("c1","c1",1800,1200);
  int ipad = 0;
  for(int i=0; i<NumOfPadTPC; i++){
    if(ipad==0){
      c1->Clear();
      c1->Divide(5,6);
    }
    c1->cd(ipad+1);
    if(hist_fadc[i]){
      gPad->SetLogz();
      hist_fadc[i]->SetTitle(Form("Pad %d",i));
      hist_fadc[i]->Draw("colz");
    }
    ipad++;
    if(ipad==30){
      c1->Print(pdfname);
      ipad=0;
    }
  }
  if(ipad!=0){
    c1->Print(pdfname);
  }
  TCanvas *cclose = new TCanvas();
  cclose->Print(pdfname+"]");
  
}
