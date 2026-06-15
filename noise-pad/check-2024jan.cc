#include "TPCPadHelper_2401.hh"

void check_2024jan()
{
  auto TPC_noise = new TH2Poly("TPC_noise","TPC_noise;Z;X",MinZ,MaxZ,MinX,MaxX);
  
  
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
      TPC_noise->AddBin(5, X, Y);
    }
  }
  /*
  for(int i=0;i<NumOfLayersTPC;i++){
    for(int j=0;j<42;j++){
      if(tpc::NoiseChannel[i][j] != -1){
	TPC_noise->SetBinContent(tpc::GetPadId(i,tpc::NoiseChannel[i][j]),1);
      }
    }
  }
  */

  
  

  /*
    int n = sizeof(tpc::padOnFrame) / sizeof(tpc::padOnFrame[0]);
  for(int i=0;i<n;i++){
    TPC_noise->SetBinContent(tpc::padOnFrame[i],1);
  }
  */

  //0 origin
  int n = sizeof(tpc::padOnCenterFrame) / sizeof(tpc::padOnCenterFrame[0]);
  for(int i=0;i<n;i++){
    TPC_noise->SetBinContent(tpc::padOnCenterFrame[i]+1,1);
  }


  //1 origin
  int n_pad = sizeof(tpc::Noise_Pad) / sizeof(tpc::Noise_Pad[0]);
  for(int i=0;i<n_pad;i++){
    TPC_noise->SetBinContent(tpc::Noise_Pad[i],1);
  }

  

    
  auto c1 = new TCanvas("c1","c1");
  TPC_noise->Draw("colz");
}
