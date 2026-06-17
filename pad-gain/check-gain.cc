#include "../TPCPadHelper_260416.hh"

void check_gain(const char* param_in  = "param_history/TPCParam_e72_20260616")
{
  TH2Poly *hGain = new TH2Poly("hGain","hGain",MinZ,MaxZ,MinX,MaxX);

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
      hGain->AddBin(5, X, Y);
    }
  }

  ifstream fin(param_in);

  string line;
  int iline = 0;
  int ipad_count = 0;

  while(getline(fin, line)){
    iline++;

    // TPC ADC data starts from line 13
    if(iline >= 13 && ipad_count < NumOfPadTPC){

      stringstream ss(line);

      int layer, row, aty;
      double p0, p1;

      if(ss >> layer >> row >> aty >> p0 >> p1){

        if(aty == 0){
          int padid = tpc::GetPadId(layer, row);
	  cout<<padid+1<<" : "<<p1<<endl;
	  hGain->SetBinContent(padid+1, 200.0 / p1);
        }

        ipad_count++;
        continue;
      }
    }
    
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  gStyle->SetPadRightMargin(0.15);
  hGain->SetTitle(";Z [mm];X [mm]");
  hGain->Draw("colz");
  c1->SaveAs("result/260616/check-gain-260616.pdf");
  TFile *f = new TFile("result/260616/check-gain-260616.root","RECREATE");
  hGain->Write();
  f->Close();

  
  
  
}
