#include "../TPCPadHelper_260416.hh"

const int runnumber = 2489;

void pedestal_shape(){
  gROOT->SetBatch(kTRUE);
  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/gain_calib_noiseoff_202401";
  TFile *file = new TFile(Form("%s/run0%d_DstTPCHelixTracking.root",dir.c_str(),runnumber));
  TTree *tree = (TTree*)file->Get("tpc");

  int nclTpc;
  vector<double>* cluster_de = nullptr;
  vector<double>* cluster_layer = nullptr;
  vector<double>* cluster_x_center = nullptr;
  vector<double>* cluster_y_center = nullptr;
  vector<double>* cluster_z_center = nullptr;
  vector<double>* cluster_row_center = nullptr;

  int ntTpc;
  vector<int>* nhtrack = nullptr;
  vector<vector<double>>* hitlayer = nullptr;
  vector<vector<double>>* theta_diff = nullptr;
  vector<vector<double>>* track_cluster_de = nullptr;
  vector<vector<double>>* track_cluster_y_center = nullptr;
  vector<vector<double>>* track_cluster_size = nullptr;
  vector<vector<double>>* track_cluster_mrow = nullptr;
  vector<vector<double>>* track_cluster_row_center = nullptr;
  

  tree->SetBranchAddress("nclTpc",&nclTpc);
  tree->SetBranchAddress("cluster_de",&cluster_de);
  tree->SetBranchAddress("cluster_layer",&cluster_layer);
  tree->SetBranchAddress("cluster_x_center",&cluster_x_center);
  tree->SetBranchAddress("cluster_y_center",&cluster_y_center);
  tree->SetBranchAddress("cluster_z_center",&cluster_z_center);
  tree->SetBranchAddress("cluster_row_center",&cluster_row_center);

  tree->SetBranchAddress("ntTpc",&ntTpc);
  tree->SetBranchAddress("nhtrack",&nhtrack);
  tree->SetBranchAddress("hitlayer",&hitlayer);
  tree->SetBranchAddress("theta_diff",&theta_diff);
  tree->SetBranchAddress("track_cluster_de",&track_cluster_de);
  tree->SetBranchAddress("track_cluster_y_center",&track_cluster_y_center);
  tree->SetBranchAddress("track_cluster_size",&track_cluster_size);
  tree->SetBranchAddress("track_cluster_mrow",&track_cluster_mrow);
  tree->SetBranchAddress("track_cluster_row_center",&track_cluster_row_center);

  auto TPC_pattern1 = new TH2Poly("TPC_pattern1","TPC_pattern1;Z;X",MinZ,MaxZ,MinX,MaxX);
  auto TPC_pattern2 = new TH2Poly("TPC_pattern2","TPC_pattern2;Z;X",MinZ,MaxZ,MinX,MaxX);
  auto TPC_pattern3 = new TH2Poly("TPC_pattern3","TPC_pattern3;Z;X",MinZ,MaxZ,MinX,MaxX);
  auto TPC_pattern4 = new TH2Poly("TPC_pattern4","TPC_pattern4;Z;X",MinZ,MaxZ,MinX,MaxX);

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
      TPC_pattern1->AddBin(5, X, Y);
      TPC_pattern2->AddBin(5, X, Y);
      TPC_pattern3->AddBin(5, X, Y);
      TPC_pattern4->AddBin(5, X, Y);
    }
  }

  TH1D *hist_y[4];
  for(int i=0;i<4;i++){
    hist_y[i] = new TH1D(Form("hist_y%d",i),Form("hist_y%d",i),1000,-500,500);
  }
  TH1D *hist_de = new TH1D("hist_de","hist_de",1000,0,500);
  
  //for(int n = 0;n<tree->GetEntries();n++){
  for(int n = 0;n<10000;n++){
    tree->GetEntry(n);
    if(n%1000 == 0)cout<<n<<endl;

    for(int i=0;i<nclTpc;i++){
      int padid = tpc::GetPadId((*cluster_layer)[i],(*cluster_row_center)[i]);

      hist_de->Fill((*cluster_de)[i]);
      if((*cluster_de)[i] < 40){
	double cnt = TPC_pattern1->GetBinContent(padid+1);
	TPC_pattern1->SetBinContent(padid+1,cnt+1);
	hist_y[0]->Fill((*cluster_y_center)[i]);
      }
      else if((*cluster_de)[i] > 40 && (*cluster_de)[i] < 60){
	double cnt = TPC_pattern2->GetBinContent(padid+1);
	TPC_pattern2->SetBinContent(padid+1,cnt+1);
	hist_y[1]->Fill((*cluster_y_center)[i]);
	
      }
      else if((*cluster_de)[i] > 60 && (*cluster_de)[i] < 80){
	double cnt = TPC_pattern3->GetBinContent(padid+1);
	TPC_pattern3->SetBinContent(padid+1,cnt+1);
	hist_y[2]->Fill((*cluster_y_center)[i]);
      }
      else if((*cluster_de)[i] > 100){
	double cnt = TPC_pattern4->GetBinContent(padid+1);
	TPC_pattern4->SetBinContent(padid+1,cnt+1);
	hist_y[3]->Fill((*cluster_y_center)[i]);
	
      }
    }
  }

  auto c1 = new TCanvas("c1","c1");

  hist_de->Draw();
  c1->Print("pedestal-shape.pdf(");

  c1->Clear();
  gPad->SetLogz();
  TPC_pattern1->Draw("colz");
  c1->Print("pedestal-shape.pdf");

  c1->Clear();
  gPad->SetLogz();
  TPC_pattern2->Draw("colz");
  c1->Print("pedestal-shape.pdf");

  c1->Clear();
  gPad->SetLogz();
  TPC_pattern3->Draw("colz");
  c1->Print("pedestal-shape.pdf");

  c1->Clear();
  gPad->SetLogz();
  TPC_pattern4->Draw("colz");
  c1->Print("pedestal-shape.pdf");

  c1->Clear();
  c1->Divide(2,2);
  for(int i=0;i<4;i++){
    c1->cd(i+1);
    hist_y[i]->Draw();
  }
  c1->Print("pedestal-shape.pdf)");
}
