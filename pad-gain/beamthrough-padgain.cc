#include "../TPCPadHelper.hh"

const int runnumber = 2489;

void beamthrough_padgain(){
  gROOT->SetBatch(kTRUE);
  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/gain_calib_noiseoff_202401_minlayer5";
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
  

  auto TPC_cluster = new TH2Poly("TPC_cluster","TPC_cluster;Z;X",MinZ,MaxZ,MinX,MaxX);
  auto TPC_tr_cluster = new TH2Poly("TPC_tr_cluster","TPC_tr_cluster;Z;X",MinZ,MaxZ,MinX,MaxX);
  auto TPC_gain = new TH2Poly("TPC_gain","TPC_gain;Z;X",MinZ,MaxZ,MinX,MaxX);
  TGraph *graph_gain = new TGraph();
  graph_gain->SetName("graph_gain");
  TH1D *hist_de[NumOfPadTPC];
  TH1D *hist_size = new TH1D("hist_size","hist_size; Cluster Size;Counts",5,-0.5,4.5);
  TH1D *hist_inner_alpha = new TH1D("hist_inner_alpha","hist_inner_alpha;#alpha [rad];Counts",100,-1,1);
  TH1D *hist_outer_alpha = new TH1D("hist_outer_alpha","hist_outer_alpha;#alpha [rad];Counts",100,-1,1);
  TH1D *hist_y_total = new TH1D("hist_y_total","hist_y_total;Y [mm];Counts",600,-300,300);
  TH1D *hist_y = new TH1D("hist_y","hist_y;Y [mm];Counts",600,-300,300);
  
  for(int i=0;i<NumOfPadTPC;i++){
    hist_de[i] = new TH1D(Form("hist_de%d",i),Form("hist_de%d",i),100,0,1000);
  }
  
  
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
      TPC_cluster->AddBin(5, X, Y);
      TPC_tr_cluster->AddBin(5, X, Y);
      TPC_gain->AddBin(5,X,Y);
    }
  }

  bool y_cut = false;
  bool alpha_cut = false;
  bool de_cut = false;
  bool cluster_cut = false;
  
  //for(int n = 0;n<tree->GetEntries();n++){
  for(int n = 0;n<100000;n++){
    tree->GetEntry(n);
    
    
    
    if(n%1000 == 0)cout<<n<<endl;
    for(int i=0;i<nclTpc;i++){
      int padid = tpc::GetPadId((*cluster_layer)[i],(*cluster_row_center)[i]);
      //if((*cluster_de)[i]>100 && (*cluster_y_center)[i] > -50 && (*cluster_y_center)[i] < 50){
      double cnt = TPC_cluster->GetBinContent(padid+1);
	TPC_cluster->SetBinContent(padid+1,cnt+1);
	//}
    }
    for(int i=0;i<ntTpc;i++){
      for(int j=0;j<(*nhtrack)[i];j++){
	//init.
	y_cut = false;
	alpha_cut = false;
	de_cut = false;
	cluster_cut = false;
	int padid = tpc::GetPadId((*hitlayer)[i][j],(*track_cluster_row_center)[i][j]);
	
	if((*track_cluster_y_center)[i][j]>-50 && (*track_cluster_y_center)[i][j]<50)y_cut = true;
	hist_y_total->Fill((*track_cluster_y_center)[i][j]);
	if(!tpc::Noise(padid)){
	  hist_y->Fill((*track_cluster_y_center)[i][j]);
	}
	
	if((*hitlayer)[i][j]<10){ //inner layer
	  if(TMath::Abs((*theta_diff)[i][j])<0.1)alpha_cut = true;
	  else if(TMath::Abs((*theta_diff)[i][j])>TMath::Pi()-0.1 && TMath::Abs((*theta_diff)[i][j])<TMath::Pi()+0.1)alpha_cut = true;
	}
	else if((*hitlayer)[i][j]>=10){ //outer layer
	  if(TMath::Abs((*theta_diff)[i][j])<0.2)alpha_cut = true;
	}
	if((*track_cluster_de)[i][j]>80)de_cut = true;
	if((*track_cluster_size)[i][j] == 1)cluster_cut = true;
	hist_size->Fill((*track_cluster_size)[i][j]);
	if((*hitlayer)[i][j]<=9 && (*hitlayer)[i][j]>=0)
	  hist_inner_alpha->Fill((*theta_diff)[i][j]);
	else if((*hitlayer)[i][j]<=29 && (*hitlayer)[i][j]>=10)
	  hist_outer_alpha->Fill((*theta_diff)[i][j]);
	//if(y_cut && alpha_cut && de_cut){
	if(y_cut && alpha_cut && cluster_cut){
	  double cnt = TPC_tr_cluster->GetBinContent(padid+1);
	  TPC_tr_cluster->SetBinContent(padid+1,cnt+1);
	  hist_de[padid]->Fill((*track_cluster_de)[i][j]);
	  
	}
      }
    }
      
  }

  TFile *f = new TFile("beamthrough-padgain.root","RECREATE");
  
  auto c1 = new TCanvas("c1","c1");
  gPad->SetLogz();
  TPC_tr_cluster->Draw("colz");
  TPC_tr_cluster->Write();
  c1->Print("beamthrough-padgain.pdf(");

  for(int ipad = 0;ipad <NumOfPadTPC;ipad++){
    TF1 fL("fL", "landau", 0, 1000);
    fL.SetParLimits(1,50,300);
    fL.SetLineColor(kRed);
    
    if(ipad%30 == 0){
      if(ipad !=0)c1->Print("beamthrough-padgain.pdf");
      c1->Clear();
      c1->Divide(6,5);
    }
    c1->cd(ipad%30+1);
    if(hist_de[ipad]->GetEntries()>100){
      hist_de[ipad]->Fit(&fL, "QR");
      auto fit_param = fL.GetParameters();
      graph_gain->AddPoint(ipad,fit_param[1]);
      TPC_gain->SetBinContent(ipad+1,fit_param[1]);
    }

    hist_de[ipad]->Draw();
    hist_de[ipad]->Write();
  }
  
  c1->Print("beamthrough-padgain.pdf");
  c1->Clear();
  
  TPC_cluster->Draw("colz");
  c1->Print("beamthrough-padgain.pdf");
  c1->Clear();
  hist_size->Draw();
  hist_size->Write();
  c1->Print("beamthrough-padgain.pdf");
  c1->Clear();
  hist_inner_alpha->Draw();
  hist_inner_alpha->Write();
  c1->Print("beamthrough-padgain.pdf");
  c1->Clear();
  hist_outer_alpha->Draw();
  hist_outer_alpha->Write();
  c1->Print("beamthrough-padgain.pdf");
  c1->Clear();
  hist_y_total->Draw();
  hist_y_total->Write();
  hist_y->SetLineColor(kRed);
  hist_y->Draw("same");
  hist_y->Write();
  c1->Print("beamthrough-padgain.pdf");
  c1->Clear();
  
  graph_gain->SetMarkerStyle(4);
  graph_gain->Draw("AP");
  graph_gain->Write();
  c1->Print("beamthrough-padgain.pdf");

  c1->Clear();
  TPC_gain->SetMinimum(0);
  TPC_gain->SetMaximum(300);
  TPC_gain->Draw("colz");
  TPC_gain->Write();
  c1->Print("beamthrough-padgain.pdf)");

  f->Close();
  
  
}
