#include "../TPCPadHelper_260416.hh"

const int runnumber = 2599;
const int cutflag = 8;
//Cut Flag
// 0 - 000 : no cut 
// 1 - 100 : cluster
// 2 - 010 : y
// 3 - 001 : mip
// 4 - 110 : cluster y
// 5 - 101 : cluster mip
// 6 - 011 : y mip
// 7 - 111 : all
void tpchelixtracking_padgain(){

  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/gain_calib_noiseoff_202401";
  TFile *file = new TFile(Form("%s/run0%d_DstTPCHelixTracking.root",dir.c_str(),runnumber));
  TTree *tree = (TTree*)file->Get("tpc");

  int nclTpc;
  int ntTpc;

  vector<double>* mom0 = nullptr;
  vector<double>* dEdx = nullptr;

  //Track Info
  vector<int>* nhtrack = nullptr;
  vector<vector<double>>* hitlayer = nullptr;
  vector<vector<double>>* theta_diff = nullptr;
  vector<vector<double>>* track_cluster_de = nullptr;
  vector<vector<double>>* track_cluster_y_center = nullptr;
  vector<vector<double>>* track_cluster_size = nullptr;
  vector<vector<double>>* track_cluster_mrow = nullptr;
  vector<vector<double>>* track_cluster_row_center = nullptr;
  
  
  tree->SetBranchAddress("ntTpc",&ntTpc);
  tree->SetBranchAddress("mom0",&mom0);
  tree->SetBranchAddress("dEdx",&dEdx);
  
  tree->SetBranchAddress("nhtrack",&nhtrack);
  tree->SetBranchAddress("hitlayer",&hitlayer);
  tree->SetBranchAddress("theta_diff",&theta_diff);
  tree->SetBranchAddress("track_cluster_de",&track_cluster_de);
  tree->SetBranchAddress("track_cluster_y_center",&track_cluster_y_center);
  tree->SetBranchAddress("track_cluster_size",&track_cluster_size);
  tree->SetBranchAddress("track_cluster_mrow",&track_cluster_mrow);
  tree->SetBranchAddress("track_cluster_row_center",&track_cluster_row_center);
  
  
  auto TPC_gain = new TH2Poly("TPC_gain", "TPC_gain;Z;X", MinZ, MaxZ, MinX, MaxX);
  TGraph *graph_gain = new TGraph();
  auto TPC_count = new TH2Poly("TPC_count", "TPC_count;Z;X", MinZ, MaxZ, MinX, MaxX);
  auto TPC_pid = new TH2D("TPC_pid","TPC_pid;p [GeV/#it{c}];dE/dx [A.U.]",100,0,1,500,0,500);
  auto TPC_y = new TH1D("TPC_y","TPC_y;Y [mm];Counts",800,-400,400);
  auto hist_dedx = new TH1D("hist_dedx","hist_dedx",500,0,100);

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
      TPC_gain->AddBin(5, X, Y);
      TPC_count->AddBin(5, X, Y);
    }
  }

  for(int i=0;i<NumOfPadTPC;i++){
    TPC_count->SetBinContent(i+1,0);
  }

  double de_sum[NumOfPadTPC]={0.};
  int count[NumOfPadTPC]={0};
  TH1D *hist_de[cutflag][NumOfPadTPC];
  
  int total[NumOfLayersTPC]={0};
  int hit_layer[NumOfLayersTPC]={0};
  int clu_size[NumOfLayersTPC]={0};
  TH1D *hist_clu = new TH1D("hist_clu","hist_clu;Cluster;Entries",5,-0.5,4.5);
  
  for(int i=0;i<NumOfPadTPC;i++){
    hist_de[0][i] = new TH1D(Form("hist_de_000_%d",i),Form("hist_de_000_%d",i),100,0,1000);
    hist_de[1][i] = new TH1D(Form("hist_de_100_%d",i),Form("hist_de_100_%d",i),100,0,1000);
    hist_de[2][i] = new TH1D(Form("hist_de_010_%d",i),Form("hist_de_010_%d",i),100,0,1000);
    hist_de[3][i] = new TH1D(Form("hist_de_001_%d",i),Form("hist_de_001_%d",i),100,0,1000);
    hist_de[4][i] = new TH1D(Form("hist_de_110_%d",i),Form("hist_de_110_%d",i),100,0,1000);
    hist_de[5][i] = new TH1D(Form("hist_de_101_%d",i),Form("hist_de_101_%d",i),100,0,1000);
    hist_de[6][i] = new TH1D(Form("hist_de_011_%d",i),Form("hist_de_011_%d",i),100,0,1000);
    hist_de[7][i] = new TH1D(Form("hist_de_111_%d",i),Form("hist_de_111_%d",i),100,0,1000);

  }

  bool clu_cut = false; //1
  bool y_cut = false; //-50 mm < y < 50 mm
  bool mip_cut = false; //dEdx < 30
  
  //for(int n=0;n<tree->GetEntries();n++){
  for(int n=0;n<100000;n++){
    tree->GetEntry(n);
    clu_cut = false;
    y_cut = false;
    mip_cut = false;
    if(n%1000==0)std::cout<<n<<std::endl;
    if(ntTpc < 1)continue;
    
    for(int ntr = 0;ntr<ntTpc;ntr++){
      TPC_pid->Fill((*mom0)[ntr],(*dEdx)[ntr]);
      hist_dedx->Fill((*dEdx)[ntr]);
      if((*dEdx)[ntr]<30)mip_cut = true;
      for(int nhit = 0;nhit<(*nhtrack)[ntr];nhit++){
	hist_clu->Fill((*track_cluster_size)[ntr][nhit]);
	
	if((*track_cluster_size)[ntr][nhit] == 1)clu_cut = true;
	if((*track_cluster_y_center)[ntr][nhit] > -50 && (*track_cluster_y_center)[ntr][nhit] < 50)y_cut = true;
	
	int layerid = (*hitlayer)[ntr][nhit];
	int padid = tpc::GetPadId(layerid,(*track_cluster_row_center)[ntr][nhit]);
	double cnt = TPC_count->GetBinContent(padid+1);
	TPC_count->SetBinContent(padid+1,cnt+1);
	double de = (*track_cluster_de)[ntr][nhit];
	hist_de[0][padid]->Fill(de);
	if(clu_cut)hist_de[1][padid]->Fill(de);
	if(y_cut)hist_de[2][padid]->Fill(de);
	if(mip_cut)hist_de[3][padid]->Fill(de);
	if(clu_cut && y_cut)hist_de[4][padid]->Fill(de);
	if(clu_cut && mip_cut)hist_de[5][padid]->Fill(de);
	if(y_cut && mip_cut)hist_de[6][padid]->Fill(de);
	if(clu_cut && y_cut && mip_cut)hist_de[7][padid]->Fill(de);
      }
    }
  }

  
  TFile *f = new TFile("tpchelixtracking-padgain.root","RECREATE");
  
  gROOT->SetBatch(kTRUE);   
  TCanvas c("c","c",900,700);
  gStyle->SetOptStat(0);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("tpchelixtracking-padgain.cc");
  p->AddText("TPC Pad Gain Calibration");
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c.Print("tpchelixtracking-padgain.pdf(");
  c.Clear();
  gPad->SetLogz();
  TPC_pid->Draw("colz");
  c.Print("tpchelixtracking-padgain.pdf");

  c.Clear();
  hist_dedx->Draw();
  c.Print("tpchelixtracking-padgain.pdf");
  
  c.Clear();
  hist_clu->Draw();
  c.Print("tpchelixtracking-padgain.pdf");
  
  for(int ipad = 0;ipad <NumOfPadTPC;ipad++){
    TF1 fL("fL", "landau", 0, 1000);
    fL.SetParLimits(1,50,300);
    
    if(ipad%30 == 0){
      if(ipad !=0)c.Print("tpchelixtracking-padgain.pdf");
      c.Clear();
      c.Divide(6,5);
    }
    c.cd(ipad%30+1);
    
    hist_de[7][ipad]->Fit(&fL, "QR0");
    auto fit_param = fL.GetParameters();

    graph_gain->AddPoint(ipad,fit_param[1]);
    TPC_gain->SetBinContent(ipad+1,fit_param[1]);
    for(int i=0;i<cutflag;i++){
      hist_de[i][ipad]->SetLineColor(i+1);
      if(i==0)hist_de[i][ipad]->Draw();
      else{
	hist_de[i][ipad]->Draw("same");
      }
      hist_de[i][ipad]->Write();
    }
  }
  TPC_pid->Write();
  graph_gain->Write();
  
  f->Close();
  c.Print("tpchelixtracking-padgain.pdf");
  c.Clear();
  graph_gain->SetMarkerStyle(4);
  graph_gain->Draw("AP");
  c.Print("tpchelixtracking-padgain.pdf");
  c.Clear();
  TPC_gain->Draw("colz");
  c.Print("tpchelixtracking-padgain.pdf");
  c.Clear();
  TPC_count->Draw("colz");
  c.Print("tpchelixtracking-padgain.pdf)");   
   
}
