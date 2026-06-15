#include "../TPCPadHelper_260416.hh"

const int runnum[2] = {2856,3849}; //933 Beam Through
//const int runnum[2] = {2370,3777}; //735 Physics
const double tTpc_cut[2] = {58.,80.};
//const double tTpc_cut[2] = {30.,140.};
const double deTpc_cut[2] = {100.,1000.};

void compare_cluster(){
  //933 Beam Through

  string dir = "/hsm/had/sks/E72/JPARC2025Nov/rootfile/tpchit/v2";
  string pdf_name = "cluster_run"+to_string(runnum[0])+"_run"+to_string(runnum[1])+"_tight.pdf";
  
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Print(Form("%s(",pdf_name.c_str()));
  c1->Divide(2);
  
  int nrun = sizeof(runnum)/sizeof(runnum[0]);
  TH1D *hist_cl[nrun][NumOfLayersTPC];
  for(int i=0;i<nrun;i++){
    for(int j=0;j<NumOfLayersTPC;j++){
      hist_cl[i][j] = new TH1D(Form("hist_cl%d%d",i,j),Form("hist_cl%d%d",i,j),20,0,20);
    }
  }
  for(int i=0;i<nrun;i++){

    auto TPC_count = new TH2Poly("TPC_count", "TPC_count;Z;X", MinZ, MaxZ, MinX, MaxX);
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
      }
    }
    TFile *file = new TFile(Form("%s/run0%d_TPCHit.root",dir.c_str(),runnum[i]));
    TTree *tree = (TTree*)file->Get("tpc");

    int nhTpc;
    vector<int> *layerTpc = nullptr;
    vector<int> *rowTpc = nullptr;
    vector<int> *padTpc = nullptr;
    vector<double> *tTpc = nullptr;
    vector<double> *deTpc = nullptr;

    tree->SetBranchAddress("nhTpc",&nhTpc);
    tree->SetBranchAddress("layerTpc",&layerTpc);
    tree->SetBranchAddress("rowTpc",&rowTpc);
    tree->SetBranchAddress("padTpc",&padTpc);
    tree->SetBranchAddress("tTpc",&tTpc);
    tree->SetBranchAddress("deTpc",&deTpc);
    
    bool beam_window_cut = false;
    bool de_cut = false;
    for(int n=0;n<tree->GetEntries();n++){
      tree->GetEntry(n);
      beam_window_cut = false;
      de_cut = false;
      int layer_cnt[NumOfLayersTPC] = {0};
      if(nhTpc == 0)continue;
      
      for(int nhit = 0; nhit < nhTpc; nhit++){
	if(tpc::IsDead((*layerTpc)[nhit],(*rowTpc)[nhit]) || tpc::Noise((*padTpc)[nhit]))continue;
	if((*tTpc)[nhit] > tTpc_cut[0] && (*tTpc)[nhit] < tTpc_cut[1])beam_window_cut = true;
	if((*deTpc)[nhit] > deTpc_cut[0] && (*deTpc)[nhit] < deTpc_cut[1])de_cut = true;
	if(beam_window_cut && de_cut){
	  //if(beam_window_cut){
	  layer_cnt[(*layerTpc)[nhit]]++;
	  double cont = TPC_count->GetBinContent((*padTpc)[nhit]+1);
	  TPC_count->SetBinContent((*padTpc)[nhit]+1,cont+1);
	}
      } //nhit
      for(int nlayer = 0;nlayer < NumOfLayersTPC;nlayer++){
	hist_cl[i][nlayer]->Fill(layer_cnt[nlayer]);
      } //nlayer
    }//n
    c1->cd(i+1);
    TPC_count->Draw("colz");
    
  }//run

  c1->Print(Form("%s",pdf_name.c_str()));
  
  auto mg = new TMultiGraph();
  TLegend *leg = new TLegend(0.15, 0.7, 0.4, 0.88);
  c1->Clear();
  c1->Divide(8,4);
  for(int n=0;n<nrun;n++){
    //c1->Clear();
    TGraph *g_cl = new TGraph();
    for(int i=0;i<NumOfLayersTPC;i++){
      c1->cd(i+1);
      hist_cl[n][i]->Scale(1.0 / hist_cl[n][i]->Integral());
      if(n==0)hist_cl[n][i]->SetLineColor(kBlack);
      else if(n==1)hist_cl[n][i]->SetLineColor(kRed);
      hist_cl[n][i]->GetYaxis()->SetRangeUser(0,0.4);
      hist_cl[n][i]->Draw("HIST same");
      double mean = hist_cl[n][i]->GetMean();
      double var  = hist_cl[n][i]->GetRMS()*hist_cl[n][i]->GetRMS();
      cout<<"mean : "<<mean<<" , var : "<<var<<endl;
      g_cl->SetPoint(i,i,mean);
    }

    c1->Print(Form("%s",pdf_name.c_str()));
    g_cl->SetMarkerStyle(8);
    //g_cl->SetMarkerSize(2);
    if(n==1)g_cl->SetMarkerColor(kRed);
    mg->Add(g_cl);
    
    if(n==0)
      leg->AddEntry(g_cl,"2025Nov","p");
    else if(n==1)
      leg->AddEntry(g_cl,"2026April","p");
  }
  c1->Clear();
  mg->SetTitle(";Layer;Cluster Mean");
  mg->Draw("AP");
  leg->Draw("sames");
  c1->Print(Form("%s)",pdf_name.c_str()));
  
}
    
