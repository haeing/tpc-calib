#include "../TPCPadHelper_260416.hh"

void htofcalib_mincde(int runnumber){
  gROOT->SetBatch(kTRUE);
  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/mincde";
  TFile *file = new TFile(Form("%s/run0%d_DstTPCCalibration.root",dir.c_str(),runnumber));
  TTree *tree = (TTree*)file->Get("tpc");
  string outpdf = Form("result/htofcalib-mincde-run0%d.pdf", runnumber);

  int nhTpc;
  vector<int>* raw_padid = nullptr;
  vector<double>* raw_de = nullptr;
  

  int ntTpc;
  vector<int>* nhtrack = nullptr;
  vector<vector<double>>* track_hit_pad = nullptr;
  vector<vector<double>>* track_hit_de = nullptr;

  tree->SetBranchAddress("nhTpc",&nhTpc);
  tree->SetBranchAddress("raw_padid",&raw_padid);
  tree->SetBranchAddress("raw_de",&raw_de);
  
  tree->SetBranchAddress("ntTpc",&ntTpc);
  tree->SetBranchAddress("nhrawtrack",&nhtrack);
  tree->SetBranchAddress("track_hit_pad",&track_hit_pad);
  tree->SetBranchAddress("track_hit_de",&track_hit_de);

  TH1D *hist_de[NumOfPadTPC];
  TH1D *hist_track_de[NumOfPadTPC];

  auto TPC_hit = new TH2Poly("TPC_hit","TPC_hit;Z [mm];X [mm]",MinZ,MaxZ,MinX,MaxX);
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
      TPC_hit->AddBin(5, X, Y);
    }
  }
  
  
  for(int i=0;i<NumOfPadTPC;i++){
    hist_de[i] = new TH1D(Form("hist_de%d",i),Form("hist_de%d",i),100,0,1000);
    hist_track_de[i] = new TH1D(Form("hist_track_de%d",i),Form("hist_track_de%d",i),100,0,1000);
  }

  for(int n = 0;n<tree->GetEntries();n++){
    tree->GetEntry(n);
    if(n%1000 == 0)cout<<n<<endl;
    for(int i=0;i<nhTpc;i++){
      hist_de[(*raw_padid)[i]]->Fill((*raw_de)[i]);
    }
    for(int i=0;i<ntTpc;i++){
      for(int j=0;j<(*nhtrack)[i];j++){
	int padid = (int)((*track_hit_pad)[i][j]);
	hist_track_de[padid]->Fill((*track_hit_de)[i][j]);
	if((*track_hit_de)[i][j]>90){
	  double cnt = TPC_hit->GetBinContent(padid+1);
	  TPC_hit->SetBinContent(padid+1,cnt+1);
	  
	}
      }
    }
  }
  auto c1 = new TCanvas("c1","c1");
  gStyle->SetOptStat(0);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("htofcalib-mincde.cc");
  p->AddText("TPC MinCDe Check");
  p->AddText(Form("run0%d",runnumber));
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());

  c1->Clear();
  TFile *f = new TFile(Form("result/htofcalib-mincde-run0%d.root",runnumber),"RECREATE");

  for(int ipad = 0;ipad <NumOfPadTPC;ipad++){
       
    if(ipad%30 == 0){
      if(ipad !=0)c1->Print(outpdf.c_str());
      c1->Clear();
      c1->Divide(6,5);
    }
    c1->cd(ipad%30+1);
    

    hist_de[ipad]->SetLineColor(kBlack);
    hist_de[ipad]->Draw();
    hist_track_de[ipad]->SetLineColor(kRed);
    hist_track_de[ipad]->Draw("same");
    
    hist_de[ipad]->Write();
    hist_track_de[ipad]->Write();
  }
  c1->Print((outpdf + ")").c_str());
  TPC_hit->Write();
  f->Close();
}
