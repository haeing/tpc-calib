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
  tree->SetBranchAddress("nhtrack",&nhtrack);
  tree->SetBranchAddress("track_hit_pad",&track_hit_pad);
  tree->SetBranchAddress("track_hit_de",&track_hit_de);

  TH1D *hist_de[NumOfPadTPC];
  TH1D *hist_track_de[NumOfPadTPC];

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
	hist_track_de[(int)((*track_hit_pad)[i][j])]->Fill((*track_hit_de)[i][j]);
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
  f->Close();
}
