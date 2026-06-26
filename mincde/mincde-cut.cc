const int mincde[] = {30,50,70,90,110};
const int nmincde = sizeof(mincde) / sizeof(mincde[0]);
const int runnumber = 2599;

void mincde_cut(){
  gROOT->SetBatch(kTRUE);
  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/mincde";
  string outpdf = Form("result/mincde-cut-run0%d.pdf", runnumber);

  TH1D *hist_chi[nmincde];
  TH1D *hist_resx[nmincde];
  TH1D *hist_resy[nmincde];
  TH1D *hist_resz[nmincde];
  TH1D *hist_nhit[nmincde];
  TH1D *hist_rawde[nmincde];
  TH1D *hist_trackde[nmincde];
  for(int i=0;i<nmincde;i++){
    hist_chi[i] = new TH1D(Form("hist_chi%d",mincde[i]),Form("hist_chi%d",mincde[i]),100,0,7);
    hist_resx[i] = new TH1D(Form("hist_resx%d",mincde[i]),Form("hist_resx%d",mincde[i]),200,-5,5);
    hist_resy[i] = new TH1D(Form("hist_resy%d",mincde[i]),Form("hist_resy%d",mincde[i]),200,-5,5);
    hist_resz[i] = new TH1D(Form("hist_resz%d",mincde[i]),Form("hist_resz%d",mincde[i]),200,-5,5);
    hist_nhit[i] = new TH1D(Form("hist_nhit%d",mincde[i]),Form("hist_nhit%d",mincde[i]),50,0,50);
    hist_rawde[i] = new TH1D(Form("hist_rawde%d",mincde[i]),Form("hist_rawde%d",mincde[i]),50,0,50);
    hist_trackde[i] = new TH1D(Form("hist_trackde%d",mincde[i]),Form("hist_trackde%d",mincde[i]),50,0,50);
  }

  auto c1 = new TCanvas("c1","c1");
  gStyle->SetOptStat(0);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("mincde-cut.cc");
  p->AddText("TPC MinCDe Check");
  p->AddText(Form("run0%d",runnumber));
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());
  c1->Clear();

  c1->Divide(3,2);
  TLegend *leg = new TLegend(0.65,0.7,0.88,0.88);
   
  for(int i=0;i<nmincde;i++){
    TFile *file = new TFile(Form("%s/run0%d_DstTPCCalibration_mincde_%d.root",dir.c_str(),runnumber,mincde[i]));
    TTree *tree = (TTree*)file->Get("tpc");
    int ntTpc;
    vector<double>* chisqr = nullptr;
    vector<vector<double>>* closeDistTpc = nullptr;
    vector<vector<double>>* residual_x = nullptr;
    vector<vector<double>>* residual_y = nullptr;
    vector<vector<double>>* residual_z = nullptr;
    
    vector<int>* nhrawtrack = nullptr;
    vector<vector<double>>* track_hit_pad = nullptr;
    vector<vector<double>>* track_hit_de = nullptr;
    
    tree->SetBranchAddress("ntTpc",&ntTpc);
    tree->SetBranchAddress("chisqr",&chisqr);
    tree->SetBranchAddress("closeDistTpc",&closeDistTpc);
    tree->SetBranchAddress("residual_x",&residual_x);
    tree->SetBranchAddress("residual_y",&residual_y);
    tree->SetBranchAddress("residual_z",&residual_z);
    tree->SetBranchAddress("nhrawtrack",&nhrawtrack);
    tree->SetBranchAddress("track_hit_pad",&track_hit_pad);
    tree->SetBranchAddress("track_hit_de",&track_hit_de);

    for(int n=0;n<tree->GetEntries();n++){
      tree->GetEntry(n);
      for(int nt = 0;nt<ntTpc;nt++){
	hist_chi[i]->Fill((*chisqr)[nt]);

	for(int ncom =0;ncom<(*residual_x)[nt].size();ncom++){
	  hist_resx[i]->Fill((*residual_x)[nt][ncom]);
	}	
	for(int ncom =0;ncom<(*residual_y)[nt].size();ncom++){
	  hist_resy[i]->Fill((*residual_y)[nt][ncom]);
	}
	for(int ncom =0;ncom<(*residual_z)[nt].size();ncom++){
	  hist_resz[i]->Fill((*residual_z)[nt][ncom]);
	}

	hist_nhit[i]->Fill((*nhrawtrack)[nt]);
      }
    }
    c1->cd(1);
    hist_chi[i]->SetLineColor(i+1);
    hist_chi[i]->Draw("same");
    c1->cd(2);
    hist_resx[i]->SetLineColor(i+1);
    hist_resx[i]->Draw("same");
    c1->cd(3);
    hist_resy[i]->SetLineColor(i+1);
    hist_resy[i]->Draw("same");
    c1->cd(4);
    hist_resz[i]->SetLineColor(i+1);
    hist_resz[i]->Draw("same");
    c1->cd(5);
    hist_nhit[i]->SetLineColor(i+1);
    hist_nhit[i]->Draw("same");
    
    leg->AddEntry(hist_nhit[i],Form("mincde : %d",mincde[i]),"l");
  }
  c1->cd(6);
  leg->Draw();

  c1->Print((outpdf + ")").c_str());
}
