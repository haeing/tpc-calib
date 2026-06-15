#include "../TPCPadHelper_260416.hh"

const int start_run = 2264;
//const int end_run = 3043;
const int end_run = 2264;

const int evt_ave = 100;

void tpc_condition(){
  gROOT->SetBatch(kTRUE);
  //string dir = "/hsm/had/sks/E72/JPARC2025Nov/rootfile/tpchit/v2";
  string dir = "/home/had/haein/data/JPARC2025Nov_root/rms";

  string outpdf_hit = Form("result/tpc-condition-hit-run0%d-run0%d.pdf",start_run,end_run);
  string outpdf_rms = Form("result/tpc-condition-rms-run0%d-run0%d.pdf",start_run,end_run);
  auto c1 = new TCanvas("c1","c1");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("tpc-condition.cc");
  p->AddText(Form("run0%d - 0%d",start_run,end_run));
  p->AddText("TPC Hit check per AsAd");
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf_hit + "(").c_str());
  

  auto c2 = new TCanvas("c2","c2");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TPaveText *p1 = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p1->AddText("tpc-condition.cc");
  p1->AddText(Form("run0%d - 0%d",start_run,end_run));
  p1->AddText("TPC RMS check per AsAd");
  p1->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p1->Draw();
  c2->Print((outpdf_rms + "(").c_str());
    

  int current_run;
  for(int run = start_run;run<=end_run;run++){
    TString fname = Form("%s/run0%d_TPCRaw.root",dir.c_str(), run);
    if(gSystem->AccessPathName(fname)){
      cout << "No file : " << fname << endl;
      continue;
    }
    current_run = run;
    TFile *file = TFile::Open(fname);
    TTree *tree = (TTree*)file->Get("tpc");

    const int total_evt = tree->GetEntries();
    const int n_evt = (total_evt + evt_ave - 1) / evt_ave;
    

    int nhTpc;
    vector<int>* padTpc = nullptr;
    vector<double>* rawrmsTpc = nullptr;

    tree->SetBranchAddress("nhrawTpc",&nhTpc);
    tree->SetBranchAddress("padTpc",&padTpc);
    tree->SetBranchAddress("rawrmsTpc",&rawrmsTpc);

    TFile *f = new TFile(Form("result/tpc-condition-run0%d.root",run),"RECREATE");
    
    c1->Clear();
    c2->Clear();
    TPad *titlePad = new TPad("titlePad","titlePad",0,0.94,1,1);
    titlePad->Draw();
    titlePad->cd();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextAlign(22);
    lat.SetTextSize(0.55);
    lat.DrawLatex(0.5,0.5,Form("run%05d",run));
    c1->cd();
    TPad *plotPad_hit = new TPad("plotPad_hit","plotPad_hit",0,0,1,0.94);
    plotPad_hit->Draw();
    plotPad_hit->cd();
    plotPad_hit->Divide(8,4);

    
    c2->cd();
    titlePad->Draw();
    titlePad->cd();
    lat.DrawLatex(0.5,0.5,Form("run%05d",run));
    c2->cd();
    TPad *plotPad_rms = new TPad("plotPad_rms","plotPad_rms",0,0,1,0.94);
    plotPad_rms->Draw();
    plotPad_rms->cd();
    plotPad_rms->Divide(8,4);

    
    TH1D* hist_hitnum[NumOfAsadTPC];
    TH2D* hist_rms[NumOfAsadTPC];
    for(int i=0;i<NumOfAsadTPC;i++){
      hist_hitnum[i] = new TH1D(Form("hist_hitnum%d",i),Form("AsAd%d;Event Id;Counts",i),n_evt,0,total_evt);
      hist_rms[i] = new TH2D(Form("hist_rms%d",i),Form("AsAd%d;Event Id;RMS",i),n_evt,0,total_evt,100,0,100);
    }

    for(int n=0;n<total_evt;n++){
      tree->GetEntry(n);
      int ibin = n / evt_ave + 1;
      
      for(int nh = 0;nh<nhTpc;nh++){
	int layer = tpc::getLayerID((*padTpc)[nh]);
	int row = tpc::getRowID((*padTpc)[nh]);
	int iasad = tpc::GetASADId(layer,row);
	
	double val = hist_hitnum[iasad]->GetBinContent(ibin);
	hist_hitnum[iasad]->SetBinContent(ibin,val+1.);
	hist_rms[iasad]->Fill(n,(*rawrmsTpc)[nh]);
      }
    }
    for(int i=0;i<NumOfAsadTPC;i++){
      hist_hitnum[i]->Write();
      hist_rms[i]->Write();
      //c1->cd(i+1);
      plotPad_hit->cd(i+1);
      hist_hitnum[i]->Draw("hist");
      plotPad_rms->cd(i+1);
      hist_rms[i]->Draw("colz");
    }
    if(run == end_run)c1->Print((outpdf_hit + ")").c_str());
    else{c1->Print(outpdf_hit.c_str());}

    if(run == end_run)c2->Print((outpdf_rms + ")").c_str());
    else{c2->Print(outpdf_rms.c_str());}
    f->Close();
  }
  if(current_run !=end_run){
    c1->Print((outpdf_hit + ")").c_str());
    c2->Print((outpdf_rms + ")").c_str());
  }
  
  
  
}
