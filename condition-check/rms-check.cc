#include "../TPCPadHelper_260416.hh"

const int run[16]={2264,2408,2449,2487,2651,2686,2818,2836,2842,2893,2940,2979,3732,3745,3777,3835};
const int num_run = 16;
void rms_check(){

  double x[num_run];
  double ex[num_run];
  double eff[NumOfAsadTPC][num_run];
  double err[NumOfAsadTPC][num_run];
  
  auto c1 = new TCanvas("c1","c1");
    
  for(int nrun = 0;nrun<16;nrun++){
    int runnumber = run[nrun];
    x[nrun] = nrun+1;
    ex[nrun] = 0.0;
    TFile *file = TFile::Open(Form("result/tpc-condition-run0%d.root",runnumber));
    TH2D* hist_rms[NumOfAsadTPC];
    c1->Clear();
    TPad *titlePad = new TPad("titlePad","titlePad",0,0.94,1,1);
    titlePad->Draw();
    titlePad->cd();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextAlign(22);
    lat.SetTextSize(0.55);
    lat.DrawLatex(0.5,0.5,Form("run%05d",runnumber));
    c1->cd();
    TPad *plotPad_rms = new TPad("plotPad_rms","plotPad_rms",0,0,1,0.94);
    plotPad_rms->Draw();
    plotPad_rms->cd();
    plotPad_rms->Divide(8,4);
    
    for(int i=0;i<NumOfAsadTPC;i++){
      hist_rms[i] = (TH2D*)file->Get(Form("hist_rms%d",i));
      plotPad_rms->cd(i+1);
      hist_rms[i]->Draw("colz");
      TH1D *hy = hist_rms[i]->ProjectionY();
      int rms_cut = hy->FindBin(5.);
      double rms0 = hy->Integral(0,rms_cut-1);
      double rms_normal = hy->Integral(rms_cut+1,hy->GetNbinsX());
      double rms_total = rms0 + rms_normal;
      eff[i][nrun] = rms_normal / rms_total;
      err[i][nrun] = sqrt(eff[i][nrun]*(1-eff[i][nrun])/rms_total);
	
      std::cout<<"run0"<<runnumber<<", AsAd"<<i<<" Efficiency : "<<(double)(rms_normal) / (double)(rms0+rms_normal) * 100 <<"%"<<endl;
      
    }

    c1->SaveAs(Form("result/rms-run0%d.png",runnumber));
    //file->Close();
    delete plotPad_rms;
    delete titlePad;
    for(int i=0;i<NumOfAsadTPC;i++){
      delete hist_rms[i];
      
    }

    
  }
  TCanvas *c2 = new TCanvas("c2","c2");
  TFile *f = new TFile("result/rms-check.root","RECREATE");
  c2->Divide(8,4);
  for(int i=0;i<NumOfAsadTPC;i++){
    auto gr = new TGraphErrors(num_run,x,eff[i],ex,err[i]);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);
    gr->Draw("AP");
    TAxis *ax = gr->GetXaxis();
    for(int i=0;i<num_run;i++){
      ax->SetBinLabel(ax->FindBin(i+1), Form("%d",run[i]));
    }
    
    c2->cd(i+1);
    gr->SetMarkerStyle(20);
    gr->Draw("AP");
    gr->SetName(Form("g_eff_asad%d",i));
    gr->SetTitle(Form("g_eff_asad%d;run number;Efficiency",i));
    gr->Write();
    
  }

  f->Close();
    
  
}
