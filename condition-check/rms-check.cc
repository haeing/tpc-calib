#include "../TPCPadHelper_260416.hh"

const int run[16]={2264,2408,2449,2487,2651,2686,2818,2836,2842,2893,2940,2979,3732,3745,3777,3835};

void rms_check(){
  
  auto c1 = new TCanvas("c1","c1");
  string outpdf_rms = "result/tpc-condition-rms.pdf";
    TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
    p->AddText("TPC RMS check per AsAd");
    TDatime now;
    p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
    p->Draw();
    c1->Print((outpdf_rms + "(").c_str());
    
  for(int nrun = 0;nrun<16;nrun++){
    int runnumber = run[nrun];
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
    }

    
    c1->Print(outpdf_rms.c_str());
    for(int i=0;i<NumOfAsadTPC;i++){
      hist_rms[i]->SetDirectory(0);
      delete hist_rms[i];
      
    }

    
  }

  c1->Clear();
  c1->Print((outpdf_rms+")").c_str());
  
  
  
}
