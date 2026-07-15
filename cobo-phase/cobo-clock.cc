#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

using namespace std;

const int NCobo = 8;
const int runnumber = 2599;
const bool param_update = false;
const bool auto_p1 = true;

Double_t ClockShiftFunc(Double_t *x, Double_t *par)
{
  return par[0] * TMath::Freq((x[0] - par[1]) / par[2])+par[3];
}

TH2D* CorrectHistY(TH2D *h, TF1 *fshift, const char *name)
{
  TH2D *hcorr = new TH2D(name, h->GetTitle(),
                         h->GetNbinsX(),
                         h->GetXaxis()->GetXmin(),
                         h->GetXaxis()->GetXmax(),
                         h->GetNbinsY(),
                         h->GetYaxis()->GetXmin(),
                         h->GetYaxis()->GetXmax());

  for(int ix = 1; ix <= h->GetNbinsX(); ix++){

    double x = h->GetXaxis()->GetBinCenter(ix);
    double shift = fshift->Eval(x);

    for(int iy = 1; iy <= h->GetNbinsY(); iy++){

      double y = h->GetYaxis()->GetBinCenter(iy);
      double cont = h->GetBinContent(ix, iy);

      if(cont <= 0) continue;

      double ycorr = y - shift;

      hcorr->Fill(x, ycorr, cont);
    }
  }

  return hcorr;
}

struct ClockJumpCandidate {
  double x;
  double dy;
  double score;
};

std::vector<ClockJumpCandidate> FindClockJumpCandidates(TH2D *h, double expectedJump, int icobo, TH1D *htrace)
{
  const int nx = h->GetNbinsX();
  const int ny = h->GetNbinsY();
  const int xGroup = 25;
  const int peakHalfWindow = 2;
  const int maxCandidates = 3;
  const double xSearchMin = -38.0;
  const double xSearchMax = 38.0;
  const double targetJump = (TMath::Abs(expectedJump) > 1.0) ? TMath::Abs(expectedJump) : 3.0;
  const double minJump = 1.5;
  const double maxJump = 6.0;

  struct TracePoint {
    double x;
    double y;
    double weight;
  };

  std::vector<TracePoint> trace;

  for(int ix0 = 1; ix0 <= nx; ix0 += xGroup){
    int ix1 = TMath::Min(nx, ix0 + xGroup - 1);
    double xsum = 0.0;
    std::vector<double> yproj(ny + 1, 0.0);

    for(int ix = ix0; ix <= ix1; ix++){
      xsum += h->GetXaxis()->GetBinCenter(ix);

      for(int iy = 1; iy <= ny; iy++){
        double cont = h->GetBinContent(ix, iy);
        if(cont <= 0) continue;
        yproj[iy] += cont;
      }
    }

    double xcenter = xsum / (ix1 - ix0 + 1);
    if(xcenter < xSearchMin || xcenter > xSearchMax) continue;

    int maxbin = 0;
    double maxCont = 0.0;
    for(int iy = 1; iy <= ny; iy++){
      if(yproj[iy] > maxCont){
        maxCont = yproj[iy];
        maxbin = iy;
      }
    }

    if(maxbin <= 0 || maxCont <= 0.0) continue;

    double sumw = 0.0;
    double sumwy = 0.0;
    for(int iy = TMath::Max(1, maxbin - peakHalfWindow);
        iy <= TMath::Min(ny, maxbin + peakHalfWindow); iy++){
      double cont = yproj[iy];
      if(cont <= 0) continue;

      sumw += cont;
      sumwy += cont * h->GetYaxis()->GetBinCenter(iy);
    }

    if(sumw <= 0.0) continue;

    double ycenter = sumwy / sumw;
    trace.push_back({xcenter, ycenter, sumw});
    htrace->SetBinContent(htrace->GetXaxis()->FindBin(xcenter), ycenter);
  }

  std::vector<ClockJumpCandidate> candidates;

  for(size_t i = 1; i < trace.size(); i++){
    double dy = trace[i].y - trace[i - 1].y;
    double absDy = TMath::Abs(dy);
    if(absDy < minJump || absDy > maxJump) continue;

    double score = absDy / (1.0 + TMath::Abs(absDy - targetJump));
    candidates.push_back({0.5 * (trace[i].x + trace[i - 1].x), dy, score});
  }

  auto preferNegativeP1 = [](int cobo){
    return cobo == 0 || cobo == 1 || cobo == 2 ||
           cobo == 5 || cobo == 6 || cobo == 7;
  };

  if(preferNegativeP1(icobo)){
    std::vector<ClockJumpCandidate> preferred;
    for(const auto &candidate : candidates){
      if(candidate.x > -20.0 && candidate.x < -10.0){
        preferred.push_back(candidate);
      }
    }

    if(!preferred.empty()){
      std::sort(preferred.begin(), preferred.end(),
                [](const ClockJumpCandidate &a, const ClockJumpCandidate &b){
                  double da = TMath::Abs(a.x + 14.5);
                  double db = TMath::Abs(b.x + 14.5);
                  if(TMath::Abs(da - db) > 1.e-6) return da < db;
                  return a.score > b.score;
                });

      if((int)preferred.size() > maxCandidates){
        preferred.resize(maxCandidates);
      }
      return preferred;
    }
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const ClockJumpCandidate &a, const ClockJumpCandidate &b){
              return a.score > b.score;
            });

  if((int)candidates.size() > maxCandidates){
    candidates.resize(maxCandidates);
  }

  return candidates;
}

void UpdateCoboParameter(const char* infile,
                         const char* outfile,
			 const double p0fit[],
                         const double p1fit[],
                         const int NCobo)
{
  std::ifstream fin(infile);
  std::ofstream fout(outfile);

  std::string line;
  bool inCobo = false;
  int idx = 0;

  while(std::getline(fin, line)){

    if(line.find("Coboparameter") != std::string::npos){
      inCobo = true;
      idx = 0;
      fout << line << "\n";
      continue;
    }

    if(inCobo && idx >= NCobo){
      inCobo = false;
    }

    if(inCobo){
      std::stringstream ss(line);

      int cobo, dummy, aty;
      double oldp0, oldp1, oldp2;

      if(ss >> cobo >> dummy >> aty >> oldp0 >> oldp1 >> oldp2){

        fout << std::setw(8)  << cobo
             << std::setw(8)  << dummy
             << std::setw(8)  << aty
             << std::setw(16) << std::fixed << std::setprecision(8) << -80
             << std::setw(16) << std::fixed << std::setprecision(8) << p1fit[idx]
             << std::setw(16) << std::fixed << std::setprecision(8) << oldp2
             << "\n";

        std::cout << "update cobo idx " << idx
                  << " : " << oldp0 << " " << oldp1
                  << " -> " << p0fit[idx] << " " << p1fit[idx]
                  << std::endl;

        idx++;
        continue;
      }
    }

    fout << line << "\n";
  }

  std::cout << "Total updated Cobo lines = " << idx << std::endl;
}


void cobo_clock(const char* result_subdir = "tpchit-test"){
  gROOT->SetBatch(kTRUE);
  const string result_dir = Form("result/%s", result_subdir);
  gSystem->mkdir(result_dir.c_str(), kTRUE);

  
  TFile *file = TFile::Open(Form("~/data/JPARC2025Nov_root/%s/run0%d_DstTPCHitBcOutTracking.root",result_subdir,runnumber),"read");
  TString histFmt = "TPCHit_ResY_vs_ClockTime_CoBo%d_RawClock";
  TString histFmt_cor = "TPCHit_ResY_vs_ClockTime_CoBo%d";

  string outpdf = Form("%s/cobo-clock-run0%d.pdf", result_dir.c_str(),runnumber);
  TFile *fout = new TFile(Form("%s/cobo-clock-run0%d.root",result_dir.c_str(),runnumber),"RECREATE");
  
  TCanvas *c1 = new TCanvas("c1","c1");
  
  gStyle->SetOptStat(0);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("cobo-clock.cc");
  p->AddText("CoBo Clocktime Phaseshift Correction");
  p->AddText(Form("run0%d",runnumber));
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());
  
  c1->Clear();
  c1->Divide(4,2);
  
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(4,2);
  TF1 *fshift[NCobo];
  TH2D *h2[NCobo];
  TH2D *h2_cor[NCobo];
  double p1[NCobo]={0.};
  double p0_fit[NCobo];
  double p1_fit[NCobo];
  for (int icobo = 0; icobo < NCobo; ++icobo) {
    cout<<icobo<<endl;
    
    TString hname = Form(histFmt, icobo);
    h2[icobo] = (TH2D*)file->Get(hname);
    h2[icobo]->RebinY(2);

    int bin1 = h2[icobo]->GetXaxis()->FindBin(-40.0);
    int bin2 = h2[icobo]->GetXaxis()->FindBin(-35.0);
    int bin3 = h2[icobo]->GetXaxis()->FindBin(30.0);
    int bin4 = h2[icobo]->GetXaxis()->FindBin(40.0);
    
    TH1D *hleft = h2[icobo]->ProjectionY(Form("hleft%d",icobo),bin1,bin2);
    TH1D *hright = h2[icobo]->ProjectionY(Form("hright%d",icobo),bin3,bin4);

    int maxbin_left = hleft->GetMaximumBin();
    double max_left = hleft->GetBinCenter(maxbin_left);

    int maxbin_right = hright->GetMaximumBin();
    double max_right = hright->GetBinCenter(maxbin_right);
    
    TF1 *fleft = new TF1(Form("fleft%d",icobo),"gaus",max_left-5,max_left+5);
    TF1 *fright = new TF1(Form("fright%d",icobo),"gaus",max_right-5,max_right+5);

    hleft->Fit(fleft,"RQ0");
    hright->Fit(fright,"RQ0");
    double yleft = fleft->GetParameter(1);
    double yright = fright->GetParameter(1);
    double p0 = yright - yleft;

    hleft->Write();
    hright->Write();
    fshift[icobo] = new TF1(Form("fshift%d",icobo),ClockShiftFunc,-60,60,4);

    double p2 = 0.01;
    double p3 = yleft;
    TH1D *hmean = new TH1D(Form("hmean_cobo%d",icobo),
                           Form("CoBo%d mean residual vs raw clock;Raw clock;Mean residual Y",icobo),
                           h2[icobo]->GetNbinsX(),
                           h2[icobo]->GetXaxis()->GetXmin(),
                           h2[icobo]->GetXaxis()->GetXmax());

    double p1_init = p1[icobo];
    double found_jump = 0.0;
    std::vector<ClockJumpCandidate> candidates;
    if(auto_p1){
      candidates = FindClockJumpCandidates(h2[icobo], p0, icobo, hmean);
    }

    if(!candidates.empty()){
      p1_init = candidates[0].x;
      found_jump = candidates[0].dy;
      cout << "auto p1 cobo " << icobo << " = " << p1_init
           << "  grouped jump = " << found_jump << "  expected = " << p0 << endl;
      for(size_t icand = 0; icand < candidates.size(); icand++){
        cout << "  candidate " << icand
             << " p1 = " << candidates[icand].x
             << " jump = " << candidates[icand].dy
             << " score = " << candidates[icand].score << endl;
      }
    }
    else{
      cout << "auto p1 failed cobo " << icobo
           << ", use fallback p1 = " << p1_init << endl;
    }

    fshift[icobo]->SetParameters(p0,p1_init,p2,p3);
    
    fshift[icobo]->SetParLimits(0,p0-0.1,p0+0.1);
    fshift[icobo]->SetParLimits(1,p1_init-2.0,p1_init+2.0);
    //fshift[icobo]->SetParLimits(2,0,p2+0.02);
    //fshift[icobo]->FixParameter(0,p0);
    fshift[icobo]->FixParameter(2,p2);
    fshift[icobo]->SetParLimits(3,p3-0.2,p3+0.2);
    
    c1->cd(icobo+1);
    h2[icobo]->Fit(fshift[icobo],"R");
    h2[icobo]->SetTitle(Form("CoBo%d",icobo));
    h2[icobo]->Draw("colz");
    
    hmean->Write();
    h2[icobo]->Write();
    fshift[icobo]->Write();
    p0_fit[icobo] = fshift[icobo]->GetParameter(0) * -1 / 0.055;
    p1_fit[icobo] = fshift[icobo]->GetParameter(1);
    
  }
  //c1->Modified();
  //c1->Update();
  c1->Print(outpdf.c_str());
  c1->Clear();
  c1->Divide(4,2);
  for(int icobo=0;icobo<NCobo;icobo++){
    
    TH2D *hcorr = CorrectHistY(h2[icobo], fshift[icobo], Form("hcorr_cobo%d",icobo));
    c1->cd(icobo+1);
    hcorr->Draw("colz");
    hcorr->Write();

  }

  
  

  //c1->Modified();
  //c1->Update();
  c1->Print(outpdf.c_str());

  c1->Clear();
  c1->Divide(4,2);
  for(int icobo=0;icobo<NCobo;icobo++){

    TString hname = Form(histFmt_cor, icobo);
    h2_cor[icobo] = (TH2D*)file->Get(hname);
    h2_cor[icobo]->RebinY(2);
    c1->cd(icobo+1);
    h2_cor[icobo]->Draw("colz");
    h2_cor[icobo]->Write();
    
  }
  //c1->Modified();
  //c1->Update();
  c1->Print((outpdf + ")").c_str());


  
  fout->Close();

  if(param_update){
    UpdateCoboParameter("param_history/TPCParam_e72_run02447_lasthit_2",
			Form("param_history/TPCParam_e72_run0%d_lasthit",runnumber),
			p0_fit, p1_fit, NCobo);
  }
}
