#include "../TPCPadHelper.hh"
#include "../TPCEventDisplayHelper.hh"

const int runnumber[] = {2447, 2449, 2450, 2451, 2452, 2453, 2454, 2456, 2457, 2458, 2459, 2460, 2462, 2463, 2465,2466, 2468};
const int nrun = sizeof(runnumber) / sizeof(runnumber[0]);

struct FitResult {
  TF1  *f = nullptr;
  TH1D *h = nullptr;
};

bool IsGoodLandauFit(TH1D* h, TF1* f)
{
  if (!h || !f) return false;
  if (h->GetEntries() < 30) return false;
  if (f->GetNDF() <= 0) return false;

  const double mpv = f->GetParameter(1);
  if (mpv < 50. || mpv > 350.) return false;

  return true;
}

FitResult FitLandauWithRetry(TH1D* hraw)
{
  if (!hraw) return {};

  std::vector<int> rebin_list = {1, 2, 4};
  for (auto rebin : rebin_list) {
    TH1D* h = (TH1D*)hraw->Clone(Form("%s_tmp_rebin%d", hraw->GetName(), rebin));
    h->SetDirectory(nullptr);
    h->Rebin(rebin);

    TF1* f = new TF1(Form("f_%s_rebin%d", hraw->GetName(), rebin), "landau", 100, 700);
    f->SetParameter(1, 200);
    f->SetParameter(2, 33);
    f->SetLineColor(kRed);

    h->Fit(f, "RQ0N");
    if (IsGoodLandauFit(h, f)) return {f, h};

    delete f;
    delete h;
  }

  return {};
}


void padgain_combine(const char* result_subdir = "260701")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const string result_dir = Form("result/%s", result_subdir);
  gSystem->mkdir(result_dir.c_str(), kTRUE);

  const string outpdf = Form("%s/padgain-combine.pdf", result_dir.c_str());
  const string outroot = Form("%s/padgain-combine.root", result_dir.c_str());

  TH1D* hist_de[NumOfPadTPC] = {nullptr};
  auto TPC_count = tpcdisp::MakeTPCPadMap("TPC_count", "TPC_count");
  auto TPC_gain = tpcdisp::MakeTPCPadMap("TPC_gain", "TPC_gain");
  

  auto graph_gain = new TGraph();
  auto graph_chi = new TGraph();
  auto graph_sigma = new TGraph();
  graph_gain->SetName("graph_gain");
  graph_chi->SetName("graph_chi");
  graph_sigma->SetName("graph_sigma");

  auto hist_chi = new TH1D("hist_chi", "hist_chi;#chi^2/NDF;Counts", 500, 0, 10);
  auto hist_sigma = new TH1D("hist_sigma", "hist_sigma;#sigma;Counts", 100, 0, 200);

  for (int irun = 0; irun < nrun; ++irun) {
    const string input_root = Form("%s/physics-padgain-run0%d.root", result_dir.c_str(), runnumber[irun]);
    TFile* file = TFile::Open(input_root.c_str(), "READ");
    if (!file || file->IsZombie()) {
      cout << "Cannot open " << input_root << endl;
      continue;
    }

    

    for (int ipad = 0; ipad < NumOfPadTPC; ++ipad) {
      //const double run_cnt = TPC_tr_cluster->GetBinContent(ipad + 1);
      //const double total_cnt = TPC_count->GetBinContent(ipad + 1);
      

      auto hist = (TH1D*)file->Get(Form("hist_de%d", ipad));
      if (!hist) continue;

      if (!hist_de[ipad]) {
        hist_de[ipad] = (TH1D*)hist->Clone(Form("hist_de%d", ipad));
        hist_de[ipad]->SetDirectory(nullptr);
      }
      else {
        hist_de[ipad]->Add(hist);
      }
    }
    
    file->Close();
  }
  for(int ipad = 0; ipad < NumOfPadTPC;++ipad){
    double total_entry = hist_de[ipad] -> Integral();
    TPC_count->SetBinContent(ipad + 1, total_entry);
    
  }
  

  TFile* fout = new TFile(outroot.c_str(), "RECREATE");
  auto c1 = new TCanvas("c1", "c1");

  TPaveText* title = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  title->AddText("padgain-combine.cc");
  title->AddText("TPC Pad Gain Calibration");
  title->AddText(Form("result directory: %s", result_dir.c_str()));
  TDatime now;
  title->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",
                      now.GetYear(), now.GetMonth(), now.GetDay(),
                      now.GetHour(), now.GetMinute(), now.GetSecond()));
  title->Draw();
  c1->Print((outpdf + "(").c_str());

  for (int ipad = 0; ipad < NumOfPadTPC; ++ipad) {
    if (ipad % 30 == 0) {
      if (ipad != 0) c1->Print(outpdf.c_str());
      c1->Clear();
      c1->Divide(6, 5);
    }

    c1->cd(ipad % 30 + 1);
    FitResult res = FitLandauWithRetry(hist_de[ipad]);
    if (res.f && res.h) {
      const double mpv = res.f->GetParameter(1);
      const double sigma = res.f->GetParameter(2);
      const double chi2ndf = res.f->GetChisquare() / res.f->GetNDF();

      hist_chi->Fill(chi2ndf);
      hist_sigma->Fill(sigma);
      graph_gain->AddPoint(ipad, mpv);
      graph_chi->AddPoint(ipad, chi2ndf);
      graph_sigma->AddPoint(ipad, sigma);
      TPC_gain->SetBinContent(ipad + 1, mpv);

      res.h->Draw();
      res.f->Draw("same");
      res.h->Write(Form("hist_de%d", ipad));
    }
    else if (hist_de[ipad]) {
      hist_de[ipad]->Draw();
      hist_de[ipad]->Write();
    }
  }
  c1->Print(outpdf.c_str());

  c1->Clear();
  TPC_gain->SetMinimum(0);
  TPC_gain->SetMaximum(400);
  TPC_gain->Draw("colz");
  TPC_gain->Write();
  c1->Print(outpdf.c_str());

    

  c1->Clear();
  auto mg = new TMultiGraph();
  graph_gain->SetMarkerStyle(4);
  graph_chi->SetMarkerStyle(4);
  graph_chi->SetMarkerColor(kRed);
  graph_sigma->SetMarkerStyle(4);
  graph_sigma->SetMarkerColor(kBlue);
  mg->Add(graph_gain);
  mg->Add(graph_chi);
  mg->Add(graph_sigma);
  mg->SetMinimum(0);
  mg->SetMaximum(400);
  mg->Draw("AP");
  graph_gain->Write();
  graph_chi->Write();
  graph_sigma->Write();
  c1->Print(outpdf.c_str());

  c1->Clear();
  gPad->SetLogz();
  TPC_count->SetMinimum(1);
  TPC_count->Draw("colz");
  TPC_count->Write();
  c1->Print((outpdf + ")").c_str());

  hist_chi->Write();
  hist_sigma->Write();
  fout->Close();
}
