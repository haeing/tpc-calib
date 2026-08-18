#include <iostream>
#include <memory>
#include <map>
#include <string>
#include <cmath>
#include <vector>

#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TH2.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

const std::vector<int> runnumbers = {2599, 2601, 2602, 2603, 2604};  // Add run numbers here.
namespace {
constexpr int kNumDriftBins = 7;
constexpr double kDriftBinWidthCm = 10.;

double GeneralizedGaussianRms(double alpha, double beta)
{
  if (alpha <= 0. || beta <= 0.) return 0.;
  return alpha * std::sqrt(TMath::Gamma(3. / beta) / TMath::Gamma(1. / beta));
}
}
void transverse_diffusion()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleH(0.055);
  const char* input_dir ="/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/calibration";

  if (runnumbers.empty()) {
    std::cerr << "No run number is configured." << std::endl;
    return;
  }
  std::vector<std::unique_ptr<TFile>> input_files;
  TString run_label;
  for (const int run : runnumbers) {
    const std::string input_name = Form("%s/run0%d_DstTPCHelixTracking.root", input_dir, run);
    std::unique_ptr<TFile> input(TFile::Open(input_name.c_str(), "READ"));
    if (!input || input->IsZombie()) {
      std::cerr << "Cannot open " << input_name << std::endl;
      continue;
    }
    input_files.push_back(std::move(input));
    if (!run_label.IsNull()) run_label += ", ";
    run_label += Form("run0%d", run);
  }
  if (input_files.empty()) return;

  std::map<std::string, std::unique_ptr<TH2>> summed_histograms;
  const auto get_hist = [&](const TString& name) -> TH2* {
    const auto cached = summed_histograms.find(name.Data());
    if (cached != summed_histograms.end()) return cached->second.get();
    std::unique_ptr<TH2> total;
    for (const auto& input : input_files) {
      auto* source = dynamic_cast<TH2*>(input->Get(name));
      if (!source) continue;
      if (!total) {
        total.reset(dynamic_cast<TH2*>(source->Clone(Form("sum_%s", name.Data()))));
        total->SetDirectory(nullptr);
      } else {
        total->Add(source);
      }
    }
    if (!total) return nullptr;
    auto* result = total.get();
    summed_histograms.emplace(name.Data(), std::move(total));
    return result;
  };

  const std::string output_dir = "result";
  gSystem->mkdir(output_dir.c_str(), kTRUE);
  const std::string output_name =
    Form("%s/transverse-diffusion-%s.pdf", output_dir.c_str(),
         runnumbers.size() == 1 ? Form("run0%d", runnumbers.front()) : "combined-runs");

  auto* canvas = new TCanvas("c_transverse_diffusion",
                             "TPC transverse diffusion", 1000, 800);
  canvas->Print((output_name + "[").c_str());

  TPaveText cover(0.10, 0.10, 0.90, 0.90, "NDC");
  cover.SetBorderSize(0);
  cover.SetFillStyle(0);
  cover.SetTextAlign(12);
  cover.SetTextSize(0.055);
  cover.AddText("transverse-diffusion.cc");
  cover.AddText("TPC transverse charge-sharing distributions");
  cover.AddText(Form("Input: %s", run_label.Data()));
  TDatime now;
  cover.AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",
                     now.GetYear(), now.GetMonth(), now.GetDay(),
                     now.GetHour(), now.GetMinute(), now.GetSecond()));
  cover.Draw();
  canvas->Print(output_name.c_str());

  const auto draw = [&](const TString& name, const TString& title) {
    auto* source = get_hist(name);
    if (!source) {
      std::cerr << "Missing histogram: " << name << std::endl;
      return;
    }
    canvas->Clear();
    canvas->SetTopMargin(0.10);
    canvas->SetLogz();
    source->SetTitle(Form("%s  [%s]", title.Data(), name.Data()));
    source->Draw("colz");
    canvas->Print(output_name.c_str());
  };

  // Overlay the ProfileX points and generalized-Gaussian fit on the TH2.
  const auto draw_profile_fit = [&](TH2* source, TProfile* profile, TF1* fit,
                                    TPaveText& result) {
    canvas->Clear();
    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    source->Draw("colz");
    profile->Draw("E1 same");
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(3);
    fit->Draw("same");

    TLegend legend(0.14, 0.76, 0.42, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(profile, "ProfileX", "lep");
    legend.AddEntry(fit, "Generalized Gaussian fit", "l");
    legend.Draw();
    result.Draw();
    canvas->Print(output_name.c_str());
  };

  const auto fit_by_cluster_size = [&](const TString& size_name, const TString& size_label) {
    std::vector<double> x, ex, y, ey;
    std::unique_ptr<TH2> all_drift_histograms;
    for (int drift_bin = 0; drift_bin < kNumDriftBins; ++drift_bin) {
      const double drift_low = drift_bin * kDriftBinWidthCm;
      const double drift_high = drift_low + kDriftBinWidthCm;
      const TString name = Form("Transverse_Diffusion_ClusterSize%s_Drift%02d_%02dcm",
                                size_name.Data(), static_cast<int>(drift_low), static_cast<int>(drift_high));
      auto* source = get_hist(name);
      if (!source || source->GetEntries() == 0.) {
        std::cerr << "Missing or empty histogram: " << name << std::endl;
        continue;
      }
      if (!all_drift_histograms) {
        all_drift_histograms.reset(dynamic_cast<TH2*>(source->Clone(
          Form("sum_Transverse_Diffusion_ClusterSize%s_AllDrift", size_name.Data()))));
        all_drift_histograms->SetDirectory(nullptr);
      } else {
        all_drift_histograms->Add(source);
      }
      std::unique_ptr<TProfile> profile(source->ProfileX(Form("p_%s", name.Data()), 1, -1, "s"));
      profile->SetDirectory(nullptr);
      source->SetTitle(Form("TPC transverse diffusion: %s, drift distance = %.0f-%.0f cm;X_{cluster}-X_{pad} [mm];<#it{A}/#it{A}_{sum}",
                            size_label.Data(), drift_low, drift_high));
      profile->SetMarkerStyle(20); profile->SetMarkerColor(kBlack); profile->SetLineColor(kBlack);
      TF1 fit(Form("f_%s", name.Data()),
              "[0]*TMath::Exp(-TMath::Power(TMath::Abs((x-[1])/[2]),[3]))", -10., 10.);
      fit.SetParNames("A", "#mu", "#alpha", "#beta");
      fit.SetParameters(profile->GetMaximum(), 0., 3., 2.);
      fit.SetParLimits(2, 0.01, 100.);
      fit.SetParLimits(3, 0.25, 10.);
      profile->Fit(&fit, "Q0R");
      const double alpha = fit.GetParameter(2);
      const double beta = fit.GetParameter(3);
      const double sigma = GeneralizedGaussianRms(alpha, beta);
      const double sigma_error = alpha > 0. ? fit.GetParError(2) * sigma / alpha : 0.;
      if (sigma > 0. && sigma_error >= 0.) {
        x.push_back(0.5 * (drift_low + drift_high));
        ex.push_back(0.5 * kDriftBinWidthCm);
        y.push_back(sigma * sigma);
        ey.push_back(2. * sigma * sigma_error);
      }
      TPaveText result(0.52, 0.72, 0.88, 0.88, "NDC");
      result.SetBorderSize(0); result.SetFillStyle(0);
      result.AddText(size_label);
      result.AddText(Form("Drift distance = %.0f-%.0f cm", drift_low, drift_high));
      result.AddText(Form("RMS = %.3f #pm %.3f mm", sigma, sigma_error));
      result.AddText(Form("#beta = %.3f #pm %.3f", beta, fit.GetParError(3)));
      draw_profile_fit(source, profile.get(), &fit, result);
    }
    if (all_drift_histograms && all_drift_histograms->GetEntries() > 0.) {
      const TString name = Form("Transverse_Diffusion_ClusterSize%s_AllDrift", size_name.Data());
      std::unique_ptr<TProfile> profile(all_drift_histograms->ProfileX(Form("p_%s", name.Data()), 1, -1, "s"));
      profile->SetDirectory(nullptr);
      all_drift_histograms->SetTitle(Form("TPC transverse diffusion: %s, all drift distances;X_{cluster}-X_{pad} [mm];<#it{A}/#it{A}_{sum}", size_label.Data()));
      profile->SetMarkerStyle(20); profile->SetMarkerColor(kBlack); profile->SetLineColor(kBlack);
      TF1 fit(Form("f_%s", name.Data()), "[0]*TMath::Exp(-TMath::Power(TMath::Abs((x-[1])/[2]),[3]))", -5., 5.);
      fit.SetParNames("A", "#mu", "#alpha", "#beta");
      fit.SetParameters(profile->GetMaximum(), 0., 3., 2.);
      fit.SetParLimits(2, 0.01, 100.);
      fit.SetParLimits(3, 0.25, 10.);
      profile->Fit(&fit, "Q0R");
      const double alpha = fit.GetParameter(2);
      const double beta = fit.GetParameter(3);
      const double rms = GeneralizedGaussianRms(alpha, beta);
      const double rms_error = alpha > 0. ? fit.GetParError(2) * rms / alpha : 0.;
      TPaveText result(0.52, 0.68, 0.88, 0.88, "NDC");
      result.SetBorderSize(0); result.SetFillStyle(0);
      result.AddText(size_label);
      result.AddText("All drift distances");
      result.AddText(Form("RMS = %.3f #pm %.3f mm", rms, rms_error));
      result.AddText(Form("#beta = %.3f #pm %.3f", beta, fit.GetParError(3)));
      draw_profile_fit(all_drift_histograms.get(), profile.get(), &fit, result);
    }
    if (x.size() < 2) return;
    TGraphErrors graph(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    graph.SetTitle(Form("Transverse diffusion: %s;Drift distance [cm];#sigma_{T}^{2} [mm^{2}]", size_label.Data()));
    graph.SetMarkerStyle(20); graph.SetMarkerSize(1.1);
    TF1 diffusion_fit(Form("f_diffusion_%s", size_name.Data()), "[0]+[1]*x",
                      0., kNumDriftBins * kDriftBinWidthCm);
    graph.Fit(&diffusion_fit, "Q0R");
    canvas->Clear(); canvas->SetLogz(0);
    graph.Draw("AP");
    diffusion_fit.SetLineColor(kRed + 1); diffusion_fit.SetLineWidth(3); diffusion_fit.Draw("same");
    const double slope = diffusion_fit.GetParameter(1);
    const double slope_error = diffusion_fit.GetParError(1);
    TPaveText result(0.50, 0.68, 0.88, 0.88, "NDC");
    result.SetBorderSize(0); result.SetFillStyle(0);
    result.AddText(size_label);
    result.AddText("#sigma_{T}^{2}(L)=#sigma_{0}^{2}+D_{T}^{2}L");
    if (slope > 0.) {
      const double diffusion = std::sqrt(slope);
      const double diffusion_error = slope_error / (2. * diffusion);
      result.AddText(Form("D_{T} = %.3f #pm %.3f mm/#sqrt{cm}", diffusion, diffusion_error));
      std::cout << "run0" << runnumbers.front() << ", " << size_label.Data()
                << ": D_T = " << diffusion << " +/- " << diffusion_error << " mm/sqrt(cm)" << std::endl;
    } else {
      result.AddText("D_{T}: fitted slope is not positive");
    }
    result.Draw();
    canvas->Print(output_name.c_str());
  };
  fit_by_cluster_size("03plus", "N_{cluster} > 2");

  canvas->Print((output_name + "]").c_str());
}
