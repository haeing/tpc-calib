#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TDatime.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TObject.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

const std::vector<int> runnumbers = {2457, 2458, 2459, 2460, 2462, 2463, 2465, 2466, 2468};
//const std::vector<int> runnumbers = {2447};

double min_abs_mom = 0.02;
double max_abs_mom = 1.5;
double min_dedx = 10.0;
double max_dedx = 300.0;
double max_chisqr = 2.0;
double max_close_dist_tpc = 5.0;
const double conversion_factor = 7388.11;

namespace {

const double me = 0.5109989461; // MeV/c2
const double mpi = 139.57061;
const double mk = 493.677;
const double mp = 938.2720813;

std::vector<double> LogBins(int nbins, double xmin, double xmax)
{
  std::vector<double> bins(nbins + 1);
  const double log_min = TMath::Log10(xmin);
  const double log_max = TMath::Log10(xmax);

  for (int i = 0; i <= nbins; ++i) {
    const double t = static_cast<double>(i) / nbins;
    bins[i] = TMath::Power(10., log_min + (log_max - log_min) * t);
  }

  return bins;
}

double SignedLog10(double x, double xmin)
{
  return TMath::Sign(TMath::Log10(TMath::Abs(x) / xmin), x);
}

double DensityEffectCorrection(double betagamma, double* par)
{
  const double c = 2. * TMath::Log(10.);
  const double X = TMath::Log10(betagamma);

  if (X <= par[2])
    return par[5] * TMath::Power(10., 2. * (X - par[2]));

  if (X < par[3])
    return c * X - par[4] + par[0] * TMath::Power(par[3] - X, par[1]);

  return c * X - par[4];
}

double BetheP10Raw(double mass, double beta)
{
  const double rho = 1.e-3 * (0.9 * 1.662 + 0.1 * 0.6672);
  const double ZoverA = 17.2 / 37.6;
  const double I = 0.9 * 188.0 + 0.1 * 41.7; // eV

  double den[6] = {
    0.9 * 0.19714 + 0.1 * 0.09253,
    0.9 * 2.9618  + 0.1 * 3.6257,
    0.9 * 1.7635  + 0.1 * 1.6263,
    0.9 * 4.4855  + 0.1 * 3.9716,
    0.9 * 11.9480 + 0.1 * 9.5243,
    0.
  };

  const double K = 0.307075;
  const double beta2 = beta * beta;
  const double gamma2 = 1. / (1. - beta2);
  const double MeVToeV = 1.e6;

  const double Wmax =
    2. * me * beta2 * gamma2 /
    (TMath::Sq(me / mass + 1.) +
     2. * (me / mass) * (TMath::Sqrt(gamma2) - 1.));

  const double delta = DensityEffectCorrection(TMath::Sqrt(beta2 * gamma2), den);

  return rho * K * ZoverA / beta2 *
    (0.5 * TMath::Log(2. * me * beta2 * gamma2
                      * Wmax * MeVToeV * MeVToeV / (I * I))
     - beta2
     - 0.5 * delta);
}

double BetheByMass(double poq, double mass_MeV)
{
  const double p_MeV = 1000. * TMath::Abs(poq);
  if (p_MeV <= 0.) return 0.;

  const double energy = TMath::Hypot(mass_MeV, p_MeV);
  const double beta = p_MeV / energy;

  return conversion_factor * BetheP10Raw(mass_MeV, beta);
}

TGraph* MakeBetheCurve(const char* name, double mass_MeV, int charge_sign,
                       int color, int n = 400)
{
  TGraph* g = new TGraph(n);
  g->SetName(name);
  g->SetLineColor(color);
  g->SetLineWidth(2);

  for (int i = 0; i < n; ++i) {
    const double t = static_cast<double>(i) / (n - 1);
    const double p = min_abs_mom * TMath::Power(max_abs_mom / min_abs_mom, t);
    const double poq = charge_sign * p;
    g->SetPoint(i, SignedLog10(poq, min_abs_mom), BetheByMass(poq, mass_MeV));
  }

  return g;
}

template <class T>
bool SetBranch(TChain& chain, const char* name, T* address, bool required = true)
{
  if (!chain.GetBranch(name)) {
    if (required) std::cerr << "pid_plot: missing branch " << name << std::endl;
    return !required;
  }

  chain.SetBranchAddress(name, address);
  return true;
}

double MinCloseDistTpc(const std::vector<std::vector<double>>* close_dist)
{
  if (!close_dist) return -1.;

  double min_dist = -1.;
  for (const auto& pair_dist : *close_dist) {
    for (double dist : pair_dist) {
      if (!std::isfinite(dist)) continue;
      if (min_dist < 0. || dist < min_dist) min_dist = dist;
    }
  }

  return min_dist;
}

bool HasTrackValue(const std::vector<double>* values, size_t itrack)
{
  return values && itrack < values->size() && std::isfinite(values->at(itrack));
}

bool HasTrackValue(const std::vector<int>* values, size_t itrack)
{
  return values && itrack < values->size();
}

std::string RunListString(const std::vector<int>& runs)
{
  std::ostringstream oss;
  for (size_t i = 0; i < runs.size(); ++i) {
    if (i) oss << ", ";
    oss << runs[i];
  }
  return oss.str();
}

void SavePDF(TCanvas* c, const char* multipage_pdf, bool first, bool last)
{
  c->Update();

  if (first && last) {
    c->Print(multipage_pdf);
  } else if (first) {
    c->Print(Form("%s(", multipage_pdf));
  } else if (last) {
    c->Print(Form("%s)", multipage_pdf));
  } else {
    c->Print(multipage_pdf);
  }
}

void DrawSignedLogXAxis(double min_abs_mom, double max_abs_mom)
{
  const std::vector<double> decades = {0.01, 0.1, 1.0};
  const double xmax_log = TMath::Log10(max_abs_mom / min_abs_mom);
  const double xmin = -xmax_log;
  const double xmax =  xmax_log;

  const double left = gPad->GetLeftMargin();
  const double right = 1. - gPad->GetRightMargin();
  const double bottom = gPad->GetBottomMargin();
  const double width = right - left;

  auto ToNDC = [&](double signed_log_x) {
    return left + width * (signed_log_x - xmin) / (xmax - xmin);
  };

  TLine line;
  line.SetLineWidth(1);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAlign(23);
  latex.SetTextSize(0.045);

  auto DrawTick = [&](double p, bool label) {
    if (p < min_abs_mom || p > max_abs_mom) return;

    const double xp = ToNDC(SignedLog10(p, min_abs_mom));
    const double xn = ToNDC(SignedLog10(-p, min_abs_mom));
    const double len = label ? 0.026 : 0.014;

    line.DrawLineNDC(xn, bottom, xn, bottom + len);
    line.DrawLineNDC(xp, bottom, xp, bottom + len);

    if (label) {
      const char* text = TMath::Abs(p - 0.1) < 1.e-9 ? "10^{-1}" : "1";
      latex.DrawLatex(xn, bottom - 0.03, Form("-%s", text));
      latex.DrawLatex(xp, bottom - 0.03, text);
    }
  };

  for (double decade : decades) {
    for (int m = 1; m < 10; ++m) {
      const double p = m * decade;
      const bool label = TMath::Abs(p - 0.1) < 1.e-9 || TMath::Abs(p - 1.0) < 1.e-9;
      DrawTick(p, label);
    }
  }
}

} // namespace

void pid_plot(const char* result_subdir = "physics-735-minlayer15",
              const char* tree_name = "tpc")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const std::string input_dir = Form("~/data/JPARC2025Nov_root/%s", result_subdir);
  const char* input_file_format = "%s/run%05d_DstTPCHelixTracking.root";
  TString output_pdf = Form("result/%s/pid-dedx-vs-mom.pdf", result_subdir);
  TString output_root = Form("result/%s/pid-dedx-vs-mom.root", result_subdir);
  if (runnumbers.size() == 1){
    output_pdf = Form("result/%s/pid-dedx-vs-mom-run%05d.pdf", result_subdir,runnumbers.front());
    output_root = Form("result/%s/pid-dedx-vs-mom-run%05d.root", result_subdir,runnumbers.front());

  }
  else{
    output_pdf = Form("result/%s/pid-dedx-vs-mom-run%05d-%05d-n%zu.pdf",
                  result_subdir, runnumbers.front(), runnumbers.back(), runnumbers.size());
    output_root = Form("result/%s/pid-dedx-vs-mom-run%05d-%05d-n%zu.root",
                  result_subdir, runnumbers.front(), runnumbers.back(), runnumbers.size());
  }

  TString outdir(output_pdf);
  const Ssiz_t slash = outdir.Last('/');
  if (slash != kNPOS) {
    outdir.Remove(slash);
    gSystem->mkdir(outdir, kTRUE);
  }

  TChain chain(tree_name);
  int nfiles = 0;
  for (int run : runnumbers) {
    TString input = Form(input_file_format, input_dir.c_str(), run);
    const int added = chain.Add(input);
    if (added > 0) {
      nfiles += added;
      std::cout << "[add] run " << run << " : " << input << std::endl;
    } else {
      std::cerr << "[skip] run " << run << " : " << input << std::endl;
    }
  }

  if (nfiles <= 0 || chain.GetEntries() <= 0) {
    std::cerr << "pid_plot: no entries found from " << input_dir
              << " with tree " << tree_name << std::endl;
    return;
  }

  std::vector<double>* mom0 = nullptr;
  std::vector<double>* dEdx = nullptr;
  std::vector<int>* charge = nullptr;
  std::vector<double>* chisqr = nullptr;
  std::vector<std::vector<double>>* closeDistTpc = nullptr;

  bool ok = true;
  ok &= SetBranch(chain, "mom0", &mom0);
  ok &= SetBranch(chain, "dEdx", &dEdx);
  ok &= SetBranch(chain, "charge", &charge);
  ok &= SetBranch(chain, "chisqr", &chisqr);
  if (chain.GetBranch("closeDistTpc")) {
    ok &= SetBranch(chain, "closeDistTpc", &closeDistTpc);
  } else {
    ok &= SetBranch(chain, "closeDistTPC", &closeDistTpc, max_close_dist_tpc > 0.);
  }
  if (!ok) return;

  const int nx = 240;
  const int ny = 240;
  const double xmax_log = TMath::Log10(max_abs_mom / min_abs_mom);
  const std::vector<double> ybins = LogBins(ny, min_dedx, max_dedx);

  TH2D* h_dedx =
    new TH2D("h_tpc_dedx_vs_mom",
             "TPC dE/dx vs momentum;#it{p/z} [GeV/#it{c}];TPC #LT#it{dE/dx}#GT (a.u.)",
             nx, -xmax_log, xmax_log, ny, ybins.data());
  h_dedx->SetDirectory(nullptr);
  h_dedx->SetMaximum(10000);
  h_dedx->GetXaxis()->SetLabelSize(0.);
  h_dedx->GetXaxis()->SetTickLength(0.);
  h_dedx->GetXaxis()->SetTitleOffset(1.55);

  Long64_t n_tracks = 0;
  Long64_t n_filled = 0;
  Long64_t n_zero_mom = 0;
  Long64_t n_cut = 0;

  const Long64_t nentries = chain.GetEntries();
  for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
    //for (Long64_t ientry = 0; ientry < 1000; ++ientry) {
    chain.GetEntry(ientry);

    const double min_close_dist = MinCloseDistTpc(closeDistTpc);
    if (max_close_dist_tpc > 0. &&
        (min_close_dist < 0. /*|| min_close_dist > max_close_dist_tpc*/)) {
      continue;
    }

    if (!mom0 || !charge || !dEdx) continue;

    const size_t ntrack = std::min({mom0->size(), charge->size(), dEdx->size()});

    for (size_t itrack = 0; itrack < ntrack; ++itrack) {
      ++n_tracks;

      if (!HasTrackValue(mom0, itrack) ||
          !HasTrackValue(charge, itrack) ||
          !HasTrackValue(dEdx, itrack) ||
          !HasTrackValue(chisqr, itrack)) {
        ++n_cut;
        continue;
      }

      const double poq = mom0->at(itrack) * charge->at(itrack);
      const double dedx = dEdx->at(itrack);
      const double abs_poq = TMath::Abs(poq);

      if (poq == 0.) {
        ++n_zero_mom;
        continue;
      }

      if (abs_poq < min_abs_mom || abs_poq > max_abs_mom ||
          dedx <= min_dedx || dedx >= max_dedx) {
        ++n_cut;
        continue;
      }

      if (max_chisqr > 0. && chisqr->at(itrack) > max_chisqr) {
        ++n_cut;
        continue;
      }

      h_dedx->Fill(SignedLog10(poq, min_abs_mom), dedx);
      ++n_filled;
    }
  }

  TGraph* g_e_neg  = MakeBetheCurve("g_e_bethe_negative",  me,  -1, kBlack);
  TGraph* g_e_pos  = MakeBetheCurve("g_e_bethe_positive",  me,   1, kBlack);
  TGraph* g_pi_neg = MakeBetheCurve("g_pi_bethe_negative", mpi, -1, kBlue);
  TGraph* g_pi_pos = MakeBetheCurve("g_pi_bethe_positive", mpi,  1, kBlue);
  TGraph* g_k_neg  = MakeBetheCurve("g_k_bethe_negative",  mk,  -1, kGreen+2);
  TGraph* g_k_pos  = MakeBetheCurve("g_k_bethe_positive",  mk,   1, kGreen+2);
  TGraph* g_p_neg  = MakeBetheCurve("g_p_bethe_negative",  mp,  -1, kRed);
  TGraph* g_p_pos  = MakeBetheCurve("g_p_bethe_positive",  mp,   1, kRed);

  TCanvas* c_info = new TCanvas("c_pid_plot_info", "PID plot information", 900, 700);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("pid-plot.cc");
  p->AddText("TPC dEdx vs p*z plot in log scale");
  p->AddText("chisqr < 2 cut");
  p->AddText(Form("conversion factor = %.2f", conversion_factor));
  if (runnumbers.size() == 1)
    p->AddText(Form("run%05d", runnumbers.front()));
  else
    p->AddText(Form("runs %05d-%05d, n=%zu", runnumbers.front(), runnumbers.back(), runnumbers.size()));
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  SavePDF(c_info, output_pdf.Data(), true, false);

  TCanvas* c = new TCanvas("c_tpc_dedx_vs_mom", "TPC dE/dx vs momentum", 1000, 700);
  c->SetLogy();
  c->SetLogz();
  c->SetBottomMargin(0.20);
  h_dedx->GetXaxis()->SetRangeUser(-xmax_log, xmax_log);
  h_dedx->GetYaxis()->SetRangeUser(min_dedx, max_dedx);
  h_dedx->Draw("colz");

  TLine* center_line = new TLine(0., min_dedx, 0., max_dedx);
  center_line->SetLineColor(kBlack);
  center_line->SetLineWidth(1);
  center_line->Draw("same");

  g_e_neg->Draw("L same");
  g_e_pos->Draw("L same");
  g_pi_neg->Draw("L same");
  g_pi_pos->Draw("L same");
  g_k_neg->Draw("L same");
  g_k_pos->Draw("L same");
  g_p_neg->Draw("L same");
  g_p_pos->Draw("L same");

  TLegend* leg = new TLegend(0.72, 0.68, 0.8, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite, 0.75);
  leg->SetTextSize(0.030);
  leg->AddEntry(g_e_pos, "#it{e}", "l");
  leg->AddEntry(g_pi_pos, "#pi", "l");
  leg->AddEntry(g_k_pos, "#it{K}", "l");
  leg->AddEntry(g_p_pos, "#it{p}", "l");
  leg->Draw();

  DrawSignedLogXAxis(min_abs_mom, max_abs_mom);
  SavePDF(c, output_pdf.Data(), false, true);

  TFile* fout = TFile::Open(output_root, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "pid_plot: cannot create " << output_root << std::endl;
    return;
  }

  fout->cd();
  h_dedx->Write("", TObject::kOverwrite);
  g_e_neg->Write("", TObject::kOverwrite);
  g_e_pos->Write("", TObject::kOverwrite);
  g_pi_neg->Write("", TObject::kOverwrite);
  g_pi_pos->Write("", TObject::kOverwrite);
  g_k_neg->Write("", TObject::kOverwrite);
  g_k_pos->Write("", TObject::kOverwrite);
  g_p_neg->Write("", TObject::kOverwrite);
  g_p_pos->Write("", TObject::kOverwrite);
  c_info->Write("", TObject::kOverwrite);
  c->Write("", TObject::kOverwrite);
  fout->Close();

  std::cout << "pid_plot: input files = " << nfiles << std::endl;
  std::cout << "pid_plot: entries     = " << nentries << std::endl;
  std::cout << "pid_plot: tracks      = " << n_tracks << std::endl;
  std::cout << "pid_plot: filled      = " << n_filled << std::endl;
  std::cout << "pid_plot: zero p/z    = " << n_zero_mom << std::endl;
  std::cout << "pid_plot: cut/skipped = " << n_cut << std::endl;
  std::cout << "pid_plot: wrote " << output_pdf << std::endl;
  std::cout << "pid_plot: wrote " << output_root << std::endl;
}
