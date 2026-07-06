#include <glob.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TObjString.h"
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

namespace {

double min_abs_mom = 0.02;
double max_abs_mom = 1.5;
double min_dedx = 10.0;
double max_dedx = 300.0;
const double conversion_factor = 7388.11;

const double me = 0.5109989461; // MeV/c2
const double mpi = 139.57061;
const double mk = 493.677;
const double mp = 938.2720813;

std::vector<TString> GlobFiles(const char* pattern)
{
  std::vector<TString> files;
  glob_t glob_result;
  const int ret = glob(pattern, 0, nullptr, &glob_result);
  if (ret == 0) {
    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
      files.emplace_back(glob_result.gl_pathv[i]);
    }
  }
  globfree(&glob_result);
  return files;
}

double SignedLog10(double x, double xmin)
{
  return TMath::Sign(TMath::Log10(TMath::Abs(x) / xmin), x);
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

void DrawGraph(TGraph* graph)
{
  if (graph) graph->Draw("L same");
}

void WriteGraph(TGraph* graph)
{
  if (graph) graph->Write("", TObject::kOverwrite);
}

} // namespace

void pid_plot_combine(const char* result_subdir = "physics-735-minlayer15",
                      const char* input_pattern = nullptr)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TString pattern = input_pattern;
  if (!input_pattern || TString(input_pattern).IsNull()) {
    pattern = Form("result/%s/pid-dedx-vs-mom-run*.root", result_subdir);
  }

  std::vector<TString> input_files = GlobFiles(pattern.Data());
  if (input_files.empty()) {
    TString fallback = Form("result/%s/pid-dedx-vs-mom.root", result_subdir);
    if (!gSystem->AccessPathName(fallback)) input_files.push_back(fallback);
  }

  if (input_files.empty()) {
    std::cerr << "No input ROOT files matched: " << pattern << std::endl;
    return;
  }

  TH2D* h_sum = nullptr;

  for (const TString& filename : input_files) {
    TFile* fin = TFile::Open(filename, "READ");
    if (!fin || fin->IsZombie()) {
      std::cerr << "[skip] cannot open " << filename << std::endl;
      if (fin) fin->Close();
      continue;
    }

    TH2D* h = dynamic_cast<TH2D*>(fin->Get("h_tpc_dedx_vs_mom"));
    if (!h) {
      std::cerr << "[skip] missing h_tpc_dedx_vs_mom in " << filename << std::endl;
      fin->Close();
      continue;
    }

    if (!h_sum) {
      h_sum = dynamic_cast<TH2D*>(h->Clone("h_tpc_dedx_vs_mom_combined"));
      h_sum->SetDirectory(nullptr);
      h_sum->Reset();
    }

    h_sum->Add(h);
    std::cout << "[add] " << filename << " entries=" << h->GetEntries() << std::endl;
    fin->Close();
  }

  if (!h_sum) {
    std::cerr << "No valid h_tpc_dedx_vs_mom histograms were loaded." << std::endl;
    return;
  }

  TGraph* g_e_neg  = MakeBetheCurve("g_e_bethe_negative",  me,  -1, kBlack);
  TGraph* g_e_pos  = MakeBetheCurve("g_e_bethe_positive",  me,   1, kBlack);
  TGraph* g_pi_neg = MakeBetheCurve("g_pi_bethe_negative", mpi, -1, kBlue);
  TGraph* g_pi_pos = MakeBetheCurve("g_pi_bethe_positive", mpi,  1, kBlue);
  TGraph* g_k_neg  = MakeBetheCurve("g_k_bethe_negative",  mk,  -1, kGreen+2);
  TGraph* g_k_pos  = MakeBetheCurve("g_k_bethe_positive",  mk,   1, kGreen+2);
  TGraph* g_p_neg  = MakeBetheCurve("g_p_bethe_negative",  mp,  -1, kRed);
  TGraph* g_p_pos  = MakeBetheCurve("g_p_bethe_positive",  mp,   1, kRed);

  TString output_pdf = Form("result/%s/pid-dedx-vs-mom-combine.pdf", result_subdir);
  TString output_root = Form("result/%s/pid-dedx-vs-mom-combine.root", result_subdir);

  TString outdir(output_pdf);
  const Ssiz_t slash = outdir.Last('/');
  if (slash != kNPOS) {
    outdir.Remove(slash);
    gSystem->mkdir(outdir, kTRUE);
  }

  TCanvas* c_info = new TCanvas("c_pid_plot_combine_info", "PID plot combine information", 900, 700);
  TPaveText* ptext = new TPaveText(0.08, 0.08, 0.92, 0.92, "NDC");
  ptext->SetFillColor(0);
  ptext->SetBorderSize(1);
  ptext->SetTextAlign(12);
  ptext->SetTextSize(0.022);
  ptext->AddText("pid-plot-combine.cc");
  ptext->AddText("Combined TPC dE/dx vs p/z plot from pid-plot ROOT files");
  ptext->AddText(Form("result_subdir: %s", result_subdir));
  ptext->AddText(Form("input pattern: %s", pattern.Data()));
  ptext->AddText(Form("input files: %zu", input_files.size()));
  for (const TString& input_file : input_files) {
    ptext->AddText(Form("  %s", gSystem->BaseName(input_file)));
  }
  ptext->AddText(Form("combined entries: %.0f", h_sum->GetEntries()));
  ptext->AddText(Form("conversion factor = %.2f", conversion_factor));
  TDatime now;
  ptext->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",
                      now.GetYear(), now.GetMonth(), now.GetDay(),
                      now.GetHour(), now.GetMinute(), now.GetSecond()));
  ptext->Draw();
  SavePDF(c_info, output_pdf.Data(), true, false);

  const double xmax_log = TMath::Log10(max_abs_mom / min_abs_mom);
  TCanvas* c = new TCanvas("c_tpc_dedx_vs_mom_combined", "Combined TPC dE/dx vs momentum", 1000, 700);
  c->SetLogy();
  c->SetLogz();
  c->SetBottomMargin(0.20);

  h_sum->SetMaximum(10000);
  h_sum->GetXaxis()->SetRangeUser(-xmax_log, xmax_log);
  h_sum->GetYaxis()->SetRangeUser(min_dedx, max_dedx);
  h_sum->GetXaxis()->SetLabelSize(0.);
  h_sum->GetXaxis()->SetTickLength(0.);
  h_sum->GetXaxis()->SetTitleOffset(1.55);
  h_sum->GetXaxis()->SetTitle("#it{p/z} [GeV/#it{c}]");
  h_sum->GetYaxis()->SetTitle("TPC #LT#it{dE/dx}#GT (a.u.)");
  h_sum->Draw("colz");

  TLine* center_line = new TLine(0., min_dedx, 0., max_dedx);
  center_line->SetLineColor(kBlack);
  center_line->SetLineWidth(1);
  center_line->Draw("same");

  //DrawGraph(g_e_neg);
  //DrawGraph(g_e_pos);
  DrawGraph(g_pi_neg);
  DrawGraph(g_pi_pos);
  DrawGraph(g_k_neg);
  DrawGraph(g_k_pos);
  DrawGraph(g_p_neg);
  DrawGraph(g_p_pos);

  TLegend* leg = new TLegend(0.72, 0.68, 0.8, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite, 0.75);
  leg->SetTextSize(0.030);
  //leg->AddEntry(g_e_pos, "#it{e}", "l");
  leg->AddEntry(g_pi_pos, "#pi", "l");
  leg->AddEntry(g_k_pos, "#it{K}", "l");
  leg->AddEntry(g_p_pos, "#it{p}", "l");
  leg->Draw();

  DrawSignedLogXAxis(min_abs_mom, max_abs_mom);
  SavePDF(c, output_pdf.Data(), false, true);

  TFile* fout = TFile::Open(output_root, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Cannot create " << output_root << std::endl;
    return;
  }

  fout->cd();
  h_sum->Write("", TObject::kOverwrite);
  WriteGraph(g_e_neg);
  WriteGraph(g_e_pos);
  WriteGraph(g_pi_neg);
  WriteGraph(g_pi_pos);
  WriteGraph(g_k_neg);
  WriteGraph(g_k_pos);
  WriteGraph(g_p_neg);
  WriteGraph(g_p_pos);
  c_info->Write("", TObject::kOverwrite);
  c->Write("", TObject::kOverwrite);
  fout->Close();

  std::cout << "pid_plot_combine: wrote " << output_pdf << std::endl;
  std::cout << "pid_plot_combine: wrote " << output_root << std::endl;
}
