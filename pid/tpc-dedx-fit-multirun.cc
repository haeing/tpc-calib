#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TObject.h"

namespace DedxFit
{
  const double conversion_factor = 12171.3;

  const double me  = 0.5109989461; // MeV/c2
  const double mpi = 139.57039;
  const double mk  = 493.677;
  const double mp  = 938.2720813;

  double DensityEffectCorrection(double betagamma, double* par)
  {
    const double constant = 2. * TMath::Log(10.);
    const double X = TMath::Log10(betagamma);

    if (X <= par[2])
      return par[5] * TMath::Power(10., 2. * (X - par[2]));

    if (par[2] < X && X < par[3])
      return constant * X - par[4] + par[0] * TMath::Power((par[3] - X), par[1]);

    return constant * X - par[4];
  }

  double HypTPCdEdxP10(double mass, double beta)
  {
    double rho = TMath::Power(10., -3) * (0.9 * 1.662 + 0.1 * 0.6672);
    double ZoverA = 17.2 / 37.6;
    double I = 0.9 * 188.0 + 0.1 * 41.7;

    double density_effect_par[6];
    density_effect_par[0] = 0.9 * 0.19714 + 0.1 * 0.09253;
    density_effect_par[1] = 0.9 * 2.9618  + 0.1 * 3.6257;
    density_effect_par[2] = 0.9 * 1.7635  + 0.1 * 1.6263;
    density_effect_par[3] = 0.9 * 4.4855  + 0.1 * 3.9716;
    density_effect_par[4] = 0.9 * 11.9480 + 0.1 * 9.5243;
    density_effect_par[5] = 0.;

    double Z = 1.;
    double K = 0.307075;
    double constant = rho * K * ZoverA;

    double beta2 = beta * beta;
    double gamma2 = 1. / (1. - beta2);
    double MeVToeV = 1.e6;
    double I2 = I * I;

    double Wmax =
      2. * me * beta2 * gamma2 /
      (TMath::Sq(me / mass + 1.) + 2. * (me / mass) * (TMath::Sqrt(gamma2) - 1.));

    double delta = DensityEffectCorrection(TMath::Sqrt(beta2 * gamma2),
                                            density_effect_par);

    double dedx =
      constant * Z * Z / beta2 *
      (0.5 * TMath::Log(2. * me * beta2 * gamma2 * Wmax * MeVToeV * MeVToeV / I2)
       - beta2
       - 0.5 * delta);

    return conversion_factor * dedx;
  }

  double Bethe(double poq, double mass)
  {
    double mom = 1000. * TMath::Abs(poq); // GeV/c -> MeV/c
    if (mom <= 0.) return 0.;

    double energy = TMath::Hypot(mass, mom);
    double beta = mom / energy;

    return HypTPCdEdxP10(mass, beta);
  }

  double Mass(const std::string& species)
  {
    if (species == "e")  return me;
    if (species == "pi") return mpi;
    if (species == "k")  return mk;
    if (species == "p")  return mp;

    std::cerr << "Unknown species: " << species << std::endl;
    return mpi;
  }
}

TH2D* LoadAndMergeTH2D(const std::vector<int>& run_list,
                       const char* input_dir,
                       const char* file_format,
                       const char* histname)
{
  TH2D* h2_sum = nullptr;

  for (int run : run_list) {
    TString filename = Form(file_format, input_dir, run);

    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "[skip] Cannot open file: " << filename << std::endl;
      if (f) f->Close();
      continue;
    }

    TH2D* h2 = dynamic_cast<TH2D*>(f->Get(histname));
    if (!h2) {
      std::cerr << "[skip] Cannot find TH2D: " << histname
                << " in " << filename << std::endl;
      f->Close();
      continue;
    }

    if (!h2_sum) {
      h2_sum = dynamic_cast<TH2D*>(h2->Clone("h2_dedx_merged"));
      h2_sum->SetDirectory(nullptr);
      h2_sum->Reset();
      h2_sum->SetTitle(Form("Merged %s", histname));
    }

    h2_sum->Add(h2);
    std::cout << "[add] run " << run << " : " << filename
              << " entries = " << h2->GetEntries() << std::endl;

    f->Close();
  }

  if (!h2_sum) {
    std::cerr << "No valid histogram was loaded." << std::endl;
    return nullptr;
  }

  std::cout << "[merged] total entries = " << h2_sum->GetEntries() << std::endl;
  return h2_sum;
}

void FitOneSpecies(TH2D* h2,
                   const std::string& species,
                   double pmin,
                   double pmax,
                   double pstep,
                   double fit_window,
                   TFile* fout)
{
  if (!h2 || !fout) return;

  const double mass = DedxFit::Mass(species);

  TGraphErrors* g_mean  = new TGraphErrors();
  TGraphErrors* g_sigma = new TGraphErrors();

  g_mean->SetName(Form("g_mean_%s", species.c_str()));
  g_sigma->SetName(Form("g_sigma_%s", species.c_str()));

  TCanvas* cfit = new TCanvas(Form("cfit_%s", species.c_str()),
                              Form("dE/dx projection fits: %s", species.c_str()),
                              900, 700);

  int ipoint = 0;

  for (double p = pmin; p < pmax; p += pstep) {
    double p1 = p;
    double p2 = p + pstep;
    double pc = 0.5 * (p1 + p2);

    int binx1 = h2->GetXaxis()->FindBin(p1);
    int binx2 = h2->GetXaxis()->FindBin(p2);

    TH1D* hy = h2->ProjectionY(Form("hy_%s_%03d", species.c_str(), ipoint),
                               binx1, binx2);

    if (hy->GetEntries() < 50) {
      delete hy;
      continue;
    }

    double expected = DedxFit::Bethe(pc, mass);
    double ymin = expected - fit_window;
    double ymax = expected + fit_window;

    TF1* gaus = new TF1(Form("gaus_%s_%03d", species.c_str(), ipoint),
                        "gaus", ymin, ymax);

    gaus->SetParameters(hy->GetMaximum(), expected, 10.);
    hy->Fit(gaus, "RQ0");

    double mean   = gaus->GetParameter(1);
    double sigma  = TMath::Abs(gaus->GetParameter(2));
    double emean  = gaus->GetParError(1);
    double esigma = gaus->GetParError(2);

    if (sigma <= 0 || sigma > 100) {
      delete hy;
      delete gaus;
      continue;
    }

    g_mean->SetPoint(ipoint, pc, mean);
    g_mean->SetPointError(ipoint, 0.5 * pstep, emean);

    g_sigma->SetPoint(ipoint, pc, sigma);
    g_sigma->SetPointError(ipoint, 0.5 * pstep, esigma);

    if (ipoint < 20) {
      cfit->cd();
      hy->Draw();
      gaus->Draw("same");
    }

    ipoint++;

    delete hy;
    delete gaus;
  }

  TF1* fsigma = new TF1(Form("fsigma_%s", species.c_str()),
                        "[0] + [1]*abs(x) + [2]*x*x + [3]*exp([4]*abs(x))",
                        pmin, pmax);

  fsigma->SetParameters(5., 0., 0., 20., -5.);
  if (g_sigma->GetN() >= 5) {
    g_sigma->Fit(fsigma, "R");
  } else {
    std::cerr << "[warning] Not enough sigma points for " << species << std::endl;
  }

  TCanvas* c1 = new TCanvas(Form("c_mean_%s", species.c_str()),
                            Form("Mean dE/dx: %s", species.c_str()),
                            900, 700);

  h2->Draw("colz");

  TF1* fbethe = new TF1(Form("fbethe_%s", species.c_str()),
                        [mass](double* x, double*) {
                          return DedxFit::Bethe(x[0], mass);
                        },
                        pmin, pmax, 0);

  fbethe->SetLineWidth(3);
  fbethe->Draw("same");
  g_mean->SetMarkerStyle(20);
  g_mean->Draw("P same");

  TCanvas* c2 = new TCanvas(Form("c_sigma_%s", species.c_str()),
                            Form("dE/dx sigma: %s", species.c_str()),
                            900, 700);

  g_sigma->SetTitle(Form("%s dE/dx resolution;|p/q| [GeV/c];#sigma_{dE/dx}", species.c_str()));
  g_sigma->SetMarkerStyle(20);
  g_sigma->Draw("AP");
  fsigma->Draw("same");

  std::cout << "\n===== sigma parameter for " << species << " =====" << std::endl;
  for (int i = 0; i < 5; ++i) {
    std::cout << "par[" << i << "] = " << fsigma->GetParameter(i)
              << " +/- " << fsigma->GetParError(i) << std::endl;
  }

  fout->cd();
  g_mean->Write("", TObject::kOverwrite);
  g_sigma->Write("", TObject::kOverwrite);
  fsigma->Write("", TObject::kOverwrite);
  fbethe->Write("", TObject::kOverwrite);
  cfit->Write("", TObject::kOverwrite);
  c1->Write("", TObject::kOverwrite);
  c2->Write("", TObject::kOverwrite);
}

void tpc_dedx_fit_multirun(const std::vector<int>& run_list = {2602, 2603, 2604},
                           const char* input_dir = "/home/had/haein/data/JPARC2025Nov_root/gain_calib_260616",
                           const char* file_format = "%s/run%05d_DstTPCHelixTracking.root",
                           const char* histname = "PID_dEdx_vs_SignedMom",
                           double pmin = 0.10,
                           double pmax = 1.20,
                           double pstep = 0.05,
                           double fit_window = 40.0,
                           const char* outname = "dedx_sigma_fit.root")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TH2D* h2 = LoadAndMergeTH2D(run_list, input_dir, file_format, histname);
  if (!h2) return;

  TFile* fout = TFile::Open(outname, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Cannot create output file: " << outname << std::endl;
    delete h2;
    return;
  }

  fout->cd();
  h2->Write("h2_dedx_merged", TObject::kOverwrite);

  const std::vector<std::string> species_list = {"pi", "k", "p"};
  for (const auto& species : species_list) {
    FitOneSpecies(h2, species, pmin, pmax, pstep, fit_window, fout);
  }

  fout->Close();
  delete h2;

  std::cout << "\nOutput saved to: " << outname << std::endl;
}
