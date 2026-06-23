// conversion-factor.cc
//
// Usage:
// root -l
// .L conversion-factor.cc
// conversion_factor();

#include <iostream>

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"

namespace {
  const double me = 0.5109989461; // MeV/c2
  const double mp = 938.2720813;  // MeV/c2
  const double mpi = 139.57061;   // MeV/c2

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

    const double delta =
      DensityEffectCorrection(TMath::Sqrt(beta2 * gamma2), den);

    const double dedx =
      rho * K * ZoverA / beta2 *
      (0.5 * TMath::Log(2. * me * beta2 * gamma2
                        * Wmax * MeVToeV * MeVToeV / (I * I))
       - beta2
       - 0.5 * delta);

    return dedx;
  }

  double ProtonBethe(double* x, double* par)
  {
    const double p_GeV = x[0];
    if (p_GeV <= 0.) return 0.;

    const double p_MeV = 1000. * p_GeV;
    const double energy = TMath::Hypot(mp, p_MeV);
    const double beta = p_MeV / energy;

    return par[0] * BetheP10Raw(mp, beta);
  }

  double PionBethe(double* x, double* par)
  {
    const double p_GeV = x[0];
    if (p_GeV <= 0.) return 0.;

    const double p_MeV = 1000. * p_GeV;
    const double energy = TMath::Hypot(mpi, p_MeV);
    const double beta = p_MeV / energy;

    return par[0] * BetheP10Raw(mpi, beta);
  }
  double ElectronBethe(double* x, double* par)
  {
    const double p_GeV = x[0];
    if (p_GeV <= 0.) return 0.;

    const double p_MeV = 1000. * p_GeV;
    const double energy = TMath::Hypot(me, p_MeV);
    const double beta = p_MeV / energy;

    return par[0] * BetheP10Raw(me, beta);
  }
  
}

void conversion_factor()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile* fin = TFile::Open("result/dedx_sigma_fit.root", "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open result/dedx_sigma_fit.root" << std::endl;
    return;
  }

  TH2D* h2_in = (TH2D*)fin->Get("h2_dedx_merged");
  if (!h2_in) {
    std::cerr << "Cannot find h2_dedx_merged" << std::endl;
    fin->Close();
    return;
  }

  TH2D* h2 = (TH2D*)h2_in->Clone("h2_dedx_merged_clone");
  h2->SetDirectory(nullptr);
  fin->Close();

  const double xmin = -2.0;
  const double xmax = 2.0;
  const double ymin = 0.0;
  const double ymax = 600.0;

  // proton band가 비교적 분리되어 보이는 구간만 사용
  //const double pstart = 0.15;
  const double pstart = 0;
  //const double pend   = 0.75;
  const double pend   = 1.5;
  const double pstep  = 0.1;

  const int min_entries = 100;

  TF1* f_init = new TF1("f_init_proton_bethe", ProtonBethe, 0.05, 2.0, 1);
  f_init->SetParameter(0, 8500.0);

  TGraphErrors* g = new TGraphErrors();
  g->SetName("g_proton_peak");
  g->SetTitle("Proton peak; q#timesp [GeV/c]; dE/dx peak [a.u.]");

  int ip = 0;

  for (double x1 = pstart; x1 < pend; x1 += pstep) {
    const double x2 = x1 + pstep;
    const double xc = 0.5 * (x1 + x2);

    const int bx1 = h2->GetXaxis()->FindBin(x1);
    const int bx2 = h2->GetXaxis()->FindBin(x2);

    TH1D* py = h2->ProjectionY(Form("py_%03d", ip), bx1, bx2);

    const double expected = f_init->Eval(xc);

    double ylow  = expected * 0.70;
    double yhigh = expected * 1.35;

    //if (ylow < 40.) ylow = 40.;
    if (yhigh > ymax) yhigh = ymax;

    const int by1 = py->GetXaxis()->FindBin(ylow);
    const int by2 = py->GetXaxis()->FindBin(yhigh);

    if (py->Integral(by1, by2) < min_entries) {
      delete py;
      continue;
    }

    TF1* fgaus = new TF1(Form("fgaus_%03d", ip), "gaus", ylow, yhigh);
    fgaus->SetParameters(py->GetMaximum(), expected, 25.0);
    fgaus->SetParLimits(1, ylow, yhigh);
    fgaus->SetParLimits(2, 5.0, 80.0);

    int status = py->Fit(fgaus, "RQ0");

    if (status != 0) {
      delete fgaus;
      delete py;
      continue;
    }

    const double mean  = fgaus->GetParameter(1);
    const double sigma = TMath::Abs(fgaus->GetParameter(2));
    const double emean = fgaus->GetParError(1);

    if (mean < ylow || mean > yhigh) {
      delete fgaus;
      delete py;
      continue;
    }

    if (mean < expected * 0.65) {
      delete fgaus;
      delete py;
      continue;
    }

    if (sigma < 5.0 || sigma > 80.0) {
      delete fgaus;
      delete py;
      continue;
    }

    g->SetPoint(ip, xc, mean);
    g->SetPointError(ip, 0.5 * pstep, emean);

    ip++;

    delete fgaus;
    delete py;
  }

  TF1* f_p = new TF1("f_proton_bethe", ProtonBethe, pstart, pend, 1);
  f_p->SetParName(0, "conversion_factor");
  f_p->SetParameter(0, 8500.);
  f_p->SetParLimits(0, 6000., 10000.);
  f_p->SetLineColor(kRed);
  f_p->SetLineWidth(3);

  g->Fit(f_p, "R");

  std::cout << "\nconversion_factor = "
            << f_p->GetParameter(0)
            << " +/- "
            << f_p->GetParError(0)
            << std::endl;

  TF1* f_pi = new TF1("f_pion_bethe", PionBethe, pstart, pend, 1);
  f_pi->SetParName(0, "conversion_factor");
  f_pi->SetParameter(0, f_p->GetParameter(0));

  TF1* f_e = new TF1("f_pion_bethe", ElectronBethe, pstart, pend, 1);
  f_e->SetParName(0, "conversion_factor");
  f_e->SetParameter(0, f_p->GetParameter(0));
  

  TCanvas* c = new TCanvas("c_conversion_factor", "conversion factor", 900, 700);
  c->SetLogz();

  h2->GetXaxis()->SetRangeUser(xmin, xmax);
  h2->GetYaxis()->SetRangeUser(ymin, ymax);
  h2->GetXaxis()->SetTitle("q#timesp [GeV/c]");
  h2->GetYaxis()->SetTitle("dE/dx (a.u.)");
  h2->Draw("colz");

  f_p->Draw("same");

  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.9);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->Draw("P same");

  f_pi->Draw("same");
  f_e->Draw("same");

  TLegend* leg = new TLegend(0.55, 0.70, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(g, "Proton peak points", "p");
  leg->AddEntry(f_p, "Proton Bethe-Bloch fit", "l");
  leg->Draw();

  TFile* fout = TFile::Open("result/conversion-factor.root", "RECREATE");
  h2->Write("h2_dedx_merged");
  g->Write();
  f_init->Write();
  f_p->Write();
  c->Write();
  fout->Close();

  c->cd();
  c->Update();
}
