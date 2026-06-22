// draw_tpc_dedx_pid.C
// Usage:
// root -l
// .L draw_tpc_dedx_pid.C
// draw_tpc_dedx_pid("input.root", "h_dedx_poq");

#include <iostream>

#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"

namespace TPCPID
{
  const Double_t conversion_factor = 12171.3;

  const Double_t sigma_dedx_pi[5] = {3.94842, 0.0138502, -0.110281, 12.6065, -10.9347};
  const Double_t sigma_dedx_k [5] = {6.24543, -3.21037, 1.52683, 127.099, -9.1004};
  const Double_t sigma_dedx_p [5] = {12.9717, -8.43799, 3.10608, 166.494, -6.56123};

  const Double_t sigma_dedx_ep = 7.8511;
  const Double_t sigma_dedx_em = 8.45029;

  const Double_t mpi = 139.57039;   // MeV/c2
  const Double_t mk  = 493.677;     // MeV/c2
  const Double_t mp  = 938.2720813; // MeV/c2
  const Double_t me  = 0.5109989461;// MeV/c2

  Double_t Beta(Double_t energy, Double_t momentum)
  {
    return momentum / energy;
  }

  Double_t DensityEffectCorrection(Double_t betagamma, Double_t *par)
  {
    Double_t constant = 2. * TMath::Log(10.);
    Double_t delta = 0.;
    Double_t X = TMath::Log10(betagamma);

    if (X <= par[2])
      delta = par[5] * TMath::Power(10., 2. * (X - par[2]));
    else if (par[2] < X && X < par[3])
      delta = constant * X - par[4] + par[0] * TMath::Power((par[3] - X), par[1]);
    else if (X >= par[3])
      delta = constant * X - par[4];

    return delta;
  }

  Double_t HypTPCdEdxP10(Double_t mass, Double_t beta)
  {
    // P10 gas
    Double_t rho = TMath::Power(10., -3) * (0.9 * 1.662 + 0.1 * 0.6672);
    Double_t ZoverA = 17.2 / 37.6;
    Double_t I = 0.9 * 188.0 + 0.1 * 41.7;

    Double_t density_effect_par[6];
    density_effect_par[0] = 0.9 * 0.19714 + 0.1 * 0.09253;
    density_effect_par[1] = 0.9 * 2.9618  + 0.1 * 3.6257;
    density_effect_par[2] = 0.9 * 1.7635  + 0.1 * 1.6263;
    density_effect_par[3] = 0.9 * 4.4855  + 0.1 * 3.9716;
    density_effect_par[4] = 0.9 * 11.9480 + 0.1 * 9.5243;
    density_effect_par[5] = 0.;

    Double_t Z = 1.;
    Double_t K = 0.307075;
    Double_t constant = rho * K * ZoverA;

    Double_t I2 = I * I;
    Double_t beta2 = beta * beta;
    Double_t gamma2 = 1. / (1. - beta2);
    Double_t MeVToeV = TMath::Power(10., 6);

    Double_t Wmax =
      2. * me * beta2 * gamma2 /
      (TMath::Sq(me / mass + 1.) + 2. * (me / mass) * (TMath::Sqrt(gamma2) - 1.));

    Double_t delta = DensityEffectCorrection(TMath::Sqrt(beta2 * gamma2),
                                              density_effect_par);

    Double_t dedx =
      constant * Z * Z / beta2 *
      (0.5 * TMath::Log(2. * me * beta2 * gamma2 * Wmax * MeVToeV * MeVToeV / I2)
       - beta2
       - 0.5 * delta);

    return dedx;
  }

  Double_t HypTPCBethe(Double_t poq, Double_t mass)
  {
    Double_t momentum = 1000. * TMath::Abs(poq); // GeV/c -> MeV/c
    if (momentum <= 0.) return 0.;

    Double_t energy = TMath::Hypot(mass, momentum);
    Double_t beta = Beta(energy, momentum);

    return conversion_factor * HypTPCdEdxP10(mass, beta);
  }

  Double_t CalcSigma(const Double_t *par, Double_t poq)
  {
    return par[0] +
           par[1] * TMath::Abs(poq) +
           par[2] * poq * poq +
           par[3] / TMath::Abs(poq) +
           par[4] / (poq * poq);
  }

  Double_t dEdxPi(Double_t poq) { return HypTPCBethe(poq, mpi); }
  Double_t dEdxK (Double_t poq) { return HypTPCBethe(poq, mk);  }
  Double_t dEdxP (Double_t poq) { return HypTPCBethe(poq, mp);  }
  Double_t dEdxE (Double_t poq) { return HypTPCBethe(poq, me);  }

  Double_t SigmaPi(Double_t poq) { return CalcSigma(sigma_dedx_pi, poq); }
  Double_t SigmaK (Double_t poq) { return CalcSigma(sigma_dedx_k,  poq); }
  Double_t SigmaP (Double_t poq) { return CalcSigma(sigma_dedx_p,  poq); }

  Double_t PiPSeparationPower(Double_t poq)
  {
    Double_t sigma_pi = SigmaPi(poq);
    Double_t sigma_p  = SigmaP(poq);

    Double_t avg_sigma = 0.5 * (sigma_pi + sigma_p);
    Double_t dedx_diff = TMath::Abs(dEdxPi(poq) - dEdxP(poq));

    return dedx_diff / avg_sigma;
  }

  Double_t PiPSeparationCut(Double_t poq)
  {
    Double_t dedx_pi = dEdxPi(poq);
    Double_t dedx_p  = dEdxP(poq);

    Double_t sigma_pi = SigmaPi(poq);
    Double_t sigma_p  = SigmaP(poq);

    Double_t separation_power = PiPSeparationPower(poq);

    if (dedx_pi < dedx_p)
      return dedx_pi + 0.5 * separation_power * sigma_pi;
    else
      return dedx_p + 0.5 * separation_power * sigma_p;
  }
}

void SetGraphStyle(TGraph* g, Int_t color, Int_t style, Int_t width = 2)
{
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(width);
}

void tpc_dedx_pid(const char* filename = "/home/had/haein/data/JPARC2025Nov_root/gain_calib_260616/run02602_DstTPCHelixTracking.root",
                       const char* histname = "PID_dEdx_vs_SignedMom",
                       Double_t xmin_user = -1.5,
                       Double_t xmax_user =  1.5)
{
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file: " << filename << std::endl;
    return;
  }

  TH2D* h = dynamic_cast<TH2D*>(f->Get(histname));
  if (!h) {
    std::cerr << "Cannot find TH2D: " << histname << std::endl;
    return;
  }

  TCanvas* c = new TCanvas("c_dedx_pid", "TPC dE/dx PID", 1000, 800);
  c->SetLogz();

  h->GetXaxis()->SetTitle("p/q [GeV/c]");
  h->GetYaxis()->SetTitle("TPC dE/dx [a.u.]");
  h->Draw("colz");

  const Int_t n = 800;

  TGraph* g_pi = new TGraph();
  TGraph* g_k  = new TGraph();
  TGraph* g_p  = new TGraph();

  TGraph* g_pi_p3 = new TGraph();
  TGraph* g_pi_m3 = new TGraph();

  TGraph* g_k_p3 = new TGraph();
  TGraph* g_k_m3 = new TGraph();

  TGraph* g_p_m4 = new TGraph();
  TGraph* g_p_p6 = new TGraph();

  TGraph* g_pip_cut = new TGraph();

  Int_t ip = 0;

  for (Int_t i = 0; i < n; ++i) {
    Double_t poq = xmin_user + (xmax_user - xmin_user) * i / (n - 1.);

    if (TMath::Abs(poq) < 0.03) continue;

    Double_t dedx_pi = TPCPID::dEdxPi(poq);
    Double_t dedx_k  = TPCPID::dEdxK(poq);
    Double_t dedx_p  = TPCPID::dEdxP(poq);

    Double_t sig_pi = TPCPID::SigmaPi(poq);
    Double_t sig_k  = TPCPID::SigmaK(poq);
    Double_t sig_p  = TPCPID::SigmaP(poq);

    g_pi->SetPoint(ip, poq, dedx_pi);
    g_k ->SetPoint(ip, poq, dedx_k);
    g_p ->SetPoint(ip, poq, dedx_p);

    g_pi_p3->SetPoint(ip, poq, dedx_pi + 3. * sig_pi);
    g_pi_m3->SetPoint(ip, poq, dedx_pi - 3. * sig_pi);

    g_k_p3->SetPoint(ip, poq, dedx_k + 3. * sig_k);
    g_k_m3->SetPoint(ip, poq, dedx_k - 3. * sig_k);

    g_p_m4->SetPoint(ip, poq, dedx_p - 4. * sig_p);
    g_p_p6->SetPoint(ip, poq, dedx_p + 6. * sig_p);

    if (TPCPID::PiPSeparationPower(poq) >= 6.)
      g_pip_cut->SetPoint(g_pip_cut->GetN(), poq, TPCPID::PiPSeparationCut(poq));

    ip++;
  }

  SetGraphStyle(g_pi, kBlue + 1, 1, 3);
  SetGraphStyle(g_k,  kGreen + 2, 1, 3);
  SetGraphStyle(g_p,  kRed + 1, 1, 3);

  SetGraphStyle(g_pi_p3, kBlue + 1, 2);
  SetGraphStyle(g_pi_m3, kBlue + 1, 2);

  SetGraphStyle(g_k_p3, kGreen + 2, 2);
  SetGraphStyle(g_k_m3, kGreen + 2, 2);

  SetGraphStyle(g_p_m4, kRed + 1, 2);
  SetGraphStyle(g_p_p6, kRed + 1, 2);

  SetGraphStyle(g_pip_cut, kMagenta + 2, 7, 3);

  g_pi->Draw("L same");
  g_k ->Draw("L same");
  g_p ->Draw("L same");

  g_pi_p3->Draw("L same");
  g_pi_m3->Draw("L same");
  g_k_p3 ->Draw("L same");
  g_k_m3 ->Draw("L same");
  g_p_m4 ->Draw("L same");
  g_p_p6 ->Draw("L same");

  g_pip_cut->Draw("L same");

  TLegend* leg = new TLegend(0.56, 0.62, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(g_pi, "pion expected", "l");
  leg->AddEntry(g_k,  "kaon expected", "l");
  leg->AddEntry(g_p,  "proton expected", "l");
  leg->AddEntry(g_pi_p3, "pion #pm3#sigma", "l");
  leg->AddEntry(g_k_p3,  "kaon #pm3#sigma", "l");
  leg->AddEntry(g_p_m4,  "proton -4#sigma / +6#sigma", "l");
  leg->AddEntry(g_pip_cut, "#pi/p separation cut", "l");
  leg->Draw();

  c->Update();
}
