namespace  {

const double me = 0.5109989461; // MeV/c2
const double mp = 938.2720813;  // MeV/c2
const double mpi = 139.57061;   // MeV/c2

enum PIDParticle {
  kmyElectron = 0,
  kmyPion     = 1,
  kmyKaon     = 2,
  kmyProton   = 3
};
  

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

  double CalcTPCdEdxSigma(const double sigma_par[5], double poq)
{
  const double abspoq = TMath::Abs(poq);
  return sigma_par[0]
       + sigma_par[1] * abspoq
       + sigma_par[2] * poq * poq
       + sigma_par[3] * TMath::Exp(sigma_par[4] * abspoq);
}

double BetheByMass(double poq, double mass_MeV, double conv)
{
  const double p_MeV = 1000. * TMath::Abs(poq);
  if (p_MeV <= 0.) return 0.;

  const double energy = TMath::Hypot(mass_MeV, p_MeV);
  const double beta = p_MeV / energy;

  return conv * BetheP10Raw(mass_MeV, beta);
}

double PIDMean(double poq, PIDParticle pid, double conv)
{
  if (pid == kmyElectron) return BetheByMass(poq, me, conv);
  if (pid == kmyPion)     return BetheByMass(poq, mpi, conv);
  if (pid == kmyProton)   return BetheByMass(poq, mp, conv);

  // kaon mass
  const double mk = 493.677; // MeV/c2
  return BetheByMass(poq, mk, conv);
}
  
double PIDSigma(double poq, PIDParticle pid)
{
  const double sigma_dedx_pi[5] = {3.94842, 0.0138502, -0.110281, 12.6065, -10.9347};
  const double sigma_dedx_k[5]  = {6.24543, -3.21037, 1.52683, 127.099, -9.1004};
  const double sigma_dedx_p[5]  = {12.9717, -8.43799, 3.10608, 166.0, -6.56123};

  if (pid == kmyPion)   return CalcTPCdEdxSigma(sigma_dedx_pi, poq);
  if (pid == kmyKaon)   return CalcTPCdEdxSigma(sigma_dedx_k, poq);
  if (pid == kmyProton) return CalcTPCdEdxSigma(sigma_dedx_p, poq);

  // electron
  return poq > 0 ? 7.8511 : 8.45029;
}

  /*

double PIDSigma(double poq, PIDParticle pid)
{
  const double abspoq = TMath::Abs(poq);

  if (pid == kmyPion) {
    return 4.0 + 1.0 * abspoq;
  }

  if (pid == kmyKaon) {
    return 6.0 + 2.0 * abspoq;
  }

  if (pid == kmyProton) {
    return 8.0 + 4.0 * abspoq;
  }

  if (pid == kmyElectron) {
    return 8.0;
  }

  return 10.0;
}
  */

double PIDBoundary(double poq, PIDParticle pid, double conv, double nsigma)
{
  return PIDMean(poq, pid, conv) + nsigma * PIDSigma(poq, pid);
}

double PiProtonSeparationCut(double poq, double conv)
{
  const double dedx_pi = PIDMean(poq, kmyPion, conv);
  const double dedx_p  = PIDMean(poq, kmyProton, conv);

  const double sigma_pi = PIDSigma(poq, kmyPion);
  const double sigma_p  = PIDSigma(poq, kmyProton);

  const double avg_sigma = 0.5 * (sigma_pi + sigma_p);
  const double separation_power = TMath::Abs(dedx_pi - dedx_p) / avg_sigma;

  // Kinematics.cc와 같은 정의
  if (dedx_pi < dedx_p)
    return dedx_pi + 0.5 * separation_power * sigma_pi;
  else
    return dedx_p + 0.5 * separation_power * sigma_p;
}

TGraph* MakePIDCurve(
    const char* name,
    PIDParticle pid,
    double conv,
    double nsigma,
    double xmin = -2.0,
    double xmax =  2.0,
    int n = 400)
{
  TGraph* g = new TGraph(n);
  g->SetName(name);

  for (int i = 0; i < n; ++i) {
    const double poq = xmin + (xmax - xmin) * i / (n - 1);
    const double y = PIDBoundary(poq, pid, conv, nsigma);
    g->SetPoint(i, poq, y);
  }

  return g;
}

TGraph* MakePiProtonCutCurve(
    const char* name,
    double conv,
    double xmin = -2.0,
    double xmax = 2.0,
    int n = 400)
{
  TGraph* g = new TGraph(n);
  g->SetName(name);

  for (int i = 0; i < n; ++i) {
    const double poq = xmin + (xmax - xmin) * i / (n - 1);
    const double y = PiProtonSeparationCut(poq, conv);
    g->SetPoint(i, poq, y);
  }

  return g;
}
  
}

void conversion_factor()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  //TFile* fin = TFile::Open("result/dedx_sigma_fit.root", "READ");
  TFile* fin = TFile::Open("~/data/JPARC2025Nov_root/physics-735/run02447_DstTPCHelixTracking.root", "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open result/dedx_sigma_fit.root" << std::endl;
    return;
  }

  //TH2D* h2_in = (TH2D*)fin->Get("h2_dedx_merged");
  TH2D* h2_in = (TH2D*)fin->Get("PID_dEdx_vs_SignedMom");
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




  const double new_conv = f_p->GetParameter(0);
  //const double new_conv = 7000;

  // PID mean curves
  TGraph* g_pi_mean = MakePIDCurve("g_pi_mean", kmyPion, new_conv, 0.0, xmin, xmax);
  TGraph* g_k_mean  = MakePIDCurve("g_k_mean",  kmyKaon, new_conv, 0.0, xmin, xmax);
  TGraph* g_p_mean  = MakePIDCurve("g_p_mean",  kmyProton, new_conv, 0.0, xmin, xmax);

  // PID boundaries
  TGraph* g_pi_low  = MakePIDCurve("g_pi_low_3sigma",  kmyPion,   new_conv, -3.0, xmin, xmax);
  TGraph* g_pi_high = MakePIDCurve("g_pi_high_3sigma", kmyPion,   new_conv,  3.0, xmin, xmax);

  TGraph* g_k_low   = MakePIDCurve("g_k_low_3sigma",   kmyKaon,   new_conv, -3.0, xmin, xmax);
  TGraph* g_k_high  = MakePIDCurve("g_k_high_3sigma",  kmyKaon,   new_conv,  3.0, xmin, xmax);

  TGraph* g_p_low   = MakePIDCurve("g_p_low_3sigma",   kmyProton, new_conv, -4.0, xmin, xmax);
  TGraph* g_p_high  = MakePIDCurve("g_p_high_3sigma",  kmyProton, new_conv,  6.0, xmin, xmax);

  // pi/proton separation cut
  TGraph* g_picut = MakePiProtonCutCurve("g_pi_proton_separation_cut", new_conv, xmin, xmax);
  

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

  //f_pi->Draw("same");
  //f_e->Draw("same");

  g_pi_mean->SetLineColor(kBlue);
  g_k_mean->SetLineColor(kGreen+2);
  g_p_mean->SetLineColor(kRed);

  g_pi_low->SetLineColor(kBlue);
  g_pi_high->SetLineColor(kBlue);
  g_k_low->SetLineColor(kGreen+2);
  g_k_high->SetLineColor(kGreen+2);
  g_p_low->SetLineColor(kRed);
  g_p_high->SetLineColor(kRed);
  g_picut->SetLineColor(kMagenta+2);
  
  g_pi_low->SetLineStyle(2);
  g_pi_high->SetLineStyle(2);
  g_k_low->SetLineStyle(2);
  g_k_high->SetLineStyle(2);
  g_p_low->SetLineStyle(2);
  g_p_high->SetLineStyle(2);
  g_picut->SetLineStyle(9);

  g_pi_mean->Draw("L same");
  g_k_mean->Draw("L same");
  g_p_mean->Draw("L same");

  g_pi_low->Draw("L same");
  g_pi_high->Draw("L same");
  g_k_low->Draw("L same");
  g_k_high->Draw("L same");
  g_p_low->Draw("L same");
  g_p_high->Draw("L same");

  g_picut->Draw("L same");
  

  TLegend* leg = new TLegend(0.60, 0.60, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(g_pi_mean, "#pi mean", "l");
  leg->AddEntry(g_pi_low,  "#pi #pm 3#sigma", "l");

  leg->AddEntry(g_k_mean,  "K mean", "l");
  leg->AddEntry(g_k_low,   "K #pm 3#sigma", "l");

  leg->AddEntry(g_p_mean,  "p mean", "l");
  leg->AddEntry(g_p_low,   "p -4#sigma / +6#sigma", "l");

  leg->AddEntry(g_picut,   "#pi-p separation cut", "l");

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
