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

enum PIDPreset {
  kPIDPresetE42 = 0,
  kPIDPresetE72 = 1
};

struct PIDCutConfig {
  const char* label;
  const char* output_tag;
  double pion_low;
  double pion_high;
  double kaon_low;
  double kaon_high;
  double proton_low;
  double pi_p_sep_threshold;
  bool use_fit_conversion;
  double fixed_conversion;
};

PIDCutConfig GetPIDCutConfig(PIDPreset preset)
{
  if (preset == kPIDPresetE72) {
    return {"E72", "e72", -3.0, 3.0, -2.0, 2.0, -1.5, 3.0, true, 0.0};
  }

  return {"E42", "e42", -3.0, 3.0, -3.0, 3.0, -4.0, 6.0, false, 12171.3};
}

const char* PIDName(PIDParticle pid)
{
  if (pid == kmyElectron) return "electron";
  if (pid == kmyPion)     return "pion";
  if (pid == kmyKaon)     return "kaon";
  if (pid == kmyProton)   return "proton";
  return "unknown";
}

const char* PIDLatex(PIDParticle pid)
{
  if (pid == kmyElectron) return "#it{e}";
  if (pid == kmyPion)     return "#pi";
  if (pid == kmyKaon)     return "#it{K}";
  if (pid == kmyProton)   return "#it{p}";
  return "";
}

int PIDColor(PIDParticle pid)
{
  if (pid == kmyElectron) return kBlack;
  if (pid == kmyPion)     return kBlue;
  if (pid == kmyKaon)     return kGreen+2;
  if (pid == kmyProton)   return kRed;
  return kGray+2;
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

  if (dedx_pi < dedx_p)
    return dedx_pi + 0.5 * separation_power * sigma_pi;
  else
    return dedx_p + 0.5 * separation_power * sigma_p;
}

double PiProtonSeparationPower(double poq, double conv)
{
  const double dedx_pi = PIDMean(poq, kmyPion, conv);
  const double dedx_p  = PIDMean(poq, kmyProton, conv);

  const double sigma_pi = PIDSigma(poq, kmyPion);
  const double sigma_p  = PIDSigma(poq, kmyProton);
  const double avg_sigma = 0.5 * (sigma_pi + sigma_p);

  if (avg_sigma <= 0.) return 0.;

  return TMath::Abs(dedx_pi - dedx_p) / avg_sigma;
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
    if(y>0)
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

void StylePIDCurve(TGraph* g, PIDParticle pid, int line_style = 1, int line_width = 2)
{
  g->SetLineColor(PIDColor(pid));
  g->SetLineStyle(line_style);
  g->SetLineWidth(line_width);
}

void DrawBasePIDPlot(TH2D* h2, double xmin, double xmax, double ymin, double ymax)
{
  h2->GetXaxis()->SetRangeUser(xmin, xmax);
  h2->GetYaxis()->SetRangeUser(ymin, ymax);
  h2->GetXaxis()->SetTitle("#it{p/z} [GeV/#it{c}]");
  h2->GetYaxis()->SetTitle("TPC #LT#it{dE/dx}#GT (a.u.)");
  h2->Draw("colz");
}

TBox* DrawSeparationPowerRegions(
    double conv,
    double threshold,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    int nscan = 800)
{
  TBox* first_box = nullptr;
  bool in_sep_region = false;
  double sep_region_start = xmin;

  for (int i = 0; i < nscan; ++i) {
    const double poq = xmin + (xmax - xmin) * i / (nscan - 1);
    const bool pass_sep = PiProtonSeparationPower(poq, conv) >= threshold;

    if (pass_sep && !in_sep_region) {
      sep_region_start = poq;
      in_sep_region = true;
    }

    const bool close_region = in_sep_region && (!pass_sep || i == nscan - 1);
    if (close_region) {
      const double sep_region_end = pass_sep
        ? poq
        : xmin + (xmax - xmin) * (i - 1) / (nscan - 1);
      TBox* box = new TBox(sep_region_start, ymin, sep_region_end, ymax);
      box->SetFillColorAlpha(kOrange+7, 0.35);
      box->SetFillStyle(1001);
      box->SetLineColor(kOrange+7);
      box->Draw("same");
      if (!first_box) first_box = box;
      in_sep_region = false;
    }
  }

  return first_box;
}


void SaveStepPDF(TCanvas* c, const char* multipage_pdf, bool first, bool last)
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

bool PassPiProtonSeparationSide(PIDParticle pid, double poq, double y, double conv)
{
  const double dedx_pi = PIDMean(poq, kmyPion, conv);
  const double dedx_p  = PIDMean(poq, kmyProton, conv);
  const double cut = PiProtonSeparationCut(poq, conv);

  if (pid == kmyPion) {
    return dedx_pi < dedx_p ? y <= cut : y >= cut;
  }

  if (pid == kmyProton) {
    return dedx_pi < dedx_p ? y >= cut : y <= cut;
  }

  return true;
}

bool PassPIDWindow(
    PIDParticle pid,
    double poq,
    double y,
    double conv,
    double nsigma_low,
    double nsigma_high,
    double sep_threshold)
{
  if ((pid == kmyPion || pid == kmyProton)
      && PiProtonSeparationPower(poq, conv) >= sep_threshold) {
    return PassPiProtonSeparationSide(pid, poq, y, conv);
  }

  const double low = PIDBoundary(poq, pid, conv, nsigma_low);
  if (pid == kmyProton) return y >= low;

  const double high = PIDBoundary(poq, pid, conv, nsigma_high);
  return y >= TMath::Min(low, high) && y <= TMath::Max(low, high);
}

TH2D* MakePIDWindowHist(
    TH2D* src,
    const char* name,
    PIDParticle pid,
    double conv,
    double nsigma_low,
    double nsigma_high,
    double sep_threshold)
{
  TH2D* h = (TH2D*)src->Clone(name);
  h->SetDirectory(nullptr);
  h->SetTitle("");

  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    const double poq = h->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      const double y = h->GetYaxis()->GetBinCenter(iy);
      if (!PassPIDWindow(pid, poq, y, conv, nsigma_low, nsigma_high, sep_threshold)) {
        h->SetBinContent(ix, iy, 0.);
        h->SetBinError(ix, iy, 0.);
      }
    }
  }

  return h;
}

  
}

void conversion_factor()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);

  //TFile* fin = TFile::Open("result/dedx_sigma_fit.root", "READ");
  TFile* fin = TFile::Open("~/data/JPARC2025Nov_root/gain_calib_260701/run02447_DstTPCHelixTracking.root", "READ");
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
  const double ymax = 400.0;

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

  TF1* f_e = new TF1("f_electron_bethe", ElectronBethe, pstart, pend, 1);
  f_e->SetParName(0, "conversion_factor");
  f_e->SetParameter(0, f_p->GetParameter(0));




  const PIDPreset pid_preset = kPIDPresetE72;
  const PIDCutConfig pid_cut = GetPIDCutConfig(pid_preset);
  const double new_conv = pid_cut.use_fit_conversion
    ? f_p->GetParameter(0)
    : pid_cut.fixed_conversion;

  std::cout << "PID cut preset = " << pid_cut.label
            << ", conversion factor = " << new_conv
            << (pid_cut.use_fit_conversion ? " (fit)" : " (fixed)")
            << ", pi-p separation threshold = " << pid_cut.pi_p_sep_threshold
            << std::endl;

  // PID mean curves
  TGraph* g_pi_mean = MakePIDCurve("g_pi_mean", kmyPion, new_conv, 0.0, xmin, xmax);
  TGraph* g_k_mean  = MakePIDCurve("g_k_mean",  kmyKaon, new_conv, 0.0, xmin, xmax);
  TGraph* g_p_mean  = MakePIDCurve("g_p_mean",  kmyProton, new_conv, 0.0, xmin, xmax);
  TGraph* g_e_mean  = MakePIDCurve("g_e_mean",  kmyElectron, new_conv, 0.0, xmin, xmax);

  // PID boundaries
  TGraph* g_pi_low  = MakePIDCurve("g_pi_low_3sigma",  kmyPion,   new_conv, pid_cut.pion_low, xmin, xmax);
  TGraph* g_pi_high = MakePIDCurve("g_pi_high_3sigma", kmyPion,   new_conv, pid_cut.pion_high, xmin, xmax);

  TGraph* g_k_low   = MakePIDCurve("g_k_low_3sigma",   kmyKaon,   new_conv, pid_cut.kaon_low, xmin, xmax);
  TGraph* g_k_high  = MakePIDCurve("g_k_high_3sigma",  kmyKaon,   new_conv, pid_cut.kaon_high, xmin, xmax);

  TGraph* g_p_low   = MakePIDCurve("g_p_low_3sigma",   kmyProton, new_conv, pid_cut.proton_low, xmin, xmax);
  TGraph* g_p_high  = MakePIDCurve("g_p_high_3sigma",  kmyProton, new_conv, 0.0, xmin, xmax);

  TGraph* g_pi_low_1sigma  = MakePIDCurve("g_pi_low_1sigma",  kmyPion,   new_conv, -1.0, xmin, xmax);
  TGraph* g_pi_high_1sigma = MakePIDCurve("g_pi_high_1sigma", kmyPion,   new_conv,  1.0, xmin, xmax);
  TGraph* g_p_low_1sigma   = MakePIDCurve("g_p_low_1sigma",   kmyProton, new_conv, -1.0, xmin, xmax);
  TGraph* g_p_high_1sigma  = MakePIDCurve("g_p_high_1sigma",  kmyProton, new_conv,  1.0, xmin, xmax);

  // pi/proton separation cut
  TGraph* g_picut = MakePiProtonCutCurve("g_pi_proton_separation_cut", new_conv, xmin, xmax);

  const double pi_p_sep_threshold = pid_cut.pi_p_sep_threshold;
  const int n_sep_scan = 800;
  

  TGraph* g_e_low_1sigma  = MakePIDCurve("g_e_low_1sigma",  kmyElectron, new_conv, -1.0, xmin, xmax);
  TGraph* g_e_high_1sigma = MakePIDCurve("g_e_high_1sigma", kmyElectron, new_conv,  1.0, xmin, xmax);
  TGraph* g_k_low_1sigma  = MakePIDCurve("g_k_low_1sigma",  kmyKaon,     new_conv, -1.0, xmin, xmax);
  TGraph* g_k_high_1sigma = MakePIDCurve("g_k_high_1sigma", kmyKaon,     new_conv,  1.0, xmin, xmax);

  StylePIDCurve(g_e_mean, kmyElectron, 1, 2);
  StylePIDCurve(g_pi_mean, kmyPion, 1, 2);
  StylePIDCurve(g_k_mean, kmyKaon, 1, 2);
  StylePIDCurve(g_p_mean, kmyProton, 1, 2);

  StylePIDCurve(g_e_low_1sigma, kmyElectron, 3, 2);
  StylePIDCurve(g_e_high_1sigma, kmyElectron, 3, 2);
  StylePIDCurve(g_pi_low_1sigma, kmyPion, 3, 2);
  StylePIDCurve(g_pi_high_1sigma, kmyPion, 3, 2);
  StylePIDCurve(g_k_low_1sigma, kmyKaon, 3, 2);
  StylePIDCurve(g_k_high_1sigma, kmyKaon, 3, 2);
  StylePIDCurve(g_p_low_1sigma, kmyProton, 3, 2);
  StylePIDCurve(g_p_high_1sigma, kmyProton, 3, 2);

  StylePIDCurve(g_pi_low, kmyPion, 2, 2);
  StylePIDCurve(g_pi_high, kmyPion, 2, 2);
  StylePIDCurve(g_k_low, kmyKaon, 2, 2);
  StylePIDCurve(g_k_high, kmyKaon, 2, 2);
  StylePIDCurve(g_p_low, kmyProton, 2, 2);
  StylePIDCurve(g_p_high, kmyProton, 2, 2);
  g_picut->SetLineColor(kMagenta+2);
  g_picut->SetLineStyle(9);
  g_picut->SetLineWidth(3);

  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.9);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);

  TString multipage_pdf = Form("result/conversion-factor_%s.pdf", pid_cut.output_tag);
  TCanvas* c = new TCanvas("c_conversion_factor", "conversion factor", 900, 700);
  c->SetLogz();

  DrawBasePIDPlot(h2, xmin, xmax, ymin, ymax);
  g_e_mean->Draw("L same");
  g_pi_mean->Draw("L same");
  g_k_mean->Draw("L same");
  g_p_mean->Draw("L same");
  g_e_low_1sigma->Draw("L same");
  g_e_high_1sigma->Draw("L same");
  g_pi_low_1sigma->Draw("L same");
  g_pi_high_1sigma->Draw("L same");
  g_k_low_1sigma->Draw("L same");
  g_k_high_1sigma->Draw("L same");
  g_p_low_1sigma->Draw("L same");
  g_p_high_1sigma->Draw("L same");
  TLegend* leg_step1 = new TLegend(0.58, 0.58, 0.8, 0.88);
  leg_step1->SetBorderSize(0);
  leg_step1->SetFillColorAlpha(kWhite,0.8);
  leg_step1->SetTextSize(0.03);
  leg_step1->AddEntry(g_e_mean, "#it{e} mean", "l");
  leg_step1->AddEntry(g_pi_mean, "#pi mean", "l");
  leg_step1->AddEntry(g_k_mean, "#it{K} mean", "l");
  leg_step1->AddEntry(g_p_mean, "#it{p} mean", "l");
  leg_step1->AddEntry(g_pi_low_1sigma, "#pm1#sigma", "l");
  leg_step1->Draw();
  SaveStepPDF(c, multipage_pdf.Data(), true, false);

  c->Clear();
  c->SetLogz();
  DrawBasePIDPlot(h2, xmin, xmax, ymin, ymax);
  TBox* sep_box_for_legend = DrawSeparationPowerRegions(new_conv, pi_p_sep_threshold, xmin, xmax, ymin, ymax, n_sep_scan);
  g_pi_mean->Draw("L same");
  g_p_mean->Draw("L same");
  g_picut->Draw("L same");
  TLegend* leg_step2 = new TLegend(0.56, 0.66, 0.8, 0.88);
  leg_step2->SetBorderSize(0);
  leg_step2->SetFillColorAlpha(kWhite,0.8);
  leg_step2->SetTextSize(0.03);
  leg_step2->AddEntry(g_pi_mean, "#pi mean", "l");
  leg_step2->AddEntry(g_p_mean, "#it{p} mean", "l");
  leg_step2->AddEntry(g_picut, "#pi-#it{p} separation cut", "l");
  if (sep_box_for_legend) leg_step2->AddEntry(sep_box_for_legend, Form("sep. power #geq %.1f", pi_p_sep_threshold), "f");
  leg_step2->Draw();
  SaveStepPDF(c, multipage_pdf.Data(), false, false);

  c->Clear();
  c->SetLogz();
  DrawBasePIDPlot(h2, xmin, xmax, ymin, ymax);
  g_pi_mean->Draw("L same");
  g_k_mean->Draw("L same");
  g_p_mean->Draw("L same");
  g_pi_low->Draw("L same");
  g_pi_high->Draw("L same");
  g_k_low->Draw("L same");
  g_k_high->Draw("L same");
  g_p_low->Draw("L same");
  TLegend* leg_step3 = new TLegend(0.55, 0.62, 0.8, 0.88);
  leg_step3->SetBorderSize(0);
  leg_step3->SetFillColorAlpha(kWhite,0.8);
  leg_step3->SetTextSize(0.03);
  leg_step3->AddEntry(g_pi_low, Form("#pi: %.1f#sigma to +%.1f#sigma", pid_cut.pion_low, pid_cut.pion_high), "l");
  leg_step3->AddEntry(g_k_low, Form("#it{K}: %.1f#sigma to +%.1f#sigma", pid_cut.kaon_low, pid_cut.kaon_high), "l");
  leg_step3->AddEntry(g_p_low, Form("#it{p}: #geq %.1f#sigma", pid_cut.proton_low), "l");
  leg_step3->Draw();
  SaveStepPDF(c, multipage_pdf.Data(), false, false);

  c->Clear();
  c->SetLogz();
  DrawBasePIDPlot(h2, xmin, xmax, ymin, ymax);
  g_pi_mean->Draw("L same");
  g_pi_low->Draw("L same");
  g_pi_high->Draw("L same");
  SaveStepPDF(c, multipage_pdf.Data(), false, false);

  c->Clear();
  c->SetLogz();
  DrawBasePIDPlot(h2, xmin, xmax, ymin, ymax);
  g_k_mean->Draw("L same");
  g_k_low->Draw("L same");
  g_k_high->Draw("L same");
  SaveStepPDF(c, multipage_pdf.Data(), false, false);

  c->Clear();
  c->SetLogz();
  DrawBasePIDPlot(h2, xmin, xmax, ymin, ymax);
  g_p_mean->Draw("L same");
  g_p_low->Draw("L same");
  SaveStepPDF(c, multipage_pdf.Data(), false, false);

  TH2D* h_pi_cut = MakePIDWindowHist(h2, "h2_pid_cut_pion", kmyPion, new_conv, pid_cut.pion_low, pid_cut.pion_high, pi_p_sep_threshold);
  TH2D* h_k_cut = MakePIDWindowHist(h2, "h2_pid_cut_kaon", kmyKaon, new_conv, pid_cut.kaon_low, pid_cut.kaon_high, pi_p_sep_threshold);
  // For proton, nsigma_high is ignored below separation-power threshold: there is no upper dE/dx limit.
  TH2D* h_p_cut = MakePIDWindowHist(h2, "h2_pid_cut_proton", kmyProton, new_conv, pid_cut.proton_low, 0.0, pi_p_sep_threshold);

  c->Clear();
  c->SetLogz();
  DrawBasePIDPlot(h_pi_cut, xmin, xmax, ymin, ymax);
  DrawSeparationPowerRegions(new_conv, pi_p_sep_threshold, xmin, xmax, ymin, ymax, n_sep_scan);
  g_pi_low->Draw("L same");
  g_pi_high->Draw("L same");
  g_picut->Draw("L same");
  SaveStepPDF(c, multipage_pdf.Data(), false, false);

  c->Clear();
  c->SetLogz();
  DrawBasePIDPlot(h_k_cut, xmin, xmax, ymin, ymax);
  g_k_low->Draw("L same");
  g_k_high->Draw("L same");
  SaveStepPDF(c, multipage_pdf.Data(), false, false);

  c->Clear();
  c->SetLogz();
  DrawBasePIDPlot(h_p_cut, xmin, xmax, ymin, ymax);
  DrawSeparationPowerRegions(new_conv, pi_p_sep_threshold, xmin, xmax, ymin, ymax, n_sep_scan);
  g_p_low->Draw("L same");
  g_picut->Draw("L same");
  SaveStepPDF(c, multipage_pdf.Data(), false, false);

  c->Clear();
  c->SetLogz();
  DrawBasePIDPlot(h2, xmin, xmax, ymin, ymax);
  DrawSeparationPowerRegions(new_conv, pi_p_sep_threshold, xmin, xmax, ymin, ymax, n_sep_scan);
  f_p->Draw("same");
  g->Draw("P same");
  g_e_mean->Draw("L same");
  g_pi_mean->Draw("L same");
  g_k_mean->Draw("L same");
  g_p_mean->Draw("L same");
  g_pi_low_1sigma->Draw("L same");
  g_pi_high_1sigma->Draw("L same");
  g_p_low_1sigma->Draw("L same");
  g_p_high_1sigma->Draw("L same");
  g_pi_low->Draw("L same");
  g_pi_high->Draw("L same");
  g_k_low->Draw("L same");
  g_k_high->Draw("L same");
  g_p_low->Draw("L same");
  g_picut->Draw("L same");
  TLegend* leg = new TLegend(0.60, 0.58, 0.8, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(kWhite,0.8);
  leg->SetTextSize(0.03);
  leg->AddEntry(g_e_mean, "#it{e}", "l");
  leg->AddEntry(g_pi_mean, "#pi", "l");
  leg->AddEntry(g_pi_low_1sigma, "#pi/#it{p} #pm1#sigma", "l");
  leg->AddEntry(g_k_mean, "#it{K}", "l");
  leg->AddEntry(g_p_mean, "#it{p}", "l");
  leg->AddEntry(g_picut, "#pi-#it{p} separation cut", "l");
  leg->Draw();
  SaveStepPDF(c, multipage_pdf.Data(), false, true);

  auto mg = new TMultiGraph("mg", "mg");
  mg->Add(g_pi_mean);
  mg->Add(g_k_mean);
  mg->Add(g_p_mean);

  mg->Add(g_pi_low);
  mg->Add(g_pi_high);
  mg->Add(g_k_low);
  mg->Add(g_k_high);
  mg->Add(g_p_low);

  TCanvas *c2 = new TCanvas("c2", "c2");
  mg->Draw("L");

  TFile* fout = TFile::Open(Form("result/conversion-factor_%s.root", pid_cut.output_tag), "RECREATE");
  h2->Write("h2_dedx_merged");
  g->Write();
  f_init->Write();
  f_p->Write();
  g_pi_mean->Write();
  g_k_mean->Write();
  g_p_mean->Write();
  g_e_mean->Write();
  g_pi_low->Write();
  g_pi_high->Write();
  g_k_low->Write();
  g_k_high->Write();
  g_p_low->Write();
  g_p_high->Write();
  g_e_low_1sigma->Write();
  g_e_high_1sigma->Write();
  g_pi_low_1sigma->Write();
  g_pi_high_1sigma->Write();
  g_k_low_1sigma->Write();
  g_k_high_1sigma->Write();
  g_p_low_1sigma->Write();
  g_p_high_1sigma->Write();
  g_picut->Write();
  h_pi_cut->Write();
  h_k_cut->Write();
  h_p_cut->Write();
  c->Write();
  c2->Write();
  fout->Close();

  c->cd();
  c->Update();
}
