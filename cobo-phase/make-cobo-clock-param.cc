#include <iostream>
#include <fstream>
#include <iomanip>

#include "TFile.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"

static constexpr int    NCobo  = 8;
static constexpr double vdrift = 0.055; // mm/ns

double StepFunc(double *x, double *par)
{
  const double xx = x[0];

  const double offset = par[0]; // residual offset [mm]
  const double slope  = par[1]; // residual linear slope [mm/ns]
  const double amp    = par[2]; // residual step amplitude [mm]
  const double phase  = par[3]; // phase offset [ns]
  const double width  = par[4]; // transition width [ns]

  const double arg = (xx - phase) / width;
  const double cdf = 0.5 * (1.0 + TMath::Erf(arg / TMath::Sqrt2()));

  return offset + slope * xx + amp * cdf;
}

void make_cobo_clock_param(const char *infile  = "input.root",
                        const char *histFmt = "h_res_cobo%d",
                        const char *outRoot = "cobo_clock_fit.root",
                        const char *outTxt  = "cobo_clock_param.txt",
                        bool fixWidth = true)
{
  TFile *fin = TFile::Open(infile, "read");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open input file: " << infile << std::endl;
    return;
  }

  TFile *fout = new TFile(outRoot, "recreate");
  std::ofstream ofs(outTxt);

  if (!ofs.is_open()) {
    std::cerr << "Cannot open output txt file: " << outTxt << std::endl;
    return;
  }

  ofs << "########################################\n";
  ofs << "# Coboparameter                         #\n";
  ofs << "# no parameter corresponds to be 50.    #\n";
  ofs << "########################################\n";
  ofs << "#Cobo dummy ATY p0 p1 p2\n";

  for (int icobo = 0; icobo < NCobo; ++icobo) {
    TString hname = Form(histFmt, icobo);
    TH2D *h2 = dynamic_cast<TH2D*>(fin->Get(hname));

    if (!h2) {
      std::cerr << "No histogram found: " << hname << std::endl;
      continue;
    }

    fout->cd();

    TH2D *h2copy = dynamic_cast<TH2D*>(h2->Clone(Form("h2_cobo%d", icobo)));
    h2copy->SetDirectory(fout);

    TProfile *prof =
      h2copy->ProfileX(Form("prof_cobo%d", icobo), 1, -1, "s");

    TF1 *fit = new TF1(Form("fit_cobo%d", icobo),
                       StepFunc,
                       h2copy->GetXaxis()->GetXmin(),
                       h2copy->GetXaxis()->GetXmax(),
                       5);

    double yLeft  = prof->Interpolate(-30.0);
    double yRight = prof->Interpolate( 35.0);
    double amp0   = yRight - yLeft;

    fit->SetParNames("offset_mm",
                     "slope_mm_per_ns",
                     "amp_mm",
                     "p1_ns",
                     "p2_ns");

    fit->SetParameters(yLeft, 0.0, amp0, 20.0, 0.01);

    fit->SetParLimits(2, -30.0, 30.0);
    fit->SetParLimits(3, -45.0, 45.0);

    if (fixWidth) {
      fit->FixParameter(4, 0.01);
    } else {
      fit->SetParLimits(4, 0.01, 10.0);
    }

    prof->Fit(fit, "RQ");

    const double amp_mm = fit->GetParameter(2);
    const double p1     = fit->GetParameter(3);
    const double p2     = fit->GetParameter(4);

    // residual step [mm] -> clock-time correction amplitude [ns]
    const double p0 = amp_mm / vdrift;

    ofs << std::setw(2) << icobo
        << "  0  3  "
        << std::fixed << std::setprecision(8)
        << std::setw(14) << p0 << "  "
        << std::setw(14) << p1 << "  "
        << std::setw(14) << p2 << "\n";

    h2copy->Write();
    prof->Write();
    fit->Write();

    std::cout << "CoBo " << icobo
              << " : p0 = " << p0
              << ", p1 = " << p1
              << ", p2 = " << p2
              << std::endl;
  }

  ofs.close();
  fout->Close();
  fin->Close();

  std::cout << "\nSaved ROOT file      : " << outRoot << std::endl;
  std::cout << "Saved parameter file : " << outTxt  << std::endl;
}
