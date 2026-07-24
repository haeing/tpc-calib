#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace {
constexpr const char *kDefaultMap =
    "param_history/"
    "ShsFieldMap_20210526_Extrapolated_HS0";

double EdgeMin(double first, double step)
{
  return first - 0.5 * std::abs(step);
}

double EdgeMax(double first, int nbin, double step)
{
  return first + (nbin - 0.5) * step;
}

void SetHistStyle(TH1D *h, int color, int style)
{
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetLineStyle(style);
  h->SetLineWidth(2);
}
} // namespace

void draw_fldmap(const char *mapPath = kDefaultMap,
                 const char *outDir = "result")
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);

  std::ifstream fin(mapPath);
  if (!fin) {
    std::cerr << "Cannot open field map: " << mapPath << std::endl;
    return;
  }

  int nx = 0;
  int ny = 0;
  int nz = 0;
  double x0 = 0.0;
  double y0 = 0.0;
  double z0 = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  fin >> nx >> ny >> nz >> x0 >> y0 >> z0 >> dx >> dy >> dz;

  if (!fin || nx <= 0 || ny <= 0 || nz <= 0 || dx == 0.0 || dy == 0.0 ||
      dz == 0.0) {
    std::cerr << "Bad field map header in: " << mapPath << std::endl;
    return;
  }

  const double zmin = EdgeMin(z0, dz);
  const double zmax = EdgeMax(z0, nz, dz);
  const double targetX = 0.0;
  const double targetY = 0.0;
  const int ix0 = std::lround((targetX - x0) / dx);
  const int iy0 = std::lround((targetY - y0) / dy);
  const double sliceX = x0 + ix0 * dx;
  const double sliceY = y0 + iy0 * dy;
  const double tolX = 0.5 * std::abs(dx);
  const double tolY = 0.5 * std::abs(dy);

  TH1D *hBx = new TH1D("hBx_z_x0_y0", "Field vs Z at X=0, Y=0;Z [cm];Field [T]",
                       nz, zmin, zmax);
  TH1D *hBy = new TH1D("hBy_z_x0_y0", "Field vs Z at X=0, Y=0;Z [cm];Field [T]",
                       nz, zmin, zmax);
  TH1D *hBz = new TH1D("hBz_z_x0_y0", "Field vs Z at X=0, Y=0;Z [cm];Field [T]",
                       nz, zmin, zmax);
  TH1D *hBmag = new TH1D("hBmag_z_x0_y0", "Field vs Z at X=0, Y=0;Z [cm];Field [T]",
                         nz, zmin, zmax);
  hBx->Sumw2(false);
  hBy->Sumw2(false);
  hBz->Sumw2(false);
  hBmag->Sumw2(false);

  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double bx = 0.0;
  double by = 0.0;
  double bz = 0.0;
  Long64_t nread = 0;
  Long64_t nslice = 0;

  while (fin >> x >> y >> z >> bx >> by >> bz) {
    ++nread;
    if (std::abs(x - sliceX) > tolX || std::abs(y - sliceY) > tolY) {
      continue;
    }

    const int bin = hBmag->FindBin(z);
    const double bmag = std::sqrt(bx * bx + by * by + bz * bz);
    hBx->SetBinContent(bin, bx);
    hBy->SetBinContent(bin, by);
    hBz->SetBinContent(bin, bz);
    hBmag->SetBinContent(bin, bmag);
    ++nslice;
  }

  std::cout << "Read " << nread << " field points from " << mapPath
            << std::endl;
  std::cout << "Grid: nx=" << nx << " ny=" << ny << " nz=" << nz
            << " start=(" << x0 << ", " << y0 << ", " << z0 << ")"
            << " step=(" << dx << ", " << dy << ", " << dz << ")"
            << std::endl;
  std::cout << "Z scan uses nearest point to X=0,Y=0: X=" << sliceX
            << " Y=" << sliceY << " with " << nslice << " points"
            << std::endl;

  SetHistStyle(hBx, kRed + 1, 2);
  SetHistStyle(hBy, kBlue + 1, 3);
  SetHistStyle(hBz, kGreen + 2, 1);
  SetHistStyle(hBmag, kBlack, 1);

  double ymin = hBx->GetMinimum();
  double ymax = hBx->GetMaximum();
  for (TH1D *h : {hBy, hBz, hBmag}) {
    ymin = std::min(ymin, h->GetMinimum());
    ymax = std::max(ymax, h->GetMaximum());
  }
  const double margin = 0.08 * (ymax - ymin == 0.0 ? 1.0 : ymax - ymin);
  hBmag->SetMinimum(ymin - margin);
  hBmag->SetMaximum(ymax + margin);
  hBmag->SetTitle(Form("Field vs Z at X=%.0f, Y=%.0f;Z [cm];Field [T]", sliceX,
                       sliceY));

  TCanvas *c = new TCanvas("c_fldmap_zscan", "Field vs Z at X=0,Y=0",
                           1000, 700);
  hBmag->Draw("HIST L");
  hBx->Draw("HIST L SAME");
  hBy->Draw("HIST L SAME");
  hBz->Draw("HIST L SAME");

  TLegend *leg = new TLegend(0.75, 0.70, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hBmag, "|B|", "l");
  leg->AddEntry(hBx, "B_{x}", "l");
  leg->AddEntry(hBy, "B_{y}", "l");
  leg->AddEntry(hBz, "B_{z}", "l");
  leg->Draw();

  c->SaveAs((std::string(outDir) + "/fldmap_zscan_x0_y0_th1d_update.pdf").c_str());
  std::cout << "Wrote " << outDir << "/fldmap_zscan_x0_y0_th1d_update.pdf"
            << std::endl;
}
