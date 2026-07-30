#include "../TPCEventDisplayHelper.hh"

#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr const char *kDefaultInputDir = "~/data/JPARC2025Nov_root/blc";

template <class Consumer>
Long64_t VisitNestedTracks(TTree *tree, Consumer consume)
{
  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<std::vector<double>>> xvpHS(reader, "xvpHS");
  TTreeReaderValue<std::vector<std::vector<double>>> zvpHS(reader, "zvpHS");
  Long64_t tracks = 0;

  while (reader.Next()) {
    const auto &xTracks = *xvpHS;
    const auto &zTracks = *zvpHS;
    const size_t nTracks = std::min(xTracks.size(), zTracks.size());
    for (size_t track = 0; track < nTracks; ++track) {
      consume(zTracks[track], xTracks[track]);
      ++tracks;
    }
  }
  return tracks;
}

Long64_t FillInterpolatedTPCPads(const std::vector<double> &zValues,
                                 const std::vector<double> &xValues,
                                 TH2Poly *padMap)
{
  std::vector<std::pair<double, double>> points;
  const size_t nValues = std::min(zValues.size(), xValues.size());
  for (size_t i = 0; i < nValues; ++i) {
    if (std::isfinite(zValues[i]) && std::isfinite(xValues[i])) {
      points.emplace_back(zValues[i], xValues[i]);
    }
  }
  if (points.size() < 2) {
    return 0;
  }
  std::sort(points.begin(), points.end());

  // Piecewise interpolation between VPHS points; no global straight-line fit.
  constexpr double sampleStep = 0.25; // mm, smaller than a TPC pad width
  std::set<int> filledBins;
  for (size_t point = 1; point < points.size(); ++point) {
    const double z0 = points[point - 1].first;
    const double x0 = points[point - 1].second;
    const double z1 = points[point].first;
    const double x1 = points[point].second;
    const int nSteps = std::max(
        1, static_cast<int>(std::ceil(std::hypot(z1 - z0, x1 - x0) / sampleStep)));
    for (int step = 0; step <= nSteps; ++step) {
      const double fraction = static_cast<double>(step) / nSteps;
      const double z = z0 + fraction * (z1 - z0);
      const double x = x0 + fraction * (x1 - x0);
      const int bin = padMap->FindBin(z, x);
      if (bin > 0) {
        filledBins.insert(bin);
      }
    }
  }
  for (const int bin : filledBins) {
    padMap->SetBinContent(bin, padMap->GetBinContent(bin) + 1.0);
  }
  return filledBins.size();
}
} // namespace

// Example:
//   root -l -q 'rk-hs.cc+(2448)'
// The input file is <inputDir>/run%05d_K18Tracking.root.
void rk_hs(int runNumber = 2448, const char *inputDir = kDefaultInputDir,
           const char *outDir = "result")
{
  TString directory(inputDir);
  gSystem->ExpandPathName(directory);
  const TString inputFile = Form("%s/run%05d_K18Tracking.root", directory.Data(),
                                 runNumber);

  TFile file(inputFile, "READ");
  if (file.IsZombie()) {
    std::cerr << "Cannot open input file: " << inputFile << std::endl;
    return;
  }

  auto *tree = dynamic_cast<TTree *>(file.Get("k18"));
  if (!tree) {
    std::cerr << "Cannot find TTree 'k18' in " << inputFile << std::endl;
    return;
  }
  TBranch *xBranch = tree->GetBranch("xvpHS");
  TBranch *zBranch = tree->GetBranch("zvpHS");
  if (!xBranch || !zBranch) {
    std::cerr << "Branches xvpHS and/or zvpHS are missing in " << inputFile
              << std::endl;
    return;
  }

  const std::string xType = xBranch->GetClassName();
  const std::string zType = zBranch->GetClassName();
  const bool nested = xType == zType &&
                      (xType.find("vector<vector") != std::string::npos ||
                       xType.find("vector<std::vector") != std::string::npos);
  if (!nested) {
    std::cerr << "Expected xvpHS/zvpHS as matching nested-vector branches; got "
              << xType << " and " << zType << std::endl;
    return;
  }

  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir, kTRUE);
  auto *h = tpcdisp::MakeTPCPadMap(
      "h_xvpHS_vs_zvpHS", Form("run%05d: VPHS-interpolated K18 tracks", runNumber));
  Long64_t nPadFills = 0;
  const Long64_t nTracks = VisitNestedTracks(
      tree, [&](const auto &zValues, const auto &xValues) {
        nPadFills += FillInterpolatedTPCPads(zValues, xValues, h);
      });

  TCanvas canvas("c_xvpHS_vs_zvpHS", "VPHS-interpolated K18 tracks", 1000, 800);
  canvas.SetRightMargin(0.14);
  h->Draw("COLZ");

  const TString outputFile =
      Form("%s/run%05d_xvpHS_vs_zvpHS.pdf", outDir, runNumber);
  canvas.SaveAs(outputFile);
  std::cout << "Interpolated " << nTracks << " VPHS tracks into " << nPadFills
            << " TPC pads from " << inputFile << '\n'
            << "Wrote " << outputFile << std::endl;
}
