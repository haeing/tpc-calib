#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {
constexpr const char *kDefaultInputDir = "~/data/JPARC2025Nov_root/physics-735";
constexpr double kZTarget = -143.0; // tpc::Z_TARGET

template <class T>
bool Bind(TTree *tree, const char *name, T *address)
{
  if (!tree->GetBranch(name)) {
    std::cerr << "Missing branch: " << name << std::endl;
    return false;
  }
  tree->SetBranchAddress(name, address);
  return true;
}

// Same z-to-position choice as DstTPCK18HelixTracking.cc::PositionOnHelixAtZ.
bool PositionOnHelixAtZ(double cx, double cy, double z0, double radius,
                        double dz, double thetaMin, double thetaMax,
                        double z, double &x, double &y)
{
  if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(z0) ||
      !std::isfinite(radius) || !std::isfinite(dz) ||
      !std::isfinite(thetaMin) || !std::isfinite(thetaMax) ||
      !std::isfinite(z) || std::abs(radius) < 1.0e-9 || thetaMin > thetaMax) {
    return false;
  }

  const double sine = (z - kZTarget - cy) / radius;
  if (sine < -1.0 || sine > 1.0) {
    return false;
  }
  const double pi = std::acos(-1.0);
  const double base = std::asin(sine);
  const double roots[] = {base, pi - base};
  const double thetaCenter = 0.5 * (thetaMin + thetaMax);
  double thetaBest = std::numeric_limits<double>::quiet_NaN();
  double bestDistance = std::numeric_limits<double>::infinity();
  const double twoPi = 2.0 * pi;
  for (const double root : roots) {
    const int nMin = static_cast<int>(std::floor((thetaMin - root) / twoPi)) - 1;
    const int nMax = static_cast<int>(std::ceil((thetaMax - root) / twoPi)) + 1;
    for (int n = nMin; n <= nMax; ++n) {
      const double theta = root + twoPi * n;
      if (theta < thetaMin || theta > thetaMax) {
        continue;
      }
      const double distance = std::abs(theta - thetaCenter);
      if (distance < bestDistance) {
        bestDistance = distance;
        thetaBest = theta;
      }
    }
  }
  if (!std::isfinite(thetaBest)) {
    return false;
  }

  // Local helix: (cx + r cos(theta), cy + r sin(theta), z0 + dz*r*theta).
  // LocalToGlobal: (-localX, localZ, localY + Z_TARGET).
  x = -(cx + radius * std::cos(thetaBest));
  y = z0 + dz * radius * thetaBest;
  return std::isfinite(x) && std::isfinite(y);
}
} // namespace

// Residual convention: TPC helix position - K18 VP-helix position.
// Events are used only when one is_k18 TPC track and one K18 VP helix exist.
void k18_tpc_residual(int runNumber = 2448,
                      const char *inputDir = kDefaultInputDir,
                      const char *outDir = "result")
{
  TString directory(inputDir);
  gSystem->ExpandPathName(directory);
  const TString inputFile =
      Form("%s/run%05d_DstTPCK18HelixTracking.root", directory.Data(), runNumber);

  TFile file(inputFile, "READ");
  if (file.IsZombie()) {
    std::cerr << "Cannot open input file: " << inputFile << std::endl;
    return;
  }
  auto *tree = dynamic_cast<TTree *>(file.Get("tpc"));
  if (!tree) {
    std::cerr << "Cannot find TTree 'tpc' in " << inputFile << std::endl;
    return;
  }

  int ntBeam = 0;
  std::vector<int> *isBeam = nullptr;
  std::vector<double> *tpcCx = nullptr;
  std::vector<double> *tpcCy = nullptr;
  std::vector<double> *tpcZ0 = nullptr;
  std::vector<double> *tpcR = nullptr;
  std::vector<double> *tpcDz = nullptr;
  std::vector<double> *tpcThetaMin = nullptr;
  std::vector<double> *tpcThetaMax = nullptr;
  std::vector<double> *k18Cx = nullptr;
  std::vector<double> *k18Cy = nullptr;
  std::vector<double> *k18Z0 = nullptr;
  std::vector<double> *k18R = nullptr;
  std::vector<double> *k18Dz = nullptr;
  std::vector<double> *k18ThetaMin = nullptr;
  std::vector<double> *k18ThetaMax = nullptr;
  std::vector<std::vector<double>> *xCalBeam = nullptr;
  std::vector<std::vector<double>> *yCalBeam = nullptr;
  std::vector<std::vector<double>> *zCalBeam = nullptr;
  std::vector<double> *tpcMom0 = nullptr;
  std::vector<double> *pBeam = nullptr;

  bool ok = true;
  ok &= Bind(tree, "ntK18", &ntBeam);
  ok &= Bind(tree, "is_k18", &isBeam);
  ok &= Bind(tree, "helix_cx", &tpcCx);
  ok &= Bind(tree, "helix_cy", &tpcCy);
  ok &= Bind(tree, "helix_z0", &tpcZ0);
  ok &= Bind(tree, "helix_r", &tpcR);
  ok &= Bind(tree, "helix_dz", &tpcDz);
  ok &= Bind(tree, "helix_theta_min", &tpcThetaMin);
  ok &= Bind(tree, "helix_theta_max", &tpcThetaMax);
  ok &= Bind(tree, "helix_cx_k18", &k18Cx);
  ok &= Bind(tree, "helix_cy_k18", &k18Cy);
  ok &= Bind(tree, "helix_z0_k18", &k18Z0);
  ok &= Bind(tree, "helix_r_k18", &k18R);
  ok &= Bind(tree, "helix_dz_k18", &k18Dz);
  ok &= Bind(tree, "helix_theta_min_k18", &k18ThetaMin);
  ok &= Bind(tree, "helix_theta_max_k18", &k18ThetaMax);
  ok &= Bind(tree, "x_cal_k18", &xCalBeam);
  ok &= Bind(tree, "y_cal_k18", &yCalBeam);
  ok &= Bind(tree, "z_cal_k18", &zCalBeam);
  ok &= Bind(tree, "mom0", &tpcMom0);
  ok &= Bind(tree, "mom0_k18", &pBeam);
  if (!ok) {
    return;
  }

  gStyle->SetOptStat("emr");
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetStatX(0.88);
  gStyle->SetStatY(0.88);
  gStyle->SetStatW(0.20);
  gStyle->SetStatH(0.16);
  gSystem->mkdir(outDir, kTRUE);
  // Every parameter residual below is defined as (TPC helix) - (K18 VP helix).
  const char *parameterDescriptions[] = {
      "#Delta c_{x}: helix-center x coordinate",
      "#Delta c_{y}: helix-center local-y coordinate",
      "#Delta z_{0}: helix longitudinal offset",
      "#Delta R: helix radius",
      "#Delta dz: helix pitch parameter",
      "#Delta#theta_{mid}: mean helix angle"};
  const char *parameterAxes[] = {
      "#Delta c_{x} = c_{x}^{TPC}-c_{x}^{Beam} [mm]",
      "#Delta c_{y} = c_{y}^{TPC}-c_{y}^{Beam} [mm]",
      "#Delta z_{0} = z_{0}^{TPC}-z_{0}^{Beam} [mm]",
      "#Delta R = R^{TPC}-R^{Beam} [mm]",
      "#Delta dz = dz^{TPC}-dz^{Beam}",
      "#Delta#theta_{mid} = #theta_{mid}^{TPC}-#theta_{mid}^{Beam} [rad]"};
  TH1D *parameterResiduals[6];
  for (int i = 0; i < 6; ++i) {
    const double xRange = (i == 3 || i == 4) ? 1.0 : 20.0;
    parameterResiduals[i] = new TH1D(
        Form("hParameterResidual%d", i),
        Form("run%05d: %s (TPC - Beam);%s;Counts", runNumber,
             parameterDescriptions[i], parameterAxes[i]),
        200, -xRange, xRange);
  }

  TH1D hMomentumResidual(
      "hMomentumResidual",
      "Momentum residual: TPC mom0 - Beam mom0_{k18};#Delta p = p_{TPC} - p_{Beam} [GeT/#it{c}];Counts",
      200, -0.5, 0.5);
  TH2D hMom0Correlation(
      "hMom0Correlation",
      "TPC mom0 versus Beam mom0_{k18};p_{Beam} [GeT/#it{c}];p_{TPC} [GeT/#it{c}]",
      100, 0.6, 0.8, 100, 0.6, 0.8);

  TH2D hDeltaXVsZ("hDeltaXVsZ",
                  "TPC X - Beam X at each Beam z;Z [mm];#DeltaX [mm]",
                  50, -300.0, 300.0, 200, -20.0, 20.0);
  TH2D hDeltaYVsZ("hDeltaYVsZ",
                  "TPC Y - Beam Y at each Beam z;Z [mm];#DeltaY [mm]",
                  50, -300.0, 300.0, 200, -20.0, 20.0);
  TH2D hDistanceVsZ("hDistanceVsZ",
                    "Transverse distance between TPC and Beam helices;Z [mm];#sqrt{#DeltaX^{2}+#DeltaY^{2}} [mm]",
                    50, -300.0, 300.0, 100, 0.0, 20.0);
  TH1D hDeltaXAllZ("hDeltaXAllZ",
                   "TPC X - Beam X at all Beam z in TPC volume;#DeltaX [mm];Counts",
                   200, -20.0, 20.0);
  TH1D hDeltaYAllZ("hDeltaYAllZ",
                   "TPC Y - Beam Y at all Beam z in TPC volume;#DeltaY [mm];Counts",
                   200, -20.0, 20.0);
  TH1D hDistanceAllZ("hDistanceAllZ",
                     "Transverse distance at all Beam z in TPC volume;#sqrt{#DeltaX^{2}+#DeltaY^{2}} [mm];Counts",
                     100, 0.0, 20.0);

  Long64_t nUnique = 0;
  Long64_t nAmbiguous = 0;
  Long64_t nPositions = 0;
  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree->GetEntry(entry);
    if (!isBeam || !tpcCx || !tpcCy || !tpcZ0 || !tpcR || !tpcDz ||
        !tpcThetaMin || !tpcThetaMax || !k18Cx || !k18Cy || !k18Z0 ||
        !k18R || !k18Dz || !k18ThetaMin || !k18ThetaMax || !xCalBeam ||
        !yCalBeam || !zCalBeam) continue;

    std::vector<size_t> tagged;
    for (size_t i = 0; i < isBeam->size(); ++i)
      if (isBeam->at(i) == 1) tagged.push_back(i);
    if (tagged.size() != 1 || ntBeam != 1 || k18Cx->size() != 1 ||
        xCalBeam->size() != 1 || yCalBeam->size() != 1 || zCalBeam->size() != 1) {
      ++nAmbiguous;
      continue;
    }

    const size_t tpc = tagged.front();
    const size_t k18 = 0;
    if (tpc >= tpcCx->size() || tpc >= tpcCy->size() || tpc >= tpcZ0->size() ||
        tpc >= tpcR->size() || tpc >= tpcDz->size() ||
        tpc >= tpcThetaMin->size() || tpc >= tpcThetaMax->size() ||
        k18 >= k18Cy->size() || k18 >= k18Z0->size() || k18 >= k18R->size() ||
        k18 >= k18Dz->size() || k18 >= k18ThetaMin->size() ||
        k18 >= k18ThetaMax->size()) continue;

    const double parameterDifferences[] = {
        tpcCx->at(tpc) - k18Cx->at(k18), tpcCy->at(tpc) - k18Cy->at(k18),
        tpcZ0->at(tpc) - k18Z0->at(k18), tpcR->at(tpc) - k18R->at(k18),
        tpcDz->at(tpc) - k18Dz->at(k18),
        0.5 * (tpcThetaMin->at(tpc) + tpcThetaMax->at(tpc) -
               k18ThetaMin->at(k18) - k18ThetaMax->at(k18))};
    for (int i = 0; i < 6; ++i)
      if (std::isfinite(parameterDifferences[i])) parameterResiduals[i]->Fill(parameterDifferences[i]);
    ++nUnique;
    if (tpc < tpcMom0->size() && k18 < pBeam->size()) {
      const double tpcMomentum = tpcMom0->at(tpc);
      const double beamMomentum = pBeam->at(k18);
      const double deltaP = tpcMomentum - beamMomentum;
      if (std::isfinite(deltaP)) {
        hMomentumResidual.Fill(deltaP);
        hMom0Correlation.Fill(beamMomentum, tpcMomentum);
      }
    }

    const auto &refX = xCalBeam->at(k18);
    const auto &refY = yCalBeam->at(k18);
    const auto &refZ = zCalBeam->at(k18);
    const size_t nVP = std::min({refX.size(), refY.size(), refZ.size()});
    for (size_t vp = 0; vp < nVP; ++vp) {
      const double z = refZ[vp];
      if (!std::isfinite(z) || !std::isfinite(refX[vp]) ||
          !std::isfinite(refY[vp])) continue;
      double tpcX = 0.0, tpcY = 0.0;
      if (!PositionOnHelixAtZ(tpcCx->at(tpc), tpcCy->at(tpc), tpcZ0->at(tpc),
                              tpcR->at(tpc), tpcDz->at(tpc),
                              tpcThetaMin->at(tpc), tpcThetaMax->at(tpc), z,
                              tpcX, tpcY)) continue;
      const double deltaX = tpcX - refX[vp];
      const double deltaY = tpcY - refY[vp];
      const double distance = std::hypot(deltaX, deltaY);
      hDeltaXVsZ.Fill(z, deltaX);
      hDeltaYVsZ.Fill(z, deltaY);
      hDistanceVsZ.Fill(z, distance);
      hDeltaXAllZ.Fill(deltaX);
      hDeltaYAllZ.Fill(deltaY);
      hDistanceAllZ.Fill(distance);
      ++nPositions;
    }
  }

  const TString outputFile =
      Form("%s/run%05d_k18_tpc_helix_residual.pdf", outDir, runNumber);
  TCanvas canvas("c_k18_tpc_residual", "Beam-TPC helix residual", 1200, 800);
  for (int i = 0; i < 6; ++i) {
    canvas.Clear();
    parameterResiduals[i]->Draw();
    canvas.Print(i == 0 ? (outputFile + "(").Data() : outputFile.Data());
  }
  canvas.Clear();
  hMomentumResidual.Draw();
  canvas.Print(outputFile.Data());


  TH2D *zResiduals[] = {&hDeltaXVsZ, &hDeltaYVsZ, &hDistanceVsZ,
                         &hMom0Correlation};
  for (TH2D *histogram : zResiduals) {
    canvas.Clear();
    histogram->Draw("COLZ");
    canvas.Print(outputFile.Data());
  }

  TH1D *allZResiduals[] = {&hDeltaXAllZ, &hDeltaYAllZ, &hDistanceAllZ};
  for (int quantity = 0; quantity < 3; ++quantity) {
    canvas.Clear();
    allZResiduals[quantity]->Draw();
    canvas.Print(quantity == 2 ? (outputFile + ")").Data() : outputFile.Data());
  }

  std::cout << "Used " << nUnique << " unambiguous events and " << nPositions
            << " TPC-volume z positions; skipped " << nAmbiguous
            << " events with multiple/no Beam candidates.\nWrote " << outputFile
            << std::endl;
}
