#include "DstTPCBranches.hh"
#include "TPC2DFrameBuilder.hh"

#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2Poly.h>
#include <TAxis.h>
#include <TAttLine.h>
#include <TKey.h>
#include <TObjString.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Options {
  std::string input;
  std::string output = "web/display.root";
  std::string treeName = "tpc";
  int firstEvent = 0;
  int count = 1;
};

struct EventSummary {
  struct TrackInfo {
    int pid = 0;
    int charge = 0;
    int isBeam = -1;
    int isAccidental = -1;
    double mom0 = 0.0;
    std::string particle;
  };

  struct VertexInfo {
    std::size_t row = 0;
    std::size_t col = 0;
    int track1 = -1;
    int track2 = -1;
  };

  int entry = 0;
  Long64_t entries = 0;
  std::size_t rawHits = 0;
  std::size_t clusters = 0;
  std::size_t tracks = 0;
  std::size_t trackClusters = 0;
  std::size_t vertexPairs = 0;
  std::vector<TrackInfo> trackInfo;
  std::vector<VertexInfo> vertexInfo;
};

void PrintUsage(const char *argv0) {
  std::cerr << "Usage: " << argv0
            << " [-i input.root] [-o web/display.root] [-t tree_name]"
            << " [-e first_entry] [-n count|all]\n";
}

bool ParseOptions(int argc, char **argv, Options &opt) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
      opt.input = argv[++i];
    } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
      opt.output = argv[++i];
    } else if ((arg == "-t" || arg == "--tree") && i + 1 < argc) {
      opt.treeName = argv[++i];
    } else if ((arg == "-e" || arg == "--event") && i + 1 < argc) {
      opt.firstEvent = std::atoi(argv[++i]);
    } else if ((arg == "-n" || arg == "--count") && i + 1 < argc) {
      const std::string value = argv[++i];
      opt.count = (value == "all") ? -1 : std::atoi(value.c_str());
    } else if (arg == "-h" || arg == "--help") {
      PrintUsage(argv[0]);
      return false;
    } else {
      std::cerr << "Unknown or incomplete option: " << arg << "\n";
      PrintUsage(argv[0]);
      return false;
    }
  }
  return true;
}

TTree *FindTree(TDirectory *dir, const std::string &treeName) {
  if (!treeName.empty()) {
    if (auto *tree = dynamic_cast<TTree *>(dir->Get(treeName.c_str()))) {
      return tree;
    }
  }

  TIter next(dir->GetListOfKeys());
  while (TKey *key = dynamic_cast<TKey *>(next())) {
    std::unique_ptr<TObject> obj(key->ReadObj());
    if (auto *tree = dynamic_cast<TTree *>(obj.get())) {
      obj.release();
      return tree;
    }
  }
  return nullptr;
}

void StyleAxis(TAxis *axis) {
  if (!axis) {
    return;
  }
  axis->SetTitleFont(132);
  axis->SetLabelFont(132);
}

void Style2DHistogram(TH1 *hist) {
  if (!hist) {
    return;
  }
  hist->SetStats(false);
  hist->SetTitle("");
  StyleAxis(hist->GetXaxis());
  StyleAxis(hist->GetYaxis());
}

void Style2DGraph(TGraph &graph) {
  graph.SetTitle(";Z [mm];X [mm]");
  StyleAxis(graph.GetXaxis());
  StyleAxis(graph.GetYaxis());
}

std::string DecodePid(int pid, int charge) {
  std::vector<std::string> labels;
  const char *sign = charge > 0 ? "+" : (charge < 0 ? "-" : "");

  if (pid & 0x1) {
    labels.emplace_back(std::string("pi") + sign);
  }
  if (pid & 0x2) {
    labels.emplace_back(std::string("K") + sign);
  }
  if (pid & 0x4) {
    labels.emplace_back(std::string("p") + sign);
  }

  if (labels.empty()) {
    return "unknown";
  }

  std::ostringstream out;
  for (std::size_t i = 0; i < labels.size(); ++i) {
    if (i != 0) {
      out << "/";
    }
    out << labels[i];
  }
  return out.str();
}

std::string TrackInfoJson(const std::vector<EventSummary::TrackInfo> &tracks) {
  std::ostringstream out;
  out << "[";
  for (std::size_t i = 0; i < tracks.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << "{\"pid\":" << tracks[i].pid
        << ",\"charge\":" << tracks[i].charge
        << ",\"is_beam\":" << tracks[i].isBeam
        << ",\"is_accidental\":" << tracks[i].isAccidental
        << ",\"particle\":\"" << tracks[i].particle << "\""
        << ",\"mom0\":" << std::fixed << std::setprecision(4) << tracks[i].mom0
        << "}";
  }
  out << "]";
  return out.str();
}

std::string VertexInfoJson(const std::vector<EventSummary::VertexInfo> &vertices) {
  std::ostringstream out;
  out << "[";
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << "{\"index\":" << i
        << ",\"row\":" << vertices[i].row
        << ",\"col\":" << vertices[i].col
        << ",\"track1\":" << vertices[i].track1
        << ",\"track2\":" << vertices[i].track2
        << "}";
  }
  out << "]";
  return out.str();
}

int WriteLine3D(TDirectory *out, const char *name, const std::vector<double> &xs,
                const std::vector<double> &ys, const std::vector<double> &zs,
                Color_t color = kGray + 2, Width_t width = 2) {
  const std::size_t npoints = std::min({xs.size(), ys.size(), zs.size()});
  if (npoints < 2) {
    return 0;
  }

  TPolyLine3D line(static_cast<Int_t>(npoints));
  line.SetLineColor(color);
  line.SetLineWidth(width);
  for (std::size_t i = 0; i < npoints; ++i) {
    line.SetPoint(static_cast<Int_t>(i), xs[i], ys[i], zs[i]);
  }

  out->cd();
  line.Write(name);
  return 1;
}

int WriteTPCWireframe(TDirectory *out) {
  constexpr double flength = 586.0;
  constexpr double fheight = 550.0;
  constexpr double pi = 3.14159265358979323846;
  const double edge = flength / (1.0 + std::sqrt(2.0));
  const double tan22 = std::tan(22.5 * pi / 180.0);

  const std::vector<double> x = {
      -flength / 2.0, -edge / 2.0, edge / 2.0, flength / 2.0,
      flength / 2.0, edge / 2.0, -edge / 2.0, -flength / 2.0,
      -flength / 2.0};
  const std::vector<double> z = {
      -tan22 * flength / 2.0, -flength / 2.0, -flength / 2.0,
      -tan22 * flength / 2.0, tan22 * flength / 2.0, flength / 2.0,
      flength / 2.0, tan22 * flength / 2.0, -tan22 * flength / 2.0};
  const std::vector<double> yLow(x.size(), -fheight / 2.0);
  const std::vector<double> yHigh(x.size(), fheight / 2.0);

  int count = 0;
  count += WriteLine3D(out, Form("tpc_frame_%d", count), x, yLow, z);
  count += WriteLine3D(out, Form("tpc_frame_%d", count), x, yHigh, z);
  for (std::size_t i = 0; i + 1 < x.size(); ++i) {
    count += WriteLine3D(out, Form("tpc_frame_%d", count),
                         {x[i], x[i]}, {-fheight / 2.0, fheight / 2.0},
                         {z[i], z[i]});
  }

  constexpr int circlePoints = 80;
  constexpr double targetRadius = 40.0;
  constexpr double targetCenterZ = -143.0;
  constexpr double targetHalfLength = 50.0;
  for (const double targetY : {-targetHalfLength, targetHalfLength}) {
    std::vector<double> tx;
    std::vector<double> ty;
    std::vector<double> tz;
    tx.reserve(circlePoints + 1);
    ty.reserve(circlePoints + 1);
    tz.reserve(circlePoints + 1);
    for (int i = 0; i <= circlePoints; ++i) {
      const double phi = 2.0 * pi * i / circlePoints;
      tx.push_back(targetRadius * std::cos(phi));
      ty.push_back(targetY);
      tz.push_back(targetCenterZ + targetRadius * std::sin(phi));
    }
    count += WriteLine3D(out, Form("tpc_frame_%d", count), tx, ty, tz, kAzure + 2, 2);
  }

  for (const double phi : {0.0, pi / 2.0, pi, 3.0 * pi / 2.0}) {
    const double tx = targetRadius * std::cos(phi);
    const double tz = targetCenterZ + targetRadius * std::sin(phi);
    count += WriteLine3D(out, Form("tpc_frame_%d", count), {tx, tx},
                         {-targetHalfLength, targetHalfLength},
                         {tz, tz}, kAzure + 2, 2);
  }

  constexpr double axisOriginX = 350.0;
  constexpr double axisOriginY = -260.0;
  constexpr double axisOriginZ = -300.0;
  constexpr double axisLength = 90.0;
  count += WriteLine3D(out, Form("tpc_frame_%d", count),
                       {axisOriginX, axisOriginX + axisLength},
                       {axisOriginY, axisOriginY},
                       {axisOriginZ, axisOriginZ}, kRed + 1, 3);
  count += WriteLine3D(out, Form("tpc_frame_%d", count),
                       {axisOriginX, axisOriginX},
                       {axisOriginY, axisOriginY + axisLength},
                       {axisOriginZ, axisOriginZ}, kGreen + 2, 3);
  count += WriteLine3D(out, Form("tpc_frame_%d", count),
                       {axisOriginX, axisOriginX},
                       {axisOriginY, axisOriginY},
                       {axisOriginZ, axisOriginZ + axisLength}, kBlue + 1, 3);
  return count;
}

void WriteGeometry(TDirectory *out) {
  out->cd();

  TPC2DFrameBuilder frame2d;
  std::unique_ptr<TH2Poly> padPlane(frame2d.TPC2DGeometry());
  padPlane->SetName("pad_plane_2d");
  Style2DHistogram(padPlane.get());
  padPlane->GetXaxis()->SetTitle("Z [mm]");
  padPlane->GetYaxis()->SetTitle("X [mm]");
  padPlane->SetLineColor(kWhite);
  padPlane->SetLineWidth(0);
  if (auto *bins = padPlane->GetBins()) {
    TIter next(bins);
    while (auto *obj = next()) {
      auto *bin = dynamic_cast<TH2PolyBin *>(obj);
      auto *line = bin ? dynamic_cast<TAttLine *>(bin->GetPolygon()) : nullptr;
      if (line) {
        line->SetLineColor(kWhite);
        line->SetLineWidth(0);
      }
    }
  }
  padPlane->Write();
  const int wireframeLines = WriteTPCWireframe(out);
  auto metadata = std::make_unique<TObjString>(Form("{\"wireframe_lines\":%d}", wireframeLines));
  metadata->Write("geometry_metadata");
}

std::size_t FlatSize(const std::vector<double> *xs,
                     const std::vector<double> *ys,
                     const std::vector<double> *zs) {
  if (!xs || !ys || !zs) {
    return 0;
  }
  return std::min({xs->size(), ys->size(), zs->size()});
}

std::size_t NestedTracks(const std::vector<std::vector<double>> *xs,
                         const std::vector<std::vector<double>> *ys,
                         const std::vector<std::vector<double>> *zs) {
  if (!xs || !ys || !zs) {
    return 0;
  }
  return std::min({xs->size(), ys->size(), zs->size()});
}

void WriteFlatPointSet(TDirectory *dir, const char *name, const char *,
                       const std::vector<double> *xs,
                       const std::vector<double> *ys,
                       const std::vector<double> *zs,
                       Color_t color, Style_t markerStyle, Size_t markerSize) {
  const std::size_t npoints = FlatSize(xs, ys, zs);
  if (npoints == 0) {
    return;
  }

  TPolyMarker3D hits3d(static_cast<Int_t>(npoints));
  hits3d.SetMarkerColor(color);
  hits3d.SetMarkerStyle(markerStyle);
  hits3d.SetMarkerSize(markerSize);

  TGraph hits2d(static_cast<Int_t>(npoints));
  hits2d.SetName(Form("%s2d", name));
  hits2d.SetMarkerColor(color);
  hits2d.SetMarkerStyle(markerStyle);
  hits2d.SetMarkerSize(markerSize);

  for (std::size_t i = 0; i < npoints; ++i) {
    hits3d.SetPoint(static_cast<Int_t>(i), xs->at(i), ys->at(i), zs->at(i));
    hits2d.SetPoint(static_cast<Int_t>(i), zs->at(i), xs->at(i));
  }
  Style2DGraph(hits2d);

  dir->cd();
  hits3d.Write(Form("%s3d", name));
  hits2d.Write(Form("%s2d", name));
}

void WriteTrackPointSet(TDirectory *dir, const char *prefix, const char *,
                        std::size_t track, const std::vector<double> &xs,
                        const std::vector<double> &ys, const std::vector<double> &zs,
                        Color_t color, Style_t markerStyle, Size_t markerSize) {
  const std::size_t npoints = std::min({xs.size(), ys.size(), zs.size()});
  if (npoints == 0) {
    return;
  }

  TPolyMarker3D hits3d(static_cast<Int_t>(npoints));
  hits3d.SetMarkerColor(color);
  hits3d.SetMarkerStyle(markerStyle);
  hits3d.SetMarkerSize(markerSize);

  TGraph hits2d(static_cast<Int_t>(npoints));
  hits2d.SetName(Form("track_%zu_%s2d", track, prefix));
  hits2d.SetMarkerColor(color);
  hits2d.SetMarkerStyle(markerStyle);
  hits2d.SetMarkerSize(markerSize);

  for (std::size_t i = 0; i < npoints; ++i) {
    hits3d.SetPoint(static_cast<Int_t>(i), xs[i], ys[i], zs[i]);
    hits2d.SetPoint(static_cast<Int_t>(i), zs[i], xs[i]);
  }
  Style2DGraph(hits2d);

  dir->cd();
  hits3d.Write(Form("track_%zu_%s3d", track, prefix));
  hits2d.Write(Form("track_%zu_%s2d", track, prefix));
}

void WriteTrackLine(TDirectory *dir, std::size_t track, const std::vector<double> &xs,
                    const std::vector<double> &ys, const std::vector<double> &zs,
                    Color_t color) {
  const std::size_t npoints = std::min({xs.size(), ys.size(), zs.size()});
  if (npoints < 2) {
    return;
  }

  TPolyLine3D line(static_cast<Int_t>(npoints));
  line.SetLineColor(color);
  line.SetLineWidth(2);

  TGraph line2d(static_cast<Int_t>(npoints));
  line2d.SetName(Form("track_%zu_line2d", track));
  line2d.SetLineColor(color);
  line2d.SetLineWidth(2);
  line2d.SetMarkerSize(0.0);
  for (std::size_t i = 0; i < npoints; ++i) {
    line.SetPoint(static_cast<Int_t>(i), xs[i], ys[i], zs[i]);
    line2d.SetPoint(static_cast<Int_t>(i), zs[i], xs[i]);
  }
  Style2DGraph(line2d);

  dir->cd();
  line.Write(Form("track_%zu_line3d", track));
  line2d.Write(Form("track_%zu_line2d", track));
}

bool WriteHelixLine(TDirectory *dir, std::size_t track, const DstTPCBranches &branches,
                    Color_t color) {
  const auto &tracks = branches.tracks;
  const auto &hits = branches.track_hits;
  if (!tracks.helix_cx || !tracks.helix_cy || !tracks.helix_z0 ||
      !tracks.helix_r || !tracks.helix_dz || !hits.helix_t ||
      track >= tracks.helix_cx->size() || track >= tracks.helix_cy->size() ||
      track >= tracks.helix_z0->size() || track >= tracks.helix_r->size() ||
      track >= tracks.helix_dz->size() || track >= hits.helix_t->size() ||
      hits.helix_t->at(track).empty()) {
    return false;
  }

  const auto &ts = hits.helix_t->at(track);
  const auto [minIt, maxIt] = std::minmax_element(ts.begin(), ts.end());
  const double tmin = *minIt;
  const double tmax = *maxIt;
  if (!std::isfinite(tmin) || !std::isfinite(tmax) || tmin == tmax) {
    return false;
  }

  constexpr int npoints = 160;
  TPolyLine3D line(npoints);
  line.SetLineColor(color);
  line.SetLineWidth(2);

  TGraph line2d(npoints);
  line2d.SetName(Form("track_%zu_line2d", track));
  line2d.SetLineColor(color);
  line2d.SetLineWidth(2);
  line2d.SetMarkerSize(0.0);

  const double cx = tracks.helix_cx->at(track);
  const double cy = tracks.helix_cy->at(track);
  const double z0 = tracks.helix_z0->at(track);
  const double r = tracks.helix_r->at(track);
  const double dz = tracks.helix_dz->at(track);
  for (int i = 0; i < npoints; ++i) {
    const double t = tmin + (tmax - tmin) * i / (npoints - 1);
    const double x = -(cx + r * std::cos(t));
    const double y = z0 + r * dz * t;
    const double z = cy + r * std::sin(t) - 143.0;
    line.SetPoint(i, x, y, z);
    line2d.SetPoint(i, z, x);
  }
  Style2DGraph(line2d);

  dir->cd();
  line.Write(Form("track_%zu_line3d", track));
  line2d.Write(Form("track_%zu_line2d", track));
  return true;
}

Color_t DefaultTrackColor(std::size_t track);

std::vector<EventSummary::VertexInfo> WriteVertexPairs(TDirectory *dir,
                                                       const DstTPCBranches &branches) {
  std::vector<EventSummary::VertexInfo> vertices;
  const auto &pairs = branches.pairs;
  if (!pairs.vtxTpc || !pairs.vtyTpc || !pairs.vtzTpc) {
    return vertices;
  }

  const std::size_t rows =
      std::min({pairs.vtxTpc->size(), pairs.vtyTpc->size(), pairs.vtzTpc->size()});
  std::size_t count = 0;
  for (std::size_t i = 0; i < rows; ++i) {
    const auto &xs = pairs.vtxTpc->at(i);
    const auto &ys = pairs.vtyTpc->at(i);
    const auto &zs = pairs.vtzTpc->at(i);
    const std::size_t cols = std::min({xs.size(), ys.size(), zs.size()});
    for (std::size_t j = i + 1; j < cols; ++j) {
      const double x = xs[j];
      const double y = ys[j];
      const double z = zs[j];
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
        continue;
      }

      EventSummary::VertexInfo info;
      info.row = i;
      info.col = j;
      info.track1 = static_cast<int>(i);
      info.track2 = static_cast<int>(j);
      if (pairs.combi_id && i < pairs.combi_id->size() && j < pairs.combi_id->at(i).size() &&
          std::isfinite(pairs.combi_id->at(i)[j])) {
        info.track2 = static_cast<int>(std::lround(pairs.combi_id->at(i)[j]));
      }

      const Color_t color = DefaultTrackColor(count);
      TPolyMarker3D marker3d(1);
      marker3d.SetMarkerColor(color);
      marker3d.SetMarkerStyle(34);
      marker3d.SetMarkerSize(1.2);
      marker3d.SetPoint(0, x, y, z);

      TGraph marker2d(1);
      marker2d.SetName(Form("vertex_pair_%zu_2d", count));
      marker2d.SetMarkerColor(color);
      marker2d.SetMarkerStyle(34);
      marker2d.SetMarkerSize(1.2);
      marker2d.SetPoint(0, z, x);
      Style2DGraph(marker2d);

      dir->cd();
      marker3d.Write(Form("vertex_pair_%zu_3d", count));
      marker2d.Write(Form("vertex_pair_%zu_2d", count));
      vertices.push_back(info);
      ++count;
    }
  }
  return vertices;
}

Color_t DefaultTrackColor(std::size_t track) {
  static const Color_t colors[] = {
      kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
      kOrange + 7, kCyan + 2, kViolet + 1, kGray + 2};
  return colors[track % (sizeof(colors) / sizeof(colors[0]))];
}

EventSummary WriteEvent(TTree *tree, int event, TDirectory *out) {
  EventSummary summary;
  summary.entry = event;
  summary.entries = tree ? tree->GetEntries() : 0;
  if (!tree) {
    return summary;
  }

  tree->ResetBranchAddresses();
  DstTPCBranches branches;
  branches.basic.SetBranchAddresses(tree, false);
  branches.raw_hits.SetBranchAddresses(tree, false);
  branches.clusters.SetBranchAddresses(tree, false);
  branches.tracks.SetBranchAddresses(tree, false);
  branches.track_hits.SetBranchAddresses(tree, false);
  branches.track_clusters.SetBranchAddresses(tree, false);
  branches.pairs.SetBranchAddresses(tree, false);
  branches.lambda.SetBranchAddresses(tree, false);

  const int selectedEvent = std::clamp(event, 0, static_cast<int>(tree->GetEntries() - 1));
  tree->GetEntry(selectedEvent);

  out->cd();
  auto *eventDir = out->mkdir(Form("event%d", selectedEvent));
  eventDir->cd();

  summary.rawHits = FlatSize(branches.raw_hits.raw_hitpos_x,
                             branches.raw_hits.raw_hitpos_y,
                             branches.raw_hits.raw_hitpos_z);
  WriteFlatPointSet(eventDir, "raw_hits", "raw hits", branches.raw_hits.raw_hitpos_x,
                    branches.raw_hits.raw_hitpos_y, branches.raw_hits.raw_hitpos_z,
                    kGray + 2, 20, 0.45);

  summary.clusters = FlatSize(branches.clusters.cluster_x, branches.clusters.cluster_y,
                              branches.clusters.cluster_z);
  WriteFlatPointSet(eventDir, "clusters", "clusters", branches.clusters.cluster_x,
                    branches.clusters.cluster_y, branches.clusters.cluster_z,
                    kOrange + 7, 21, 0.65);

  summary.tracks = NestedTracks(branches.track_hits.hitpos_x, branches.track_hits.hitpos_y,
                                branches.track_hits.hitpos_z);
  summary.trackInfo.reserve(summary.tracks);
  for (std::size_t track = 0; track < summary.tracks; ++track) {
    EventSummary::TrackInfo info;
    if (branches.tracks.pid && track < branches.tracks.pid->size()) {
      info.pid = branches.tracks.pid->at(track);
    }
    if (branches.tracks.charge && track < branches.tracks.charge->size()) {
      info.charge = branches.tracks.charge->at(track);
    }
    if (branches.tracks.is_beam && track < branches.tracks.is_beam->size()) {
      info.isBeam = branches.tracks.is_beam->at(track);
    }
    if (branches.tracks.is_accidental && track < branches.tracks.is_accidental->size()) {
      info.isAccidental = branches.tracks.is_accidental->at(track);
    }
    if (branches.tracks.mom0 && track < branches.tracks.mom0->size()) {
      info.mom0 = branches.tracks.mom0->at(track);
    }
    info.particle = DecodePid(info.pid, info.charge);
    summary.trackInfo.push_back(info);

    const auto &xs = branches.track_hits.hitpos_x->at(track);
    const auto &ys = branches.track_hits.hitpos_y->at(track);
    const auto &zs = branches.track_hits.hitpos_z->at(track);
    const Color_t color = DefaultTrackColor(track);
    WriteTrackPointSet(eventDir, "hits", "hits", track, xs, ys, zs, color, 20, 0.55);
    if (!WriteHelixLine(eventDir, track, branches, color)) {
      WriteTrackLine(eventDir, track, xs, ys, zs, color);
    }
  }

  const std::size_t calTracks = NestedTracks(branches.track_hits.calpos_x,
                                             branches.track_hits.calpos_y,
                                             branches.track_hits.calpos_z);
  for (std::size_t track = 0; track < calTracks; ++track) {
    WriteTrackPointSet(eventDir, "cal", "calculated positions", track,
                       branches.track_hits.calpos_x->at(track),
                       branches.track_hits.calpos_y->at(track),
                       branches.track_hits.calpos_z->at(track),
                       DefaultTrackColor(track), 24, 0.45);
  }

  const std::size_t clusterTracks = NestedTracks(branches.track_clusters.track_cluster_x_center,
                                                branches.track_clusters.track_cluster_y_center,
                                                branches.track_clusters.track_cluster_z_center);
  summary.trackClusters = clusterTracks;
  for (std::size_t track = 0; track < clusterTracks; ++track) {
    WriteTrackPointSet(eventDir, "clusters", "track clusters", track,
                       branches.track_clusters.track_cluster_x_center->at(track),
                       branches.track_clusters.track_cluster_y_center->at(track),
                       branches.track_clusters.track_cluster_z_center->at(track),
                       DefaultTrackColor(track), 22, 0.7);
  }

  summary.vertexInfo = WriteVertexPairs(eventDir, branches);
  summary.vertexPairs = summary.vertexInfo.size();

  auto metadata = std::make_unique<TObjString>(
      Form("{\"entry\":%d,\"entries\":%lld,\"raw_hits\":%zu,\"clusters\":%zu,"
           "\"tracks\":%zu,\"track_clusters\":%zu,\"vertex_pairs\":%zu,"
           "\"vertex_info\":%s}",
           selectedEvent, static_cast<Long64_t>(summary.entries), summary.rawHits,
           summary.clusters, summary.tracks, summary.trackClusters, summary.vertexPairs,
           VertexInfoJson(summary.vertexInfo).c_str()));
  metadata->Write("metadata");
  return summary;
}

void WriteDisplayMetadata(TDirectory *out, TTree *tree, int firstEvent, int count,
                          int writtenEvents, const EventSummary &lastSummary) {
  out->cd();
  auto metadata = std::make_unique<TObjString>(
      Form("{\"source\":\"%s\",\"tree\":\"%s\",\"first_event\":%d,"
           "\"exported_events\":%d,\"requested_events\":%d,\"entries\":%lld,"
           "\"tracks\":%zu,\"raw_hits\":%zu,\"clusters\":%zu,"
           "\"track_clusters\":%zu,\"vertex_pairs\":%zu}",
           tree && tree->GetCurrentFile() ? tree->GetCurrentFile()->GetName() : "",
           tree ? tree->GetName() : "", firstEvent, writtenEvents, count,
           tree ? static_cast<Long64_t>(tree->GetEntries()) : 0,
           lastSummary.tracks, lastSummary.rawHits, lastSummary.clusters,
           lastSummary.trackClusters, lastSummary.vertexPairs));
  metadata->Write("display_metadata");
}

} // namespace

int main(int argc, char **argv) {
  Options opt;
  if (!ParseOptions(argc, argv, opt)) {
    return 1;
  }

  const std::string outputDir = gSystem->DirName(opt.output.c_str());
  if (!outputDir.empty() && outputDir != ".") {
    gSystem->mkdir(outputDir.c_str(), true);
  }

  TFile out(opt.output.c_str(), "RECREATE");
  if (out.IsZombie()) {
    std::cerr << "Cannot create output file: " << opt.output << "\n";
    return 1;
  }

  WriteGeometry(&out);

  bool wroteEvent = false;
  int writtenEvents = 0;
  EventSummary lastSummary;
  std::unique_ptr<TFile> input;
  if (!opt.input.empty()) {
    input.reset(TFile::Open(opt.input.c_str(), "READ"));
    if (!input || input->IsZombie()) {
      std::cerr << "Cannot open input file: " << opt.input << "\n";
      return 1;
    }
    TTree *tree = FindTree(input.get(), opt.treeName);
    if (!tree) {
      std::cerr << "No TTree found in " << opt.input << "\n";
      return 1;
    }

    const auto entries = tree->GetEntries();
    if (entries <= 0) {
      std::cerr << "Tree " << tree->GetName() << " has no entries.\n";
      return 1;
    }
    const int firstEvent = std::clamp(opt.firstEvent, 0, static_cast<int>(entries - 1));
    const int requestedEvents =
        opt.count < 0 ? static_cast<int>(entries) - firstEvent : std::max(opt.count, 1);
    const int lastEvent =
        std::min(static_cast<int>(entries), firstEvent + requestedEvents);

    for (int event = firstEvent; event < lastEvent; ++event) {
      lastSummary = WriteEvent(tree, event, &out);
      wroteEvent = true;
      ++writtenEvents;
    }
    WriteDisplayMetadata(&out, tree, firstEvent, requestedEvents, writtenEvents, lastSummary);
  }

  out.Write();
  out.Close();

  std::cout << "Wrote " << opt.output;
  if (!wroteEvent) {
    std::cout << " (geometry only)";
  } else {
    std::cout << " (" << writtenEvents << " event entries"
              << ", tracks=" << lastSummary.tracks
              << ", raw_hits=" << lastSummary.rawHits
              << ", clusters=" << lastSummary.clusters
              << ", track_clusters=" << lastSummary.trackClusters
              << ", vertex_pairs=" << lastSummary.vertexPairs
              << ", vertex_info=" << VertexInfoJson(lastSummary.vertexInfo)
              << ", track_info=" << TrackInfoJson(lastSummary.trackInfo) << ")";
  }
  std::cout << "\n";
  return 0;
}
