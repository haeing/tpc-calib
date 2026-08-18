#include "../TPCPadHelper.hh"
#include "../TPCEventDisplayHelper.hh"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

const int runnumber = 2682;

void y_offset_calib(const char* result_subdir = "physics-755")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const string result_dir = Form("result/%s", result_subdir);
  const string param_dir = Form("param_history/%s", result_subdir);
  gSystem->mkdir(result_dir.c_str(), kTRUE);
  gSystem->mkdir(param_dir.c_str(), kTRUE);

  string outpdf = Form("result/%s/y-offset-calib-run0%d_1.pdf", result_subdir, runnumber);
  string outroot = Form("result/%s/y-offset-calib-run0%d_1.root", result_subdir, runnumber);
  string param_in = Form("param_history/%s/TPCParam_e72_run0%d_1",result_subdir, runnumber);
  string param_out = Form("param_history/%s/TPCParam_e72_run0%d_2",result_subdir, runnumber);
  string summary_out = Form("result/%s/y-offset-fit-run0%d_2.txt", result_subdir, runnumber);

  vector<double> res_offset(NumOfLayersTPC, 0.);
  vector<double> res_slope(NumOfLayersTPC, 0.);
  vector<double> res_offset_err(NumOfLayersTPC, 0.);
  vector<double> res_slope_err(NumOfLayersTPC, 0.);
  vector<bool> fit_ok(NumOfLayersTPC, false);

  auto c1 = new TCanvas("c1", "c1");

  TPaveText* title = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  title->AddText("y-offset-calib.cc");
  title->AddText("TPC Y calibration using Residual Y = offset + slope*Y");
  title->AddText(Form("run0%d", runnumber));
  title->AddText(Form("result directory: %s", result_dir.c_str()));
  title->AddText(Form("input param : %s", param_in.c_str()));
  title->AddText(Form("output param: %s", param_out.c_str()));
  TDatime now;
  title->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",
                      now.GetYear(), now.GetMonth(), now.GetDay(),
                      now.GetHour(), now.GetMinute(), now.GetSecond()));
  title->Draw();
  c1->Print((outpdf + "(").c_str());

  TFile* file = new TFile(Form("/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/%s/run0%d_DstTPCK18HelixTracking.root",
                               result_subdir, runnumber));
  if(!file || file->IsZombie()){
    Error("y_offset_calib", "cannot open input ROOT file");
    return;
  }

  TH2D* hist[NumOfLayersTPC];
  TF1* f_off[NumOfLayersTPC];
  TFile* fout = new TFile(outroot.c_str(), "RECREATE");

  c1->Clear();
  c1->Divide(8, 4);
  for(int i = 0; i < NumOfLayersTPC; i++){
    hist[i] = (TH2D*)file->Get(Form("TPCTrk_ResY_vs_Y_Layer%02d", i));
    if(!hist[i]){
      Warning("y_offset_calib", "missing histogram: TPCTrk_ResY_vs_Y_Layer%02d", i);
      continue;
    }

    f_off[i] = new TF1(Form("f_resy_vs_y_layer%02d", i), "[0] + [1]*x", -100, 100);

    c1->cd(i + 1);
    hist[i]->GetYaxis()->SetRangeUser(-10, 10);
    hist[i]->Fit(f_off[i], "R");

    res_offset[i] = f_off[i]->GetParameter(0);
    res_slope[i] = f_off[i]->GetParameter(1);
    res_offset_err[i] = f_off[i]->GetParError(0);
    res_slope_err[i] = f_off[i]->GetParError(1);
    fit_ok[i] = (hist[i]->GetEntries() > 0 && f_off[i]->GetNDF() > 0);

    fout->cd();
    hist[i]->Write();
    f_off[i]->Write();
  }
  c1->Print(outpdf.c_str());

  ofstream summary(summary_out.c_str());
  summary << "# run layer offset offset_err slope slope_err fit_ok\n";
  for(int i = 0; i < NumOfLayersTPC; i++){
    summary << runnumber << " "
            << i << " "
            << setprecision(12) << res_offset[i] << " "
            << setprecision(12) << res_offset_err[i] << " "
            << setprecision(12) << res_slope[i] << " "
            << setprecision(12) << res_slope_err[i] << " "
            << fit_ok[i] << "\n";
  }
  summary.close();

  ifstream fin(param_in.c_str());
  if(!fin.is_open()){
    Error("y_offset_calib", "cannot open input parameter file: %s", param_in.c_str());
    fout->Close();
    c1->Print((outpdf + ")").c_str());
    return;
  }

  ofstream param(param_out.c_str());
  if(!param.is_open()){
    Error("y_offset_calib", "cannot open output parameter file: %s", param_out.c_str());
    fout->Close();
    c1->Print((outpdf + ")").c_str());
    return;
  }

  bool in_tpc_y = false;
  bool saw_tpc_y_data = false;
  int updated_rows = 0;
  string line;

  while(getline(fin, line)){
    if(line.find("# TPC Y") != string::npos){
      in_tpc_y = true;
      param << line << "\n";
      continue;
    }

    if(in_tpc_y){
      istringstream iss(line);
      int layer, row, aty;
      double p0, p1;
      if((iss >> layer >> row >> aty >> p0 >> p1)){
        saw_tpc_y_data = true;
        if(5 <= layer && layer < NumOfLayersTPC && aty == 2 && fit_ok[layer]){
          const double scale = 1.0 - res_slope[layer];
          const double new_p1 = p1 * scale;
          const double new_p0 = p0 + res_offset[layer] / new_p1;
          param << setw(6) << layer
                << setw(7) << row
                << setw(7) << aty
                << setw(18) << fixed << setprecision(8) << new_p0
                << setw(18) << fixed << setprecision(8) << new_p1
                << "\n";
          updated_rows++;
          continue;
        }
      }
      else if(saw_tpc_y_data && line.find('#') != string::npos){
        in_tpc_y = false;
      }
    }

    param << line << "\n";
  }

  fin.close();
  param.close();

  c1->Clear();
  TPaveText* done = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  done->AddText("TPC Y parameter update");
  done->AddText("ResidualY fit: offset + slope*Y");
  done->AddText("Applied: p1_new = (1 - slope) * p1");
  done->AddText("Applied: p0_new = p0 + offset / p1_new");
  done->AddText(Form("updated rows: %d", updated_rows));
  done->AddText(Form("output param: %s", param_out.c_str()));
  done->AddText(Form("fit summary : %s", summary_out.c_str()));
  done->Draw();
  c1->Print((outpdf + ")").c_str());

  fout->Close();
}
