#include "../TPCPadHelper.hh"
#include "../TPCEventDisplayHelper.hh"

namespace {

const double me = 0.5109989461; // MeV/c2
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

  return rho * K * ZoverA / beta2 *
    (0.5 * TMath::Log(2. * me * beta2 * gamma2
                      * Wmax * MeVToeV * MeVToeV / (I * I))
     - beta2
     - 0.5 * delta);
}

double PionMean(double poq, double conv)
{
  const double p_MeV = 1000. * TMath::Abs(poq);
  if (p_MeV <= 0.) return 0.;

  const double energy = TMath::Hypot(mpi, p_MeV);
  const double beta = p_MeV / energy;

  return conv * BetheP10Raw(mpi, beta);
}

double CalcPionSigma(double poq)
{
  const double sigma_par[5] = {3.94842, 0.0138502, -0.110281, 12.6065, -10.9347};
  const double abspoq = TMath::Abs(poq);

  return sigma_par[0]
       + sigma_par[1] * abspoq
       + sigma_par[2] * poq * poq
       + sigma_par[3] * TMath::Exp(sigma_par[4] * abspoq);
}

bool PassPionOneSigma(double poq, double dedx, double conv)
{
  const double mean = PionMean(poq, conv);
  const double sigma = CalcPionSigma(poq);
  return dedx >= mean - sigma && dedx <= mean + sigma;
}

bool PassMipPionCut(double poq, double dedx, double conv,
                    double pmin, double pmax, double dedx_min, double dedx_max)
{
  const double absp = TMath::Abs(poq);
  return absp >= pmin
      && absp < pmax
      && dedx >= dedx_min
      && dedx <= dedx_max
      && PassPionOneSigma(poq, dedx, conv);
}

TGraph* MakePionCurve(const char* name, double conv, double nsigma,
                      double xmin = -2.0, double xmax = 2.0, int n = 400)
{
  TGraph* g = new TGraph(n);
  g->SetName(name);

  for (int i = 0; i < n; ++i) {
    const double poq = xmin + (xmax - xmin) * i / (n - 1);
    g->SetPoint(i, poq, PionMean(poq, conv) + nsigma * CalcPionSigma(poq));
  }

  return g;
}

TH2D* MakeMipPionHist(TH2D* src, const char* name, double conv,
                     double pmin, double pmax, double dedx_min, double dedx_max)
{
  TH2D* h = (TH2D*)src->Clone(name);
  h->SetDirectory(nullptr);
  h->SetTitle("MIP pion PID cut;#it{p/z} [GeV/#it{c}];TPC #LT#it{dE/dx}#GT (a.u.)");

  for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
    const double poq = h->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
      const double dedx = h->GetYaxis()->GetBinCenter(iy);
      if (!PassMipPionCut(poq, dedx, conv, pmin, pmax, dedx_min, dedx_max)) {
        h->SetBinContent(ix, iy, 0.);
        h->SetBinError(ix, iy, 0.);
      }
    }
  }

  return h;
}

void DrawPionOneSigmaCurves(double conv, double xmin, double xmax)
{
  TGraph* g_pi_mean = MakePionCurve("g_pi_mean", conv, 0.0, xmin, xmax);
  TGraph* g_pi_low = MakePionCurve("g_pi_low_1sigma", conv, -1.0, xmin, xmax);
  TGraph* g_pi_high = MakePionCurve("g_pi_high_1sigma", conv, 1.0, xmin, xmax);

  g_pi_mean->SetLineColor(kBlue);
  g_pi_mean->SetLineStyle(1);
  g_pi_mean->SetLineWidth(2);
  g_pi_low->SetLineColor(kBlue);
  g_pi_low->SetLineStyle(3);
  g_pi_low->SetLineWidth(2);
  g_pi_high->SetLineColor(kBlue);
  g_pi_high->SetLineStyle(3);
  g_pi_high->SetLineWidth(2);

  g_pi_mean->Draw("L same");
  g_pi_low->Draw("L same");
  g_pi_high->Draw("L same");
}

void DrawMipPionCutLines(double pmin, double pmax, double dedx_min, double dedx_max,
                         double ymin, double ymax)
{
  const double xs[4] = {-pmax, -pmin, pmin, pmax};
  for (int i = 0; i < 4; ++i) {
    TLine* line = new TLine(xs[i], ymin, xs[i], ymax);
    line->SetLineColor(kMagenta+2);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    line->Draw("same");
  }

  TLine* low = new TLine(-pmax, dedx_min, pmax, dedx_min);
  TLine* high = new TLine(-pmax, dedx_max, pmax, dedx_max);
  low->SetLineColor(kMagenta+2);
  high->SetLineColor(kMagenta+2);
  low->SetLineStyle(7);
  high->SetLineStyle(7);
  low->SetLineWidth(2);
  high->SetLineWidth(2);
  low->Draw("same");
  high->Draw("same");
}

}

void physics_padgain(int runnumber, Long64_t max_entries = -1){
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  string dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/physics-735";
  
  TFile *file = new TFile(Form("%s/run0%d_DstTPCHelixTracking.root",dir.c_str(),runnumber));
  if (!file || file->IsZombie()) {
    cout << "Cannot open input file for run0" << runnumber << endl;
    return;
  }
  
  TTree *tree = (TTree*)file->Get("tpc");
  if (!tree) {
    cout << "Cannot find tree: tpc" << endl;
    return;
  }

  string dir_save = "260701";
  gSystem->mkdir(Form("result/%s", dir_save.c_str()), kTRUE);
  string outpdf = Form("result/%s/physics-padgain-run0%d.pdf", dir_save.c_str(),runnumber);
  string outroot = Form("result/%s/physics-padgain-run0%d.root", dir_save.c_str(),runnumber);
  const double pion_conversion_factor = 7388.11;
  const double pid_xmin = -2.0;
  const double pid_xmax = 2.0;
  const double pid_ymin = 0.0;
  const double pid_ymax = 400.0;
  const double pion_mip_pmin = 0.20;
  const double pion_mip_pmax = 1.00;
  const double pion_mip_dedx_min = 10.0;
  const double pion_mip_dedx_max = 35.0;

  vector<double>* mom0 = nullptr;
  vector<double>* dEdx = nullptr;
  vector<int>* charge = nullptr;

  int ntTpc;
  vector<int>* nhtrack = nullptr;
  vector<vector<double>>* hitlayer = nullptr;
  vector<vector<double>>* theta_diff = nullptr;
  vector<vector<double>>* track_cluster_de = nullptr;
  vector<vector<double>>* track_cluster_y_center = nullptr;
  vector<vector<double>>* track_cluster_size= nullptr;
  vector<vector<double>>* track_cluster_mrow = nullptr;
  vector<vector<double>>* track_cluster_row_center = nullptr;

  tree->SetBranchAddress("ntTpc",&ntTpc);
  tree->SetBranchAddress("mom0",&mom0);
  tree->SetBranchAddress("dEdx",&dEdx);
  tree->SetBranchAddress("charge",&charge);
  tree->SetBranchAddress("nhtrack",&nhtrack);

  tree->SetBranchAddress("hitlayer",&hitlayer);
  tree->SetBranchAddress("theta_diff",&theta_diff);
  tree->SetBranchAddress("track_cluster_de",&track_cluster_de);
  tree->SetBranchAddress("track_cluster_y_center",&track_cluster_y_center);
  tree->SetBranchAddress("track_cluster_size",&track_cluster_size);
  tree->SetBranchAddress("track_cluster_mrow",&track_cluster_mrow);
  tree->SetBranchAddress("track_cluster_row_center",&track_cluster_row_center);

  TH2D *h_pid = (TH2D*)file->Get("PID_dEdx_vs_SignedMom");
  if (!h_pid) {
    cout << "Cannot find histogram: PID_dEdx_vs_SignedMom" << endl;
    return;
  }
  h_pid->SetDirectory(nullptr);
  TH2D* h_pion_1sigma = MakeMipPionHist(h_pid, "h_mip_pion_cut", pion_conversion_factor,
                                        pion_mip_pmin, pion_mip_pmax,
                                        pion_mip_dedx_min, pion_mip_dedx_max);

  auto TPC_hit = tpcdisp::MakeTPCPadMap("TPC_hit","TPC Hit Count");
  auto TPC_gain = tpcdisp::MakeTPCPadMap("TPC_gain","TPC Gain");
  TGraph *graph_gain = new TGraph();
  graph_gain->SetName("graph_gain");
  TH1D *hist_de[NumOfPadTPC];
  TH1D *hist_size = new TH1D("hist_size","hist_size;Cluster Size;Counts",5,-0.5,4.5);
  for(int i=0;i<NumOfPadTPC;i++){
    hist_de[i] = new TH1D(Form("hist_de%d",i),Form("hist_de%d;cluster dE [A.U.];Counts",i),100,0,1000);
  }

  bool y_cut = false;
  bool alpha_cut = false;
  bool cluster_cut = false;
  bool pid_cut = false;

  const Long64_t entries = max_entries > 0 ? TMath::Min(max_entries, tree->GetEntries()) : tree->GetEntries();
  for(Long64_t n = 0;n<entries;n++){
    tree->GetEntry(n);

    if(n%1000 == 0)cout<<n<<endl;

    for(int i=0;i<ntTpc;i++){
      const double poq = (*mom0)[i] * (*charge)[i];
      pid_cut = PassMipPionCut(poq, (*dEdx)[i], pion_conversion_factor,
                                pion_mip_pmin, pion_mip_pmax,
                                pion_mip_dedx_min, pion_mip_dedx_max);

      for(int j=0;j<(*nhtrack)[i];j++){
        y_cut = false;
        alpha_cut = false;
        cluster_cut = false;

        int padid = tpc::GetPadId((*hitlayer)[i][j],(*track_cluster_row_center)[i][j]);

        if((*track_cluster_y_center)[i][j]>-50 && (*track_cluster_y_center)[i][j]<50)y_cut = true;
        if((*hitlayer)[i][j]<10){
          if(TMath::Abs((*theta_diff)[i][j])<0.1)alpha_cut = true;
          else if(TMath::Abs((*theta_diff)[i][j])>TMath::Pi()-0.1 && TMath::Abs((*theta_diff)[i][j])<TMath::Pi()+0.1)alpha_cut = true;
        }
        else if((*hitlayer)[i][j]>=10){
          if(TMath::Abs((*theta_diff)[i][j])<0.2)alpha_cut = true;
        }
        if((*track_cluster_size)[i][j] == 1)cluster_cut = true;

        hist_size->Fill((*track_cluster_size)[i][j]);

        if(y_cut && alpha_cut && cluster_cut && pid_cut){
          double cnt = TPC_hit->GetBinContent(padid+1);
          TPC_hit->SetBinContent(padid+1,cnt+1);
          hist_de[padid]->Fill((*track_cluster_de)[i][j]);
        }
      }
    }
  }

  TFile *fout = new TFile(outroot.c_str(),"RECREATE");
  auto c1 = new TCanvas("c1","c1",900,700);

  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("physics-padgain.cc");
  p->AddText("TPC Pad Gain Calibration");
  p->AddText("MIP pion cut: mean #pm 1#sigma");
  p->AddText(Form("%.2f <= |p/z| < %.2f GeV/c", pion_mip_pmin, pion_mip_pmax));
  p->AddText(Form("%.1f <= dE/dx <= %.1f A.U.", pion_mip_dedx_min, pion_mip_dedx_max));
  p->AddText(Form("conversion factor = %.1f", pion_conversion_factor));
  p->AddText(Form("run0%d",runnumber));
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());

  c1->Clear();
  gPad->SetLogz();
  h_pid->GetXaxis()->SetRangeUser(pid_xmin, pid_xmax);
  h_pid->GetYaxis()->SetRangeUser(pid_ymin, pid_ymax);
  h_pid->GetXaxis()->SetTitle("#it{p/z} [GeV/#it{c}]");
  h_pid->GetYaxis()->SetTitle("TPC #LT#it{dE/dx}#GT (a.u.)");
  h_pid->Draw("colz");
  DrawPionOneSigmaCurves(pion_conversion_factor, pid_xmin, pid_xmax);
  DrawMipPionCutLines(pion_mip_pmin, pion_mip_pmax,
                      pion_mip_dedx_min, pion_mip_dedx_max, pid_ymin, pid_ymax);
  c1->Print(outpdf.c_str());

  c1->Clear();
  gPad->SetLogz();
  h_pion_1sigma->GetXaxis()->SetRangeUser(pid_xmin, pid_xmax);
  h_pion_1sigma->GetYaxis()->SetRangeUser(pid_ymin, pid_ymax);
  h_pion_1sigma->Draw("colz");
  DrawPionOneSigmaCurves(pion_conversion_factor, pid_xmin, pid_xmax);
  DrawMipPionCutLines(pion_mip_pmin, pion_mip_pmax,
                      pion_mip_dedx_min, pion_mip_dedx_max, pid_ymin, pid_ymax);
  h_pion_1sigma->Write();
  c1->Print(outpdf.c_str());

  c1->Clear();
  gPad->SetLogz();
  TPC_hit->Draw("colz");
  TPC_hit->Write();
  c1->Print(outpdf.c_str());

  c1->Clear();
  gPad->SetLogy(0);
  
  for(int ipad = 0;ipad <NumOfPadTPC;ipad++){
    TF1 fL("fL", "landau", 0, 1000);
    fL.SetParLimits(1,50,300);
    fL.SetLineColor(kRed);

    if(ipad%30 == 0){
      if(ipad !=0)c1->Print(outpdf.c_str());
      c1->Clear();
      c1->Divide(6,5);
    }
    c1->cd(ipad%30+1);
    if(hist_de[ipad]->GetEntries()>100){
      hist_de[ipad]->Fit(&fL, "QR");
      auto fit_param = fL.GetParameters();
      graph_gain->AddPoint(ipad,fit_param[1]);
      TPC_gain->SetBinContent(ipad+1,fit_param[1]);
    }

    hist_de[ipad]->Draw();
    hist_de[ipad]->Write();
  }

  c1->Print(outpdf.c_str());
  c1->Clear();

  hist_size->Draw();
  hist_size->Write();
  c1->Print(outpdf.c_str());

  c1->Clear();
  graph_gain->SetMarkerStyle(4);
  graph_gain->Draw("AP");
  graph_gain->Write();
  c1->Print(outpdf.c_str());

  c1->Clear();
  TPC_gain->SetMinimum(0);
  TPC_gain->SetMaximum(300);
  TPC_gain->Draw("colz");
  TPC_gain->Write();
  c1->Print((outpdf + ")").c_str());

  fout->Close();
}
