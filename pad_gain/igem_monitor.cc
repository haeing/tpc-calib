void igem_monitor()
{
  TString filename = "../../../data/JPARC2025Nov/monitor-tmp/caenhv-hyptpc.txt";

  std::ifstream fin(filename.Data());
  if (!fin.is_open()) {
    std::cerr << "Cannot open file: " << filename << std::endl;
    return;
  }

  std::vector<double> t;
  std::vector<double> vgemSet;
  std::vector<double> vgemMon;
  std::vector<double> igemMon;

  double unix_time;
  double VCatSet, VCATMon, ICatMon;
  double VGEMSet, VGEMMon, IGEMMon;
  double VG0Set, VG0Mon, IG0Mon;
  double VGPSet, VGPMon, IGPMon;
  double VGMSet, VGMMon, IGMMon;

  while (fin >> unix_time
             >> VCatSet >> VCATMon >> ICatMon
             >> VGEMSet >> VGEMMon >> IGEMMon
             >> VG0Set  >> VG0Mon  >> IG0Mon
             >> VGPSet  >> VGPMon  >> IGPMon
             >> VGMSet  >> VGMMon  >> IGMMon) {

    t.push_back(unix_time);
    vgemSet.push_back(VGEMSet);
    vgemMon.push_back(VGEMMon);
    igemMon.push_back(IGEMMon);
  }

  fin.close();

  int n = t.size();
  if (n == 0) {
    std::cerr << "No data found." << std::endl;
    return;
  }

  auto grVGEMSet = new TGraph(n, t.data(), vgemSet.data());
  auto grVGEMMon = new TGraph(n, t.data(), vgemMon.data());
  auto grIGEMMon = new TGraph(n, t.data(), igemMon.data());

  grVGEMSet->SetName("grVGEMSet");
  grVGEMMon->SetName("grVGEMMon");
  grIGEMMon->SetName("grIGEMMon");

  grVGEMSet->SetTitle("VGEM voltage;Unix time;Voltage [V]");
  grVGEMMon->SetTitle("VGEM voltage;Unix time;Voltage [V]");

  grVGEMSet->SetMarkerStyle(20);
  grVGEMMon->SetMarkerStyle(21);

  grVGEMSet->SetLineColor(kRed);
  grVGEMMon->SetLineColor(kBlue);
  grVGEMSet->SetMarkerColor(kRed);
  grVGEMMon->SetMarkerColor(kBlue);

  auto c1 = new TCanvas("c1", "VGEM voltage", 1000, 700);
  grVGEMSet->Draw("APL");
  grVGEMMon->Draw("PL same");

  auto leg1 = new TLegend(0.70, 0.75, 0.90, 0.90);
  leg1->AddEntry(grVGEMSet, "VGEMSet", "lp");
  leg1->AddEntry(grVGEMMon, "VGEMMon", "lp");
  leg1->Draw();

  auto c2 = new TCanvas("c2", "IGEM current", 1000, 700);
  grIGEMMon->SetTitle("IGEM current;Unix time;Current");
  grIGEMMon->SetMarkerStyle(20);
  grIGEMMon->SetLineColor(kGreen + 2);
  grIGEMMon->SetMarkerColor(kGreen + 2);
  grIGEMMon->GetYaxis()->SetRangeUser(155, 170);
  grIGEMMon->Draw("APL");

  /*
  auto fout = new TFile("gem_hv_graph.root", "RECREATE");
  grVGEMSet->Write();
  grVGEMMon->Write();
  grIGEMMon->Write();
  fout->Close();
  */
}
