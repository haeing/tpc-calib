namespace {

struct K18Result {
  Double_t p = 0.;
  Double_t deltaPercent = 0.;
  Double_t chi2ndf = 0.;
  Int_t ndf = 0;
  Double_t blc1Chi2 = 0.;
  Double_t blc2Chi2 = 0.;
  Double_t xOut = 0.;
  Double_t yOut = 0.;
  Double_t uOut = 0.;
  Double_t vOut = 0.;
  Double_t xBcOut = 0.;
  Double_t yBcOut = 0.;
  Double_t uBcOut = 0.;
  Double_t vBcOut = 0.;
};

ULong64_t EventKey(UInt_t run, UInt_t event)
{
  return (static_cast<ULong64_t>(run) << 32) | event;
}

bool HasOne(const std::vector<Double_t>* value)
{
  return value && value->size() == 1;
}

} // namespace

void CompareD5K18SingleLocal()
{
  const char* outputName = "d5_k18_single_local.root";
  const char* d5Name = "~/data/JPARC2025Nov_root/blc/run02448_D5Tracking.root";
  const char* k18Name = "~/data/JPARC2025Nov_root/blc/run02448_K18Tracking.root";
  TFile d5File(d5Name, "READ");
  TFile k18File(k18Name, "READ");
  if (d5File.IsZombie() || k18File.IsZombie()) {
    std::cerr << "Cannot open input file(s)." << std::endl;
    return;
  }

  auto* d5 = dynamic_cast<TTree*>(d5File.Get("d5"));
  auto* k18 = dynamic_cast<TTree*>(k18File.Get("k18"));
  if (!d5 || !k18) {
    std::cerr << "Expected trees named d5 and k18." << std::endl;
    return;
  }

  UInt_t kRun = 0, kEvent = 0;
  Int_t ntBcIn = 0, ntBcOut = 0, ntK18 = 0;
  std::vector<Double_t>* pk18 = nullptr;
  std::vector<Double_t>* deltaK18 = nullptr;
  std::vector<Double_t>* chi2K18 = nullptr;
  std::vector<Int_t>* ndfK18 = nullptr;
  std::vector<Double_t>* chi2BcIn = nullptr;
  std::vector<Double_t>* chi2BcOut = nullptr;
  std::vector<Double_t>* xin = nullptr;
  std::vector<Double_t>* yin = nullptr;
  std::vector<Double_t>* uin = nullptr;
  std::vector<Double_t>* vin = nullptr;
  std::vector<Double_t>* xout = nullptr;
  std::vector<Double_t>* yout = nullptr;
  std::vector<Double_t>* uout = nullptr;
  std::vector<Double_t>* vout = nullptr;
  std::vector<Double_t>* x0BcOut = nullptr;
  std::vector<Double_t>* y0BcOut = nullptr;
  std::vector<Double_t>* u0BcOut = nullptr;
  std::vector<Double_t>* v0BcOut = nullptr;
  
  
  k18->SetBranchAddress("run_number", &kRun);
  k18->SetBranchAddress("event_number", &kEvent);
  k18->SetBranchAddress("ntBcIn", &ntBcIn);
  k18->SetBranchAddress("ntBcOut", &ntBcOut);
  k18->SetBranchAddress("ntK18", &ntK18);
  k18->SetBranchAddress("pk18", &pk18);
  k18->SetBranchAddress("delta", &deltaK18);
  k18->SetBranchAddress("chisqrK18", &chi2K18);
  k18->SetBranchAddress("ndfK18", &ndfK18);
  k18->SetBranchAddress("chisqrBcIn", &chi2BcIn);
  k18->SetBranchAddress("chisqrBcOut", &chi2BcOut);
  k18->SetBranchAddress("xout", &xout);
  k18->SetBranchAddress("yout", &yout);
  k18->SetBranchAddress("uout", &uout);
  k18->SetBranchAddress("vout", &vout);
  k18->SetBranchAddress("x0BcOut", &x0BcOut);
  k18->SetBranchAddress("y0BcOut", &y0BcOut);
  k18->SetBranchAddress("u0BcOut", &u0BcOut);
  k18->SetBranchAddress("v0BcOut", &v0BcOut);
  

  std::unordered_map<ULong64_t, K18Result> k18ByEvent;
  for (Long64_t entry = 0; entry < k18->GetEntries(); ++entry) {
    k18->GetEntry(entry);
    if (ntBcIn != 1 || ntBcOut != 1 || ntK18 != 1 || !HasOne(pk18) ||
        !HasOne(deltaK18) || !HasOne(chi2K18) || !ndfK18 || ndfK18->size() != 1 ||
        !HasOne(chi2BcIn) || !HasOne(chi2BcOut) ||
        !HasOne(xout) || !HasOne(yout) || !HasOne(uout) || !HasOne(vout) ||
        !HasOne(x0BcOut) || !HasOne(y0BcOut) || !HasOne(u0BcOut) || !HasOne(v0BcOut))
      continue;
    k18ByEvent.emplace(EventKey(kRun, kEvent),
                        K18Result{pk18->at(0), 100.*deltaK18->at(0), chi2K18->at(0), ndfK18->at(0),
                                  chi2BcIn->at(0), chi2BcOut->at(0),
                                  xout->at(0), yout->at(0), uout->at(0), vout->at(0),
                                  x0BcOut->at(0), y0BcOut->at(0), u0BcOut->at(0), v0BcOut->at(0)});
  }

  UInt_t dRun = 0, dEvent = 0;
  Int_t ntBLC1 = 0, ntBLC2 = 0, ntD5 = 0;
  std::vector<Double_t>* pD5 = nullptr;
  std::vector<Double_t>* deltaD5 = nullptr;
  std::vector<Double_t>* chi2D5 = nullptr;
  std::vector<Double_t>* chi2NdfD5 = nullptr;
  std::vector<Int_t>* ndfD5 = nullptr;
  std::vector<Int_t>* pairRank = nullptr;
  std::vector<Double_t>* blc1Chi2 = nullptr;
  std::vector<Double_t>* blc2Chi2 = nullptr;
  std::vector<Double_t>* blc2U = nullptr;
  std::vector<Double_t>* blc2V = nullptr;
  std::vector<Double_t>* d5BcOutX = nullptr;
  std::vector<Double_t>* d5BcOutY = nullptr;
  std::vector<Double_t>* d5MtxOutX = nullptr;
  std::vector<Double_t>* d5MtxOutU = nullptr;
  std::vector<Double_t>* d5MtxOutY = nullptr;
  std::vector<Double_t>* d5MtxOutV = nullptr;
  d5->SetBranchAddress("run_number", &dRun);
  d5->SetBranchAddress("event_number", &dEvent);
  d5->SetBranchAddress("ntrack_blc1", &ntBLC1);
  d5->SetBranchAddress("ntrack_blc2", &ntBLC2);
  d5->SetBranchAddress("ntrack_d5", &ntD5);
  d5->SetBranchAddress("d5_momentum", &pD5);
  d5->SetBranchAddress("d5_delta", &deltaD5);
  d5->SetBranchAddress("d5_fit_chi2", &chi2D5);
  d5->SetBranchAddress("d5_fit_chi2_ndf", &chi2NdfD5);
  d5->SetBranchAddress("d5_fit_ndf", &ndfD5);
  d5->SetBranchAddress("d5_pair_rank", &pairRank);
  d5->SetBranchAddress("blc1_chi2", &blc1Chi2);
  d5->SetBranchAddress("blc2_chi2", &blc2Chi2);
  d5->SetBranchAddress("blc2_u", &blc2U);
  d5->SetBranchAddress("blc2_v", &blc2V);
  d5->SetBranchAddress("d5_blc2out_x", &d5BcOutX);
  d5->SetBranchAddress("d5_blc2out_y", &d5BcOutY);
  d5->SetBranchAddress("d5_mtxout_x", &d5MtxOutX);
  d5->SetBranchAddress("d5_mtxout_u", &d5MtxOutU);
  d5->SetBranchAddress("d5_mtxout_y", &d5MtxOutY);
  d5->SetBranchAddress("d5_mtxout_v", &d5MtxOutV);

  TFile output(outputName, "RECREATE");
  TH1D hPD5("h_p_d5", ";p_{D5} [GeV/c];Counts", 200, 0.7, 0.8);
  TH1D hPK18("h_p_k18", ";p_{K18} [GeV/c];Counts", 200, 0.7, 0.8);
  TH1D hDp("h_dp", ";p_{D5}-p_{K18} [GeV/c];Counts", 200, -5.e-6, 5.e-6);
  TH2D hDp2D("h_dp2D", ";p_{D5} [GeV/c];p_{K18} [GeV/c]", 200, 0.7, 0.8,200,0.7,0.8);
  TH1D hChi2NdfD5("h_chi2ndf_d5", ";(#chi^{2}/ndf)_{D5};Counts", 200, 0., 10.);
  TH1D hChi2NdfK18("h_chi2ndf_k18", ";(#chi^{2}/ndf)_{K18};Counts", 200, 0., 10.);
  TH1D hDChi2Ndf("h_dchi2ndf", ";(#chi^{2}/ndf)_{D5}-(#chi^{2}/ndf)_{K18};Counts", 200, -0.0001, 0.0001);
  TH1D hDChi2("h_dchi2", ";#chi^{2}_{D5}-#chi^{2}_{K18};Counts", 200, -5., 5.);
  TH1D hDeltaD5("h_delta_d5", ";#delta_{D5} [%];Counts", 200, -10., 10.);
  TH1D hDeltaK18("h_delta_k18", ";#delta_{K18} [%];Counts", 200, -10., 10.);
  TH1D hDDelta("h_ddelta", ";#delta_{D5}-#delta_{K18} [%];Counts", 200, -0.001, 0.001);
  TH1D hDBlc1Chi2("h_dblc1_chi2", ";#chi^{2}_{BLC1,D5}-#chi^{2}_{BLC1,K18};Counts", 200, -0.002, 0.002);
  TH1D hDBlc2Chi2("h_dblc2_chi2", ";#chi^{2}_{BLC2,D5}-#chi^{2}_{BLC2,K18};Counts", 200, -0.002, 0.002);
  TH1D hD5OutMinusBcOutX("h_d5out_minus_bcout_x", ";x_{D5 global}-x_{BCout} [mm];Counts", 200, -1., 1.);
  TH1D hD5OutMinusBcOutY("h_d5out_minus_bcout_y", ";y_{D5 global}-y_{BCout} [mm];Counts", 200, -1., 1.);
  TH1D hD5OutMinusBcOutU("h_d5out_minus_bcout_u", ";u_{D5 global}-u_{BCout};Counts", 200, -0.02, 0.02);
  TH1D hD5OutMinusBcOutV("h_d5out_minus_bcout_v", ";v_{D5 global}-v_{BCout};Counts", 200, -0.02, 0.02);
  TH1D hK18OutMinusBcOutX("h_k18out_minus_bcout_x", ";x_{K18 global}-x_{BCout} [mm];Counts", 200, -1., 1.);
  TH1D hK18OutMinusBcOutY("h_k18out_minus_bcout_y", ";y_{K18 global}-y_{BCout} [mm];Counts", 200, -1., 1.);
  TH1D hK18OutMinusBcOutU("h_k18out_minus_bcout_u", ";u_{K18 global}-u_{BCout};Counts", 200, -0.02, 0.02);
  TH1D hK18OutMinusBcOutV("h_k18out_minus_bcout_v", ";v_{K18 global}-v_{BCout};Counts", 200, -0.02, 0.02);
  TH1D hD5OutMinusK18OutX("h_d5out_minus_k18out_x", ";x_{D5 global}-x_{K18 global} [mm];Counts", 200, -0.005, 0.005);
  TH1D hD5OutMinusK18OutY("h_d5out_minus_k18out_y", ";y_{D5 global}-y_{K18 global} [mm];Counts", 200, -0.005, 0.005);
  TH1D hD5OutMinusK18OutU("h_d5out_minus_k18out_u", ";u_{D5 global}-u_{K18 global};Counts", 200, -0.0001, 0.0001);
  TH1D hD5OutMinusK18OutV("h_d5out_minus_k18out_v", ";v_{D5 global}-v_{K18 global};Counts", 200, -0.0001, 0.0001);
  for (auto* hist : {&hChi2NdfD5, &hChi2NdfK18, &hDChi2Ndf, &hDDelta, &hD5OutMinusK18OutX,&hD5OutMinusK18OutY,&hD5OutMinusK18OutU,&hD5OutMinusK18OutV}) {
    hist->GetXaxis()->SetNdivisions(505);
    //hist->GetXaxis()->SetLabelSize(0.035);
  }

  UInt_t run = 0, event = 0;
  Double_t momentumD5 = 0., momentumK18 = 0., dp = 0.;
  Double_t deltaD5Percent = 0., deltaK18Percent = 0., ddeltaPercent = 0.;
  Double_t chi2D5Out = 0., chi2K18Out = 0., dchi2 = 0.;
  Double_t chi2NdfD5Out = 0., chi2NdfK18Out = 0., dchi2ndf = 0.;
  Int_t ndfD5Out = 0, ndfK18Out = 0;
  Double_t dBlc1Chi2 = 0., dBlc2Chi2 = 0.;
  Double_t d5OutMinusBcOutX = 0., d5OutMinusBcOutY = 0.;
  Double_t d5OutMinusBcOutU = 0., d5OutMinusBcOutV = 0.;
  Double_t k18OutMinusBcOutX = 0., k18OutMinusBcOutY = 0.;
  Double_t k18OutMinusBcOutU = 0., k18OutMinusBcOutV = 0.;
  Double_t d5OutMinusK18OutX = 0., d5OutMinusK18OutY = 0.;
  Double_t d5OutMinusK18OutU = 0., d5OutMinusK18OutV = 0.;
  TTree comparison("comparison", "D5--K18 comparison for one-local-track events");
  comparison.Branch("run_number", &run);
  comparison.Branch("event_number", &event);
  comparison.Branch("p_d5", &momentumD5);
  comparison.Branch("p_k18", &momentumK18);
  comparison.Branch("dp", &dp);
  comparison.Branch("delta_d5_percent", &deltaD5Percent);
  comparison.Branch("delta_k18_percent", &deltaK18Percent);
  comparison.Branch("ddelta_percent", &ddeltaPercent);
  comparison.Branch("chi2_d5", &chi2D5Out);
  comparison.Branch("chi2_k18", &chi2K18Out);
  comparison.Branch("dchi2", &dchi2);
  comparison.Branch("chi2ndf_d5", &chi2NdfD5Out);
  comparison.Branch("chi2ndf_k18", &chi2NdfK18Out);
  comparison.Branch("dchi2ndf", &dchi2ndf);
  comparison.Branch("ndf_d5", &ndfD5Out);
  comparison.Branch("ndf_k18", &ndfK18Out);
  comparison.Branch("dblc1_chi2", &dBlc1Chi2);
  comparison.Branch("dblc2_chi2", &dBlc2Chi2);
  comparison.Branch("d5out_minus_bcout_x", &d5OutMinusBcOutX);
  comparison.Branch("d5out_minus_bcout_y", &d5OutMinusBcOutY);
  comparison.Branch("d5out_minus_bcout_u", &d5OutMinusBcOutU);
  comparison.Branch("d5out_minus_bcout_v", &d5OutMinusBcOutV);
  comparison.Branch("k18out_minus_bcout_x", &k18OutMinusBcOutX);
  comparison.Branch("k18out_minus_bcout_y", &k18OutMinusBcOutY);
  comparison.Branch("k18out_minus_bcout_u", &k18OutMinusBcOutU);
  comparison.Branch("k18out_minus_bcout_v", &k18OutMinusBcOutV);
  comparison.Branch("d5out_minus_k18out_x", &d5OutMinusK18OutX);
  comparison.Branch("d5out_minus_k18out_y", &d5OutMinusK18OutY);
  comparison.Branch("d5out_minus_k18out_u", &d5OutMinusK18OutU);
  comparison.Branch("d5out_minus_k18out_v", &d5OutMinusK18OutV);

  Long64_t matched = 0;
  for (Long64_t entry = 0; entry < d5->GetEntries(); ++entry) {
    d5->GetEntry(entry);
    if (ntBLC1 != 1 || ntBLC2 != 1 || ntD5 != 1 ||
        !HasOne(pD5) || !HasOne(deltaD5) || !HasOne(chi2D5) || !HasOne(chi2NdfD5) ||
        !ndfD5 || ndfD5->size() != 1 || !pairRank || pairRank->size() != 1 ||
        !HasOne(blc1Chi2) || !HasOne(blc2Chi2) || !HasOne(blc2U) || !HasOne(blc2V) ||
        !HasOne(d5BcOutX) || !HasOne(d5BcOutY) || !HasOne(d5MtxOutX) ||
        !HasOne(d5MtxOutU) || !HasOne(d5MtxOutY) || !HasOne(d5MtxOutV) ||
        pairRank->at(0) != 0)
      continue;

    const auto found = k18ByEvent.find(EventKey(dRun, dEvent));
    if (found == k18ByEvent.end()) continue;

    run = dRun;
    event = dEvent;
    momentumD5 = pD5->at(0);
    momentumK18 = found->second.p;
    dp = momentumD5 - momentumK18;
    hDp2D.Fill(momentumD5,momentumK18);
    deltaD5Percent = deltaD5->at(0);
    deltaK18Percent = found->second.deltaPercent;
    ddeltaPercent = deltaD5Percent - deltaK18Percent;
    chi2D5Out = chi2D5->at(0);
    chi2K18Out = found->second.chi2ndf * found->second.ndf;
    dchi2 = chi2D5Out - chi2K18Out;
    chi2NdfD5Out = chi2NdfD5->at(0);
    chi2NdfK18Out = found->second.chi2ndf;
    dchi2ndf = chi2NdfD5Out - chi2NdfK18Out;
    ndfD5Out = ndfD5->at(0);
    ndfK18Out = found->second.ndf;
    dBlc1Chi2 = blc1Chi2->at(0) - found->second.blc1Chi2;
    dBlc2Chi2 = blc2Chi2->at(0) - found->second.blc2Chi2;
    d5OutMinusBcOutX = d5MtxOutX->at(0) - d5BcOutX->at(0);
    d5OutMinusBcOutY = d5MtxOutY->at(0) - d5BcOutY->at(0);
    d5OutMinusBcOutU = d5MtxOutU->at(0) - blc2U->at(0);
    d5OutMinusBcOutV = d5MtxOutV->at(0) - blc2V->at(0);
    k18OutMinusBcOutX = found->second.xOut - found->second.xBcOut;
    k18OutMinusBcOutY = found->second.yOut - found->second.yBcOut;
    k18OutMinusBcOutU = found->second.uOut - found->second.uBcOut;
    k18OutMinusBcOutV = found->second.vOut - found->second.vBcOut;
    d5OutMinusK18OutX = d5MtxOutX->at(0) - found->second.xOut;
    d5OutMinusK18OutY = d5MtxOutY->at(0) - found->second.yOut;
    d5OutMinusK18OutU = d5MtxOutU->at(0) - found->second.uOut;
    d5OutMinusK18OutV = d5MtxOutV->at(0) - found->second.vOut;
    hPD5.Fill(momentumD5);
    hPK18.Fill(momentumK18);
    hDp.Fill(dp);
    hDeltaD5.Fill(deltaD5Percent);
    hDeltaK18.Fill(deltaK18Percent);
    hDDelta.Fill(ddeltaPercent);
    hChi2NdfD5.Fill(chi2NdfD5Out);
    hChi2NdfK18.Fill(chi2NdfK18Out);
    hDChi2Ndf.Fill(dchi2ndf);
    hDChi2.Fill(dchi2);
    hDBlc1Chi2.Fill(dBlc1Chi2);
    hDBlc2Chi2.Fill(dBlc2Chi2);
    hD5OutMinusBcOutX.Fill(d5OutMinusBcOutX);
    hD5OutMinusBcOutY.Fill(d5OutMinusBcOutY);
    hD5OutMinusBcOutU.Fill(d5OutMinusBcOutU);
    hD5OutMinusBcOutV.Fill(d5OutMinusBcOutV);
    hK18OutMinusBcOutX.Fill(k18OutMinusBcOutX);
    hK18OutMinusBcOutY.Fill(k18OutMinusBcOutY);
    hK18OutMinusBcOutU.Fill(k18OutMinusBcOutU);
    hK18OutMinusBcOutV.Fill(k18OutMinusBcOutV);
    hD5OutMinusK18OutX.Fill(d5OutMinusK18OutX);
    hD5OutMinusK18OutY.Fill(d5OutMinusK18OutY);
    hD5OutMinusK18OutU.Fill(d5OutMinusK18OutU);
    hD5OutMinusK18OutV.Fill(d5OutMinusK18OutV);
    comparison.Fill();
    ++matched;
  }

  std::cout << "Matched one-local-track events: " << matched << std::endl;
  hPD5.Write();
  hPK18.Write();
  hDp.Write();
  hDp2D.Write();
  hDeltaD5.Write();
  hDeltaK18.Write();
  hDDelta.Write();
  hChi2NdfD5.Write();
  hChi2NdfK18.Write();
  hDChi2Ndf.Write();
  hDChi2.Write();
  hDBlc1Chi2.Write();
  hDBlc2Chi2.Write();
  hD5OutMinusBcOutX.Write();
  hD5OutMinusBcOutY.Write();
  hD5OutMinusBcOutU.Write();
  hD5OutMinusBcOutV.Write();
  hK18OutMinusBcOutX.Write();
  hK18OutMinusBcOutY.Write();
  hK18OutMinusBcOutU.Write();
  hK18OutMinusBcOutV.Write();
  hD5OutMinusK18OutX.Write();
  hD5OutMinusK18OutY.Write();
  hD5OutMinusK18OutU.Write();
  hD5OutMinusK18OutV.Write();
  comparison.Write();

  TString pdfName(outputName);
  if (pdfName.EndsWith(".root")) pdfName.ReplaceAll(".root", ".pdf");
  else pdfName += ".pdf";
  TCanvas canvas("comparison_canvas", "D5--K18 comparison", 900, 700);
  auto printPage = [&](TH1& hist) {
    canvas.Clear();
    hist.Draw();
    canvas.Print(pdfName);
  };
  canvas.Print(pdfName + "[");
  printPage(hPD5);
  printPage(hPK18);
  printPage(hDp);
  printPage(hDp2D);
  printPage(hDeltaD5);
  printPage(hDeltaK18);
  printPage(hDDelta);
  printPage(hChi2NdfD5);
  printPage(hChi2NdfK18);
  printPage(hDChi2Ndf);
  printPage(hDChi2);
  printPage(hDBlc1Chi2);
  printPage(hDBlc2Chi2);
  printPage(hD5OutMinusBcOutX);
  printPage(hD5OutMinusBcOutY);
  printPage(hD5OutMinusBcOutU);
  printPage(hD5OutMinusBcOutV);
  printPage(hK18OutMinusBcOutX);
  printPage(hK18OutMinusBcOutY);
  printPage(hK18OutMinusBcOutU);
  printPage(hK18OutMinusBcOutV);
  printPage(hD5OutMinusK18OutX);
  printPage(hD5OutMinusK18OutY);
  printPage(hD5OutMinusK18OutU);
  printPage(hD5OutMinusK18OutV);
  canvas.Print(pdfName + "]");
  canvas.Write();
  output.Close();
  std::cout << "Wrote " << outputName << " and " << pdfName << std::endl;
}
