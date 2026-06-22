double DensityEffectCorrection(Double_t betagamma, Double_t *par){

    //reference : Sternheimer’s parameterizatio(PDG)
    //notation : par[0] : a, par[1] : k, par[2] : x0, par[3] : x1, par[4] : _C, par[5] : delta0
    Double_t constant = 2*TMath::Log(10);
    Double_t delta = 0.;
    Double_t X = log10(betagamma);
    if(X<=par[2]) delta = par[5]*TMath::Power(10., 2*(X - par[2]));
    else if(par[2]<X && X<par[3]) delta = constant*X - par[4] + par[0]*pow((par[3] - X), par[1]);
    else if(X>=par[3]) delta = constant*X - par[4];

  return delta;

}


double Caldedx(double mass, double beta){
  Double_t rho=TMath::Power(10.,-3)*(0.9*1.662 + 0.1*0.6672); //[g cm-3] density of the material
  Double_t ZoverA = 17.2/37.6; //[mol g-1]
  
  Double_t I = 0.9*188.0 + 0.1*41.7; //[eV] Mean excitation energy : 10 eV * Z
  Double_t density_effect_par[6]={0.}; //Sternheimer’s parameterization
  density_effect_par[0] = 0.9*0.19714 + 0.1*0.09253;
  density_effect_par[1] = 0.9*2.9618 + 0.1*3.6257;
  density_effect_par[2] = 0.9*1.7635 + 0.1*1.6263;
  density_effect_par[3] = 0.9*4.4855 + 0.1*3.9716;
  density_effect_par[4] = 0.9*11.9480 + 0.1*9.5243;
  density_effect_par[5] = 0.;

  Double_t Z = 1;
  Double_t me = 0.5109989461; //[MeV]
  Double_t K = 0.307075; //[MeV mol-1 cm2]
  Double_t constant = rho*K*ZoverA; //[MeV cm-1]

  Double_t I2 = I*I;
  Double_t beta2 = beta*beta;
  Double_t gamma2 = 1./(1.-beta2);
  Double_t MeVToeV = TMath::Power(10.,6);
  Double_t Wmax = 2*me*beta2*gamma2/((me/mass+1.)*(me/mass+1.)+2*(me/mass)*(TMath::Sqrt(gamma2)-1)); //[MeV]
  Double_t delta = DensityEffectCorrection(TMath::Sqrt(beta2*gamma2), density_effect_par);
  Double_t dedx = constant*Z*Z/beta2*(0.5*TMath::Log(2*me*beta2*gamma2*Wmax*MeVToeV*MeVToeV/I2)-beta2-0.5*delta);

  Double_t MeVTokeV = TMath::Power(10.,3);
  return dedx * MeVTokeV; // [keV cm-1]
}


void Bethe_P10(){
  auto can = new TCanvas("c1","c1",2000,2000);
  auto gr_p = new TGraph();
  auto gr_k = new TGraph();
  auto gr_pi = new TGraph();
  auto gr_e = new TGraph();
  double me = 0.510; //MeV/c2
  double mp = 938.282; //MeV/c2
  double mk = 493.677; //MeV/c2
  double mpi = 139.570; //MeV/c2
  double mom_start = 50; //MeV/c
  double mom_end = 1100; //MeV/c
  int n = 100;
  double step = (mom_end - mom_start) / n;
  for(int i=0;i<n;i++){
    double mom = mom_start + i*step;
    double beta_p = mom/TMath::Sqrt(mom*mom+mp*mp);
    double beta_k = mom/TMath::Sqrt(mom*mom+mk*mk);
    double beta_pi= mom/TMath::Sqrt(mom*mom+mpi*mpi);
    double beta_e = mom/TMath::Sqrt(mom*mom+me*me);

    gr_p->AddPoint(mom,Caldedx(mp,beta_p));
    gr_k->AddPoint(mom,Caldedx(mk,beta_k));
    gr_pi->AddPoint(mom,Caldedx(mpi,beta_pi));
    gr_e->AddPoint(mom,Caldedx(me,beta_e));
  }

  auto mg = new TMultiGraph();
  gr_p->SetLineColor(kBlue);
  gr_p->SetLineWidth(4);
  gr_k->SetLineColor(kBlack);
  gr_k->SetLineWidth(4);
  gr_pi->SetLineColor(kRed);
  gr_pi->SetLineWidth(4);
  gr_e->SetLineColor(kMagenta);
  gr_e->SetLineWidth(4);
  
  
  mg->Add(gr_p);
  mg->Add(gr_k);
  mg->Add(gr_pi);
  mg->Add(gr_e);
  

  
  TLegend *le = new TLegend(0.8,0.5,0.48,0.6);
  le->AddEntry(gr_e,"e^{-}");
  le->AddEntry(gr_p,"p");
  le->AddEntry(gr_k,"K^{#pm}");
  le->AddEntry(gr_pi,"#pi^{#pm}");
  
  
  can->cd();
  //gPad->SetLogx();
  gPad->SetLogy();
  mg->SetTitle("Bethe-bloch formula;p [MeV/c];Energy loss inside the P10 gas [keV/cm]");
  mg->GetXaxis()->SetRangeUser(60,1000);
  mg->Draw("a");
  le->Draw("same");
}
