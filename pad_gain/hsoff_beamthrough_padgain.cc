const int runnumber[7] = {2570, 2581, 2585, 2587, 2589, 2590, 2592};

const double z_tgt = -143.; //mm
const double z_tpccenter_to_ff = -150.; //mm

const Int_t NumOfPadTPC       = 5768;
const Int_t NumOfLayersTPC    = 32;

const Double_t MinX = -300;
const Double_t MaxX = 300;
const Double_t MinZ = -300;
const Double_t MaxZ = 300;

const Double_t ZTarget = -143.;
//_____________________________________________________________________________
//#OfPad #division #radius padLength
static const Double_t padParameter[NumOfLayersTPC][6] =
{{0, 48,    14.75, 48, 0,  9.},
 {1, 48,    24.25, 48, 0,  9.},
 {2, 72,    33.75, 72, 0,  9.},
 {3, 96,    43.25, 96, 0,  9.},
 {4, 120,    52.75,120,0,   9.},
 {5, 144,    62.25,144,0,   9.},
 {6, 168,    71.75,168,0,   9.},
 {7, 192,    81.25,192,0,   9.},
 {8, 216,    90.75,216,0,   9.},
 {9, 240,    100.25,240,0,  9.},
 {10,208,    111.5,241, 0,  12.5},
 {11,218,    124.5,271, 0,  12.5},
 {12,230,    137.5,300, 0,  12.5},
 {13,214,    150.5,330, 0,  12.5},
 {14,212,    163.5,360, 0,  12.5},
 {15,214,    176.5,390, 0,  12.5},
 {16,220,    189.5,420, 0,  12.5},
 {17,224,    202.5,449, 0,  12.5},
 {18,232,    215.5,479, 0,  12.5},
 {19,238,    228.5,509, 0,  12.5},
 {20,244,    241.5,539, 0,  12.5},
 {21,232,    254.5,569, 0,  12.5},
 {22,218,    267.5,599, 0,  12.5},
 {23,210,    280.5,628, 0,  12.5},
 {24,206,    293.5,658, 0,  12.5},
 {25,202,    306.5,688, 0,  12.5},
 {26,200,    319.5,718, 0,  12.5},
 {27,196,    332.5,748, 0,  12.5},
 {28,178,    345.5,777, 0,  12.5},
 {29,130,    358.5,807, 0,  12.5},
 {30,108,    371.5,837, 0,  12.5},
 {31,90,     384.5,867, 0, 12.5}};

//_____________________________________________________________________________
Int_t GetPadId(Int_t layer_id, Int_t row_id)
{
  // Return -1 (invalid) for negative row_id to prevent crashes.
  // This is primarily to handle the INT_MIN(-2147483648) value resulting from (Int_t)NaN
  // when initializing TPCHit with NaN in TPCCluster as below;
  //  m_mean_hit(new TPCHit(layer, TMath::QuietNaN())) <- this
  if (row_id < 0) return -1;

  Int_t pad_id = 0;
  for (Int_t layer = 0; layer < layer_id; layer++)
    pad_id += static_cast<Int_t>(padParameter[layer][1]);
  pad_id += row_id;
  return pad_id;
}

//_____________________________________________________________________________
static const Int_t padOnSectionFrame[] =
{
//Gem section 1 frame 0
1503 ,1504, 1505 ,1506 ,1507 ,1508 ,1509 ,1510 ,1511 ,1512 ,1513 ,1514 ,1515 ,1728 ,1729 ,1730,

//Gem section 1 frame 1
72 ,134 ,135 ,222 ,332 , 408,  467 ,626 ,809 ,1016, 59 ,110 ,185 ,284,555 ,726, 921 ,1140 ,1363, 1564 ,1565 ,1778,

//Gem section 2 frame 0, frame 1 + adjacent pads
1628,1629,1630,1817,1818,1819,1820,1849,1850,1851,2035,2036,2037,2038,2067,2068,2069,2244,2245,2246,2247,2248,2276,2277,2278,2454,2455,2456,2457,2485,2486,2487,2488,2667,2668,2669,2670,2698,2699,2700,2701,2886,2887,2888,2917,2918,2919,3110,3111,3112,3113,3141,3142,3143,3341,3342,3343,3344,3373,3374,3375,3578,3579,3580,3581,3610,3611,3612,3813,3814,3815,3844,3845,3846,4034,4035,4036,4065,4066,4067,4276,4277,4278,4480,4481,4482,4680,4681,4682,4878,4879,4880,5072,5073,5074,

//Gem section 3 frame 0
4720 ,4721 ,4722 ,4723 ,4724 ,4725 ,4726 ,4727 ,4728 ,4929 ,4930 ,4931 ,4932 ,4933 ,4934 ,4935 ,4936 ,5134 ,5135 ,5136 ,5137 ,5138 ,5139 ,5140 ,5328 ,5329 ,5330 ,5331 ,5332 ,5333 ,5487 ,5488 ,5489 ,5490 ,5491 ,5492 ,5612 ,5613 ,5614 ,5615 ,5616 ,5716 ,5717 ,5718 ,5719 ,5720, 5532 ,5533, 5612 ,5613 ,5614 ,5615 ,5616 ,5654 ,5655 ,5716 ,5717 ,5718 ,5719 ,5720 ,5757 ,5758,

//Gem section 3 frame 1
3413 ,3414 ,3415 ,3416 ,3417 ,3418 ,3419 ,3660 ,3661 ,3662 ,3663 ,3664 ,3665 ,3904 ,3905 ,3906 ,3907 ,3908 ,4134 ,4135 ,4136 ,4137 ,4353 ,4354 ,4355 ,4356 ,4565 ,4566 ,4567 ,4568 ,4774 ,4775 ,4776 ,4981 ,4982 ,5183 ,5184 ,5372 ,5373 ,5374 ,5530 ,5531 ,5532 ,5533 ,5654 ,5655 ,5757 ,5758,

//Gem section 4 frame 0
4409 ,4410 ,4411 ,4412 ,4413 ,4414 ,4415 ,4416 ,4417 ,4418 ,4419 ,4420 ,4421 ,4422 ,4423 ,4424 ,4425 ,4426 ,4427 ,4428 ,4429 ,4430 ,4431 ,4432 ,4433 ,4434 ,4435 ,4436 ,4437 ,4438 ,4439 ,4440 ,4441 ,4442 ,4443 ,4444 ,4445 ,4446 ,4447 ,4448 ,4449 ,4450 ,4451 ,4452 ,4611 ,4612 ,4613 ,4614 ,4615 ,4616 ,4617 ,4618 ,4619,

//Gem section 4 frame 1
2794 ,2795 ,2796 ,2797 ,2798 ,2799 ,2800 ,2801 ,2802 ,2803 ,2804 ,2805 ,2806 ,2807 ,2808 ,2809 ,3001 ,3002 ,3003 ,3004 ,3005 ,3006 ,3007 ,3008 ,3009 ,3010 ,3011 ,3012 ,3013 ,3014 ,3015 ,3016 ,3017, 3018, 3037, 3038, 3039 ,3040 ,3041 ,3042 ,3043 ,3044 ,3045 ,3046 ,3047 ,3048 ,3049 ,3050 ,3051 ,3052 ,3053 ,3054 ,3227 ,3228 ,3229 ,3230 ,3231 ,3289 ,3290 ,3291 ,3292 ,3293 ,3294 ,3295 ,3296 ,3297 ,3540 ,3541 ,3542 ,3543 ,3544 ,3545 ,3546
};

//_____________________________________________________________________________
inline Bool_t Noise(Int_t pad_id){
  Bool_t noise = std::find(std::begin(padOnSectionFrame), std::end(padOnSectionFrame), pad_id) != std::end(padOnSectionFrame);
  if(noise) return true;
  else return false;
}

double DCA(double xbcout, double ybcout, double ubcout, double vbcout,
           double xtpc, double ytpc, double utpc, double vtpc)
{
  const TVector3 p1(xbcout, ybcout, 0.0);
  //TPC x0, y0 is at z = TGT -> Change to z = FF;
  const TVector3 p2(xtpc + utpc*(-z_tgt -z_tpccenter_to_ff) , ytpc + vtpc*(-z_tgt - z_tpccenter_to_ff), 0.0);
  
  const TVector3 v1vec(ubcout, vbcout, 1.0);
  const TVector3 v2vec(utpc, vtpc, 1.0);

  const TVector3 w0 = p1 - p2;

  const double a = v1vec.Dot(v1vec);
  const double b = v1vec.Dot(v2vec);
  const double c = v2vec.Dot(v2vec);
  const double d = v1vec.Dot(w0);
  const double e = v2vec.Dot(w0);

  const double denom = a*c - b*b;
  const double eps = 1e-12;

  double s = 0.0;
  double t = 0.0;

  if (TMath::Abs(denom) > eps) {
    // general (not parallel)
    s = (b*e - c*d) / denom;
    t = (a*e - b*d) / denom;
  } else {
    // nearly parallel: take t=0 and minimize distance to line2 point p2
    // s = argmin |(p1 + s*v1) - p2|^2
    s = -d / a;
    t = 0.0;
  }

  const TVector3 cp1 = p1 + s * v1vec;
  const TVector3 cp2 = p2 + t * v2vec;

  /*
  if (c1) *c1 = cp1;
  if (c2) *c2 = cp2;
  if (s_out) *s_out = s;
  if (t_out) *t_out = t;
  */

  return (cp1 - cp2).Mag();
}

double cos_theta(double ubcout, double vbcout, double utpc, double vtpc)
{
    const TVector3 d1(ubcout, vbcout, 1.0);
    const TVector3 d2(utpc, vtpc, 1.0);

    const double n1 = d1.Mag();
    const double n2 = d2.Mag();
    if (n1 == 0.0 || n2 == 0.0) return -2.0; // invalid marker

    return d1.Dot(d2) / (n1 * n2);
}


void hsoff_beamthrough_padgain(){

  TH1D *hist_dca = new TH1D("hist_dca","hist_dca",1000,0,1000);
  TH1D *hist_cos = new TH1D("hist_cos","hist_cos",100,-1,1);
  TH2D *hist_dca_cos = new TH2D("hist_dca_cos","hist_dca_cos",1000,0,1000,100,-1,1);



  TH1D *hist_res_x = new TH1D("hist_res_x","hist_res_x",100,-50,50);
  TH1D *hist_res_y = new TH1D("hist_res_y","hist_res_y",100,40,100);
  TH2D *hist_res = new TH2D("hist_res","hist_res",100,-100,100,100,0,120);

  TF1 *f_res_x = new TF1("f_res_x","gaus(0)",-10,10);
  TF1 *f_res_y = new TF1("f_res_y","gaus(0)",50,100);

  TH1D *hist_all_track = new TH1D("hist_all_track","hist_all_track",6,-0.5,5.5);
  TH1D *hist_cor_track = new TH1D("hist_cor_track","hist_cor_track",6,-0.5,5.5);
  
  TH2D *hist_bcout_pass = new TH2D("hist_bcout_pass","hist_bcout_pass",300,-100,100,300,-100,100);
  TH2D *hist_bcout_fail = new TH2D("hist_bcout_fail","hist_bcout_fail",300,-100,100,300,-100,100);

  TH2D *hist_bcout_uv_pass = new TH2D("hist_bcout_uv_pass","hist_bcout_uv_pass",300,-0.5,0.5,300,-0.5,0.5);
  TH2D *hist_bcout_uv_fail = new TH2D("hist_bcout_uv_fail","hist_bcout_uv_fail",300,-0.5,0.5,300,-0.5,0.5);

  auto TPC_gain = new TH2Poly("TPC_gain", "TPC_gain;Z;X", MinZ, MaxZ, MinX, MaxX);
    auto TPC_count = new TH2Poly("TPC_count", "TPC_count;Z;X", MinZ, MaxZ, MinX, MaxX);
    
    double l = (586./2.)/(1+sqrt(2.));
    Double_t px[9] = {-l*(1+sqrt(2.)),-l,l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.))};
    Double_t py[9] = {l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.)),-l,l};
  
  auto pLine = new TPolyLine(9,px,py);

  Double_t X[5];
  Double_t Y[5];
  for (Int_t l=0; l<NumOfLayersTPC; ++l) {
    Double_t pLength = padParameter[l][5];
    Double_t st      = (180.-(360./padParameter[l][3]) *
                        padParameter[l][1]/2.);
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./padParameter[l][3])/180.*TMath::Pi();
    Double_t cRad    = padParameter[l][2];
    Int_t    nPad    = padParameter[l][1];
    for (Int_t j=0; j<nPad; ++j) {
      X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[0] = X[4];
      Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[0] = Y[4];
      for (Int_t k=0; k<5; ++k) X[k] += ZTarget;
      TPC_gain->AddBin(5, X, Y);
      TPC_count->AddBin(5, X, Y);
    }
  }

  for(int i=0;i<NumOfPadTPC;i++){
    TPC_count->SetBinContent(i+1,0);
  }
  
  for (size_t i = 0; i < 1; ++i) {
    int run = runnumber[i];
    //for(int run : runnumber){
  
    string dir = "/gpfs/group/had/sks/Users/haein/JPARC2025Nov_root/hsoff_beamthrough";
    TFile *file_hodo = new TFile(Form("%s/run0%d_Hodoscope.root",dir.c_str(),run));
    TFile *file_bcout = new TFile(Form("%s/run0%d_BcOutTracking.root",dir.c_str(),run));
    TFile *file_tpc = new TFile(Form("%s/run0%d_DstTPCTracking.root",dir.c_str(),run));

    TTree *tree_hodo = (TTree*)file_hodo->Get("hodo");
    TTree *tree_bcout = (TTree*)file_bcout->Get("bcout");
    TTree *tree_tpc = (TTree*)file_tpc->Get("tpc");

    //Hodo
    double btof0;
    tree_hodo->SetBranchAddress("btof0",&btof0);

    //BcOut
    int ntrack;
    vector<double>* chisqr = nullptr;
    vector<double>* x0 = nullptr;
    vector<double>* y0 = nullptr;
    vector<double>* u0 = nullptr;
    vector<double>* v0 = nullptr;
    tree_bcout->SetBranchAddress("ntrack",&ntrack);
    tree_bcout->SetBranchAddress("chisqr",&chisqr);
    tree_bcout->SetBranchAddress("x0",&x0);
    tree_bcout->SetBranchAddress("y0",&y0);
    tree_bcout->SetBranchAddress("u0",&u0);
    tree_bcout->SetBranchAddress("v0",&v0);

    //TPCTracking
    int nclTpc;
    int ntTpc;
    vector<double>* cluster_x = nullptr;
    vector<double>* cluster_y = nullptr;
    vector<double>* cluster_z = nullptr;
    vector<double>* cluster_de = nullptr;
    vector<int>* cluster_size = nullptr;
    vector<int>* cluster_layer = nullptr;
    vector<int>* cluster_row_center = nullptr;
    vector<double>* cluster_mrow = nullptr;

    vector<int>* nhtrack = nullptr;
    vector<double>* x0Tpc = nullptr;
    vector<double>* y0Tpc = nullptr;
    vector<double>* u0Tpc = nullptr;
    vector<double>* v0Tpc = nullptr;

    vector<vector<double>>* hitlayer = nullptr;
    vector<vector<double>>* track_cluster_de = nullptr;
    vector<vector<double>>* track_cluster_size = nullptr;
    vector<vector<double>>* track_cluster_mrow = nullptr;
    
    
    


    tree_tpc->SetBranchAddress("nclTpc",&nclTpc);
    tree_tpc->SetBranchAddress("ntTpc",&ntTpc);
    tree_tpc->SetBranchAddress("cluster_x",&cluster_x);
    tree_tpc->SetBranchAddress("cluster_y",&cluster_y);
    tree_tpc->SetBranchAddress("cluster_z",&cluster_z);
    tree_tpc->SetBranchAddress("cluster_de",&cluster_de);
    tree_tpc->SetBranchAddress("cluster_size",&cluster_size);
    tree_tpc->SetBranchAddress("cluster_layer",&cluster_layer);
    tree_tpc->SetBranchAddress("cluster_row_center",&cluster_row_center);
    tree_tpc->SetBranchAddress("cluster_mrow",&cluster_mrow);

    tree_tpc->SetBranchAddress("nhtrack",&nhtrack);
    tree_tpc->SetBranchAddress("x0Tpc",&x0Tpc);
    tree_tpc->SetBranchAddress("y0Tpc",&y0Tpc);
    tree_tpc->SetBranchAddress("u0Tpc",&u0Tpc);
    tree_tpc->SetBranchAddress("v0Tpc",&v0Tpc);
    
    tree_tpc->SetBranchAddress("hitlayer",&hitlayer);
    tree_tpc->SetBranchAddress("track_cluster_de",&track_cluster_de);
    tree_tpc->SetBranchAddress("track_cluster_size",&track_cluster_size);
    tree_tpc->SetBranchAddress("track_cluster_mrow",&track_cluster_mrow);


    bool IsPion = false;
    bool IsTPCTrack = false;
    bool IsBcoutTrack = false;
      
    for(int n=0;n<tree_hodo->GetEntries();n++){
    //for(int n=0;n<100;n++){
      tree_hodo->GetEntry(n);
      tree_bcout->GetEntry(n);
      tree_tpc->GetEntry(n);

      IsPion = false;
      IsTPCTrack = false;
      IsBcoutTrack = false;
      int IsPionTrack = 0;
      vector<int> pion_track;
      
      if(btof0 > -3)IsPion = true;
      //Tight cut to BcOut 
      if(ntrack == 1){
	if((*chisqr)[0] < 3.)
	  IsBcoutTrack = true;
      }
      if(ntTpc > 0)IsTPCTrack = true;

      if(IsPion && IsBcoutTrack && IsTPCTrack){
	hist_all_track->Fill(ntTpc);
	//Find Beam Track
	IsPionTrack = 0;
	for(int nbcouttr = 0;nbcouttr<ntrack;nbcouttr++){
	  for(int ntpctr = 0;ntpctr<ntTpc;ntpctr++){
	    double min_dis = DCA((*x0)[nbcouttr],(*y0)[nbcouttr],(*u0)[nbcouttr],(*v0)[nbcouttr],(*x0Tpc)[ntpctr],(*y0Tpc)[ntpctr],(*u0Tpc)[ntpctr],(*v0Tpc)[ntpctr]);
	    hist_dca->Fill(min_dis);
	    double track_angle = cos_theta((*u0)[nbcouttr],(*v0)[nbcouttr],(*u0Tpc)[ntpctr],(*v0Tpc)[ntpctr]);
	    hist_cos->Fill(track_angle);
	    hist_dca_cos->Fill(min_dis,track_angle);


	    double xbcout_ff = (*x0)[nbcouttr];
	    double ybcout_ff = (*y0)[nbcouttr];
	    double xtpc_ff = (*x0Tpc)[ntpctr] + (*u0Tpc)[ntpctr]*(-z_tgt -z_tpccenter_to_ff);
	    double ytpc_ff = (*y0Tpc)[ntpctr] + (*v0Tpc)[ntpctr]*(-z_tgt -z_tpccenter_to_ff);

	    hist_res_x->Fill(xbcout_ff - xtpc_ff);
	    hist_res_y->Fill(ybcout_ff - ytpc_ff);
	    hist_res->Fill(xbcout_ff - xtpc_ff,ybcout_ff - ytpc_ff);

	    double res_x = xbcout_ff - xtpc_ff;
	    double res_y = ybcout_ff - ytpc_ff;

	    if(res_x > -10 && res_x < 10 && res_y > 55 && res_y < 80){
	      IsPionTrack++;
	      pion_track.push_back(ntpctr);
	    }
	  }
	}
	hist_cor_track->Fill(IsPionTrack);
	if(IsPionTrack ==0){
	  hist_bcout_fail->Fill((*x0)[0],(*y0)[0]);
	  hist_bcout_uv_fail->Fill((*u0)[0],(*v0)[0]);
	}
	else if(IsPionTrack >0){
	  hist_bcout_pass->Fill((*x0)[0],(*y0)[0]);
	  hist_bcout_uv_pass->Fill((*u0)[0],(*v0)[0]);
	}

	if(IsPionTrack == 1){
	  int pitr = pion_track[0];
	  for(int i=0;i<(*nhtrack)[pitr];i++){
	    int padid = GetPadId((*hitlayer)[pitr][i],std::round((*track_cluster_mrow)[pitr][i]));
	    int count = TPC_count->GetBinContent(padid+1);
	    count +=1;
	    TPC_count->SetBinContent(padid+1,count);
	  }
	}
      }
    }


      
      
      
            
  }


  auto c1 = new TCanvas("c1","c1");
  c1->Divide(3);
  c1->cd(1);
  hist_dca->Draw();
  c1->cd(2);
  hist_cos->Draw();
  c1->cd(3);
  hist_dca_cos->Draw("colz");

  auto c2 = new TCanvas("c2","c2");
  c2->Divide(2,2);
  c2->cd(1);
  hist_res_x->Draw();
  hist_res_x->Fit("f_res_x");
  c2->cd(2);
  hist_res_y->Draw();
  hist_res_y->Fit("f_res_y");
  c2->cd(3);
  hist_res->Draw("colz");
  c2->cd(4);
  
  hist_cor_track->SetLineColor(kRed);
  hist_cor_track->Draw();
  hist_all_track->Draw("same");

  auto c3 = new TCanvas("c3","c3");
  c3->Divide(2,2);
  c3->cd(1);
  hist_bcout_pass->Draw("colz");
  c3->cd(2);
  hist_bcout_fail->Draw("colz");
  c3->cd(3);
  hist_bcout_uv_pass->Draw("colz");
  c3->cd(4);
  hist_bcout_uv_fail->Draw("colz");

  auto c4 = new TCanvas("c4","c4");
  c4->cd();
  TPC_count->Draw("colz");
  
    
}
