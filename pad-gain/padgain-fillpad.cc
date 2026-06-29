#include "../TPCPadHelper_260416.hh"

TH2Poly* TPC_gain;

bool checkGain(int padid, double &gain){
  gain = -1.;
  
  if(padid<0 || padid >= NumOfPadTPC)return false;
  
  double val = TPC_gain->GetBinContent(padid+1);
  if(val > 100 && val < 350){
    gain = val;
    return true;
  }
  else{return false;}
}


bool IsFrameLowEdge(int layer, int row)
{
  if(layer < 0 || layer >= NumOfLayersTPC) return false;
  for(int i = 0; i < 5; i++){
    if(tpc::frameLowEdge[layer][i] < 0) continue;
    if(row == tpc::frameLowEdge[layer][i])
      return true;
  }
  return false;
}

bool IsFrameHighEdge(int layer, int row)
{
  if(layer < 0 || layer >= NumOfLayersTPC) return false;
  if(row==0)return false;
  for(int i = 0; i < 5; i++){
    if(tpc::frameHighEdge[layer][i] < 0) continue;
    if(row == tpc::frameHighEdge[layer][i])
      return true;
  }
  return false;
}

void padgain_fillpad(){
  TFile *file = new TFile("result/260626/htofcalib-padgain-combine-260626.root");
  TPC_gain = (TH2Poly*)file->Get("TPC_gain");

  TH2Poly *TPC_gain_all = new TH2Poly("TPC_gain_all","TPC_gain_all;Z [mm];X [mm]",MinZ,MaxZ,MinX,MaxX);

  double l = (586./2.)/(1+sqrt(2.));
  Double_t px[9] = {-l*(1+sqrt(2.)),-l,l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.))};
  Double_t py[9] = {l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.)),-l,l};
  
  auto pLine = new TPolyLine(9,px,py);

  Double_t X[5];
  Double_t Y[5];
  for (Int_t l=0; l<NumOfLayersTPC; ++l) {
    Double_t pLength = tpc::padParameter[l][5];
    Double_t st      = (180.-(360./tpc::padParameter[l][3]) *
                        tpc::padParameter[l][1]/2.);
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./tpc::padParameter[l][3])/180.*TMath::Pi();
    Double_t cRad    = tpc::padParameter[l][2];
    Int_t    nPad    = tpc::padParameter[l][1];
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
      TPC_gain_all->AddBin(5, X, Y);
    }
  }


  //Make reference value
  //Inner layer for target region
  double gain_sum = 0.;
  int gain_count = 0;
  for(int ipad = 0;ipad<NumOfPadTPC;ipad++){
    if(tpc::getLayerID(ipad) == 4){
      if(double gain_layer4;checkGain(ipad,gain_layer4)){
	gain_sum += gain_layer4;
	gain_count++;
      }
    }
  }
  const double gain_target = gain_sum / (double)gain_count;

  vector<bool> pad_done(NumOfPadTPC,false);
  
  for(int ilayer = 0;ilayer<NumOfLayersTPC;ilayer++){
    for(int irow = 0;irow<tpc::padParameter[ilayer][1];irow++){
      int ipad = tpc::GetPadId(ilayer,irow);

      if(tpc::IsDead(ipad)){
	TPC_gain_all->SetBinContent(ipad+1,0);
	pad_done[ipad] = true;
      }

      int section = tpc::GetSection(ipad);
      
      double fGain = TPC_gain->GetBinContent(ipad+1);
      if(double gain_check;checkGain(ipad,gain_check)){
	TPC_gain_all->SetBinContent(ipad+1,fGain);
	pad_done[ipad] = true;
      }

      //gain not defined

      else if(double gain_check;!checkGain(ipad,gain_check)){
	//Target Region
	
	if(ilayer<=4){
	  fGain = gain_target;
	  TPC_gain_all->SetBinContent(ipad+1,fGain);
	  pad_done[ipad] = true;
	}
	
	else if(ilayer>4){
	  //statistics not enough
	  double gain_pass_before = 0.;
	  double gain_pass_after = 0.;
	  if(pad_done[ipad] == false){
	    int row_false = 1;
	    double gain_false;
	    while(!checkGain(tpc::GetPadId(ilayer,irow+row_false),gain_false)){
	      row_false++; //go to next row
	    }

	    //pad before
	    {
	      if(irow == 0){
		gain_pass_before = 0.;
	      }
	      else if(IsFrameLowEdge(ilayer,irow)){
		gain_pass_before = 0.;
	      }
	      else{
		gain_pass_before = TPC_gain->GetBinContent(tpc::GetPadId(ilayer,irow-1)+1);
	      }
	    }


	    //pad after
	    {
	      if(IsFrameHighEdge(ilayer,irow+row_false-1)){
		gain_pass_after = 0.;
	      }
	      else{
		gain_pass_after = TPC_gain->GetBinContent(tpc::GetPadId(ilayer,irow+row_false)+1);
	      }
	    }


	    //get average gain
	    {
	      if(gain_pass_before == 0. && gain_pass_after == 0.){
		fGain = 0.;
	      }
	      else if(gain_pass_before == 0.){
		fGain = gain_pass_after;
	      }
	      else if(gain_pass_after == 0.){
		fGain = gain_pass_before;
	      }
	      else{
		fGain = (gain_pass_before + gain_pass_after) / 2.;
	      }
	    }
	    for(int i=0;i<row_false;i++){
	      int fpad = tpc::GetPadId(ilayer,irow+i);
	      pad_done[fpad] = true;
	      TPC_gain_all->SetBinContent(fpad+1,fGain);
	    }
	  }
	}
      }
    }//row
  }//layer
  


  TFile *f = new TFile("result/260626/padgain-fillpad-260626.root","RECREATE");
  TPC_gain->Write();
  TPC_gain_all->Write();
  f->Close();

}
