#include "TPC2DFrameBuilder.hh"
#include "TPCPadHelper.hh"

#include <TH2Poly.h>
#include <TMath.h>


TPC2DFrameBuilder::TPC2DFrameBuilder(){
}

TH2Poly* TPC2DFrameBuilder::TPC2DGeometry(){

  TH2Poly *TPC2DPoly = new TH2Poly("TPC2DPoly", "TPC2DPoly", MinZ, MaxZ, MinX, MaxX);
  
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
      Bool_t dead = false;
      for(int i=0;i<sizeof(padOnFrame)/sizeof(*padOnFrame);i++){
	if(GetPadId(l,j) == padOnFrame[i]){
	  dead=true;
	}
	
      }

      if(!dead){
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
      }
      TPC2DPoly->AddBin(5, X, Y);

    }
  }

  return TPC2DPoly;
  
}
