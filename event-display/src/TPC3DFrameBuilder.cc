#include "TPC3DFrameBuilder.hh"

#include <TMath.h>
#include <TGeoTube.h>
#include <TGeoBBox.h>
#include <TGeoTrd1.h>
#include <TGeoMatrix.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoXtru.h>

TPC3DFrameBuilder::TPC3DFrameBuilder(){
}

TGeoVolume* TPC3DFrameBuilder::TPC3DGeometry(){
  Double_t flength = 586;
  Double_t fheight = 550;
  Double_t edge = flength/(1+sqrt(2));
  Double_t tan = TMath::Tan(22.5*TMath::Pi()/180.);

  TGeoManager *geom = new TGeoManager("HS","HS");
  TGeoMaterial *mat_world = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium *med_world = new TGeoMedium("med_world",1,mat_world);
  TGeoVolume *world = geom->MakeBox("world",med_world,flength*2,flength*2,flength*2);
  geom->SetTopVolume(world);
  geom->SetVisLevel(0);
  world->SetVisibility(kFALSE);
  
  
  TGeoMaterial *mat_tpc = new TGeoMaterial("Al",26.98,13,2.7);
  TGeoMedium *med_tpc = new TGeoMedium("med_tpc",1,mat_tpc);
  TGeoXtru *xtru_tpc_inside = new TGeoXtru(2);
  Double_t zin[] = { -tan*flength/2, -flength/2, -flength/2, -tan*flength/2,
		     tan*flength/2, flength/2, flength/2, tan*flength/2};
  Double_t xin[] = { -flength/2, -edge/2, edge/2, flength/2,
		     flength/2, edge/2, -edge/2, -flength/2};  
  Double_t yin[] = { -fheight/2, fheight/2};
  Double_t scale[] = {1., 1.};
  Double_t x0[] = {0., 0.};
  Double_t z0[] = {0., 0.};

  Int_t nxz = sizeof(xin)/sizeof(Double_t);
  xtru_tpc_inside->DefinePolygon(nxz,xin,zin);
  
  Int_t n;
  Int_t ny = sizeof(yin)/sizeof(Double_t);
  for (n = 0 ; n < ny ; n++)
    xtru_tpc_inside->DefineSection(n,yin[n],x0[n],z0[n],scale[n]);  

  TGeoVolume *tpc_inside = new TGeoVolume("tpc_inside",xtru_tpc_inside,med_tpc);

  ////////////////////////////////////////////////////////////////////////////
  flength = 586 - 10;
  edge = flength/(1+sqrt(2));

  TGeoXtru *xtru_tpc_outside = new TGeoXtru(2);
  Double_t z[] = { -tan*flength/2, -flength/2, -flength/2, -tan*flength/2,
		   tan*flength/2, flength/2, flength/2, tan*flength/2};
  Double_t x[] = { -flength/2, -edge/2, edge/2, flength/2,
		   flength/2, edge/2, -edge/2, -flength/2};  
  Double_t y[] = { -fheight/2, fheight/2};
  
  xtru_tpc_outside->DefinePolygon(nxz,x,z);  
  for (n = 0 ; n < ny ; n++)
    xtru_tpc_outside->DefineSection(n,y[n],x0[n],z0[n],scale[n]);  
  
  TGeoVolume *tpc_outside = new TGeoVolume("tpc_outside",xtru_tpc_outside,med_tpc);

  //E72 LH2 Target

  TGeoMaterial *mat_tgt = new TGeoMaterial("H",24,12,1.3);
  TGeoMedium *med_tgt = new TGeoMedium("med_tgt",1,mat_tgt);

  TGeoTube *tube_tgt = new TGeoTube("tube_tgt",0,80/2,100/2);
  TGeoVolume *tgt = new TGeoVolume("tgt",tube_tgt,med_tgt);


  //E42 CH2 Target
  /*
  TGeoMaterial *mat_tgt = new TGeoMaterial("C",24,12,1.3);
  TGeoMedium *med_tgt = new TGeoMedium("med_tgt",1,mat_tgt);

  TGeoBBox *box_tgt = new TGeoBBox("box_tgt",0.5*30.2,0.5*20.2,0.5*20.2);
  TGeoVolume *tgt = new TGeoVolume("tgt",box_tgt,med_tgt);
  */

  TGeoVolumeAssembly *tpc = new TGeoVolumeAssembly("tpc");  
  //tpc->AddNode(tgt,1,new TGeoTranslation(-143.,0.,0.));
  tpc->AddNode(tgt,1,new TGeoTranslation(0.,-143.,0.));
  tpc->AddNode(tpc_outside,1);
  tpc->AddNode(tpc_inside,1);
  
  world->AddNode(tpc,1);

  geom->CloseGeometry();
  geom->SetVisLevel(4);

  return world;
}
