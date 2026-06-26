#include "../TPCPadHelper.hh"
#include "../TPCEventDisplayHelper.hh"

void tpc_2d(){
  auto TPC_evt2d = tpcdisp::MakeTPCPadMap("TPC_evt2d", "TPC evt2d");
  auto c1 = new TCanvas("c1","c1");
  TPC_evt2d->GetZaxis()->SetTitle("Layer Id");
  TPC_evt2d->Draw("colz");
  
  tpcdisp::DrawTarget();
  tpcdisp::DrawTargetHolder();

  for(int i=0;i<NumOfPadTPC;i++){
    TPC_evt2d->SetBinContent(i+1,tpc::getLayerID(i));
  }

  c1->SaveAs("tpc-2d.pdf");

  
}
