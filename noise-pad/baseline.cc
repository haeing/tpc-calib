#include "../TPCPadHelper_260416.hh"

const Int_t NumOfCoBo        = 8;
const Int_t NumOfAsAd        = 31; // #32 is unused
const Int_t NumOfAsAdPerCoBo = 4;
const Int_t NumOfPad         = 5768;
const Int_t NumOfAGETPerAsAd = 4;
const Int_t NumOfChannelAGET = 68;

Int_t     channel_map[NumOfAGETPerAsAd][NumOfAsAd][NumOfChannelAGET];

void baseline(int runnumber){
  gROOT->SetBatch(kTRUE);
  string dir = "/hsm/had/sks/E72/JPARC2025Nov/rootfile/tpcthreshold";
  TFile *file = new TFile(Form("%s/TpcThreshold_%d.root",dir.c_str(),runnumber));
  string outpdf = Form("%s/baseline_run0%d.pdf",dir.c_str(),runnumber);
  TFile *f = new TFile(Form("%s/baseline_run0%d.root",dir.c_str(),runnumber),"RECREATE");
  
  //start : Make Channel Map
  {
    std::ifstream ifs("channel_map_20180522.param");
    if( !ifs.is_open() ){
      std::cerr<<" --> file open fail ..."<<std::endl;
      std::exit( EXIT_FAILURE );
    }


    while( ifs.good() ){
      std::string buf;
      std::getline(ifs,buf);
      if(buf[0]=='#' || buf.empty())continue;
      std::istringstream is(buf);
      std::istream_iterator<std::string> issBegin(is);
      std::istream_iterator<std::string> issEnd;
      std::vector<std::string> param(issBegin,issEnd);
      if(param.empty()||param[0].empty())continue;
      Int_t aget  = std::atoi( param[0].c_str() );
      Int_t asad  = std::atoi( param[1].c_str() );
      Int_t ch2   = std::atoi( param[2].c_str() );
      Int_t layer = std::atoi( param[3].c_str() );
      Int_t row   = std::atoi( param[4].c_str() );
      Int_t pid   = std::atoi( param[5].c_str() );
      if( aget!=0 ){
	Int_t ch;
	if(ch2>=1 && ch2<=11) ch = ch2-1;
	else if(ch2>=12 && ch2<=21) ch = ch2;
	else if(ch2>=22 && ch2<=43) ch = ch2+1;
	else if(ch2>=44 && ch2<=53) ch = ch2+2;
	else if(ch2>=54 && ch2<=64) ch = ch2+3;
	channel_map[aget-1][asad-1][ch] = pid;
      }
    }
  }
  //end :  make channel map

  auto c1 = new TCanvas("c1","c1");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TPaveText *p = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
  p->AddText("baseline.cc");
  p->AddText("TPC Pad Baseline check");
  TDatime now;
  p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",now.GetYear(),now.GetMonth(),now.GetDay(),now.GetHour(),now.GetMinute(),now.GetSecond()));
  p->Draw();
  c1->Print((outpdf + "(").c_str());
  
  int count = 0;
  for( Int_t cobo=0; cobo<NumOfCoBo; ++cobo ){
    for( Int_t asad=0; asad<NumOfAsAdPerCoBo; ++asad ){
      for( Int_t aget=0; aget<NumOfAGETPerAsAd; ++aget ){
	for( Int_t ch=0; ch<NumOfChannelAGET; ++ch ){
	  int padid = channel_map[aget][cobo*4+asad][ch] - 1;
	  if(padid <0 || padid >=NumOfPad)continue;

	  if(count%30 == 0){
	    if(count !=0)c1->Print(outpdf.c_str());
	    c1->Clear();
	    c1->Divide(6,5);
	  }
	  c1->cd(count%30+1);
	  auto h  = (TH2F*)file->Get(Form("h_fadc_%d%d%d%02d", cobo, asad, aget, ch));
	  if(h){
	    h->SetName(Form("h_fadc%d",padid));
	    h->SetTitle(Form("h_fadc%d",padid));
	    h->GetYaxis()->SetRangeUser(200,1000);
	    gPad->SetTopMargin(0.08);
	    h->Draw("colz");
	    h->Write();
	  }
	  count++;
	}
      }
    }
  }

  c1->Print(outpdf.c_str());
  auto poly = (TH2Poly*)file->Get("poly_rms");
  poly->Write();
  c1->Print((outpdf + ")").c_str());
  f->Close();
  
  
}
