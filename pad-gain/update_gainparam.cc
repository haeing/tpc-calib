#include "../TPCPadHelper_260416.hh"

void update_gainparam(
    const char* param_in  = "param_history/TPCParam_e72_20260626",
    const char* rootfile  = "result/260629/padgain-fillpad-260629.root",
    const char* histname  = "TPC_gain_all",
    const char* param_out = "param_history/TPCParam_e72_20260629")
{
  TFile *f = new TFile(rootfile, "READ");
  TH2Poly *hGain = (TH2Poly*)f->Get(histname);

  if(!hGain){
    cout << "Cannot find histogram : " << histname << endl;
    return;
  }

  ifstream fin(param_in);
  ofstream fout(param_out);

  string line;
  int iline = 0;
  int ipad_count = 0;

  while(getline(fin, line)){
    iline++;

    // TPC ADC data starts from line 13
    if(iline >= 13 && ipad_count < NumOfPadTPC){

      stringstream ss(line);

      int layer, row, aty;
      double p0, p1;

      if(ss >> layer >> row >> aty >> p0 >> p1){

        if(aty == 0){
          int padid = tpc::GetPadId(layer, row);
          double mpv = hGain->GetBinContent(padid + 1);

          if(mpv > 0.){
            p0 = 0.;
            p1 *= 200. / mpv;
          }
        }

        fout
          << setw(4) << layer << " "
          << setw(4) << row   << " "
          << setw(4) << aty   << " "
          << fixed << setprecision(8)
          << setw(14) << p0 << " "
          << setw(14) << p1
          << endl;

        ipad_count++;
        continue;
      }
    }

    // everything else: copy original line exactly
    fout << line << endl;
  }

  cout << "Updated pads = " << ipad_count << endl;
  cout << "Saved : " << param_out << endl;
}
