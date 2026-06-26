struct P2 {
  double x;
  double y;
};

double Cross(const P2& a, const P2& b)
{
  return a.x*b.y - a.y*b.x;
}

P2 Sub(const P2& a, const P2& b)
{
  return {a.x-b.x, a.y-b.y};
}

// infinite line and segment intersection
bool LineSegIntersection(
    const P2& p0,
    const P2& dir,
    const P2& a,
    const P2& b,
    double& t)
{
  P2 v = Sub(b,a);
  P2 w = Sub(a,p0);

  double denom = Cross(dir,v);

  if (std::abs(denom) < 1e-12) return false;

  t = Cross(w,v)/denom;
  double u = Cross(w,dir)/denom;

  return (u >= -1e-9 && u <= 1.+1e-9);
}

double PathLength(
    const std::vector<P2>& poly,
    const P2& p0,
    const P2& dir)
{
  std::vector<double> ts;

  int n = poly.size();

  for(int i=0;i<n;i++) {

    double t;

    const P2& a = poly[i];
    const P2& b = poly[(i+1)%n];

    if(LineSegIntersection(p0,dir,a,b,t))
      ts.push_back(t);
  }

  if(ts.size()<2) return 0.;

  std::sort(ts.begin(),ts.end());

  return std::abs(ts.back()-ts.front());
}

void path_length()
{
  // ==========================================
  // pad shape
  // ==========================================

  // width : phi direction
  // height: radial direction

  //inner layer
  /*
  double width_top    = 2.52;   // mm
  double width_bottom = 1.34;   // mm
  double height       = 9.0;  // mm
  */

  
  //outer layer
  
  double width_top    = 3.0;   // mm
  double width_bottom = 2.7;   // mm
  double height       = 12.5;  // mm


  // trapezoid centered at (0,0)

  std::vector<P2> pad;

  pad.push_back({-width_bottom/2., -height/2.});
  pad.push_back({ width_bottom/2., -height/2.});
  pad.push_back({ width_top/2.,     height/2.});
  pad.push_back({-width_top/2.,     height/2.});

  // track passes through center
  P2 p0 = {0.,0.};

  TGraph* gr = new TGraph();

  int ip = 0;

  for(double angle_deg=-85; angle_deg<=85; angle_deg+=1.) {

    // 0 degree = vertical
    double th = angle_deg*TMath::DegToRad();

    P2 dir = {
      TMath::Sin(th),
      TMath::Cos(th)
    };

    double len = PathLength(pad,p0,dir);

    gr->SetPoint(ip,angle_deg*TMath::DegToRad(),len);

    ip++;
  }


  TCanvas* c1 = new TCanvas("c1","path length",800,600);

  gr->SetTitle("Path length in one pad;#alpha [rad];Path length in pad [mm]");

  gr->SetLineWidth(2);
  gr->SetLineColor(kBlue+1);
  gr->GetXaxis()->SetLimits(-1.5,1.5);

  gr->Draw("AL");

  // ==========================================
  // draw example pad
  // ==========================================

  TCanvas* c2 = new TCanvas("c2","pad geometry",600,800);

  TH2D* frame = new TH2D(
      "frame","Example pad;X [mm];Y [mm]",
      100,-10,10,
      100,-10,10);

  frame->Draw();

  double px[5], py[5];

  for(int i=0;i<4;i++) {
    px[i] = pad[i].x;
    py[i] = pad[i].y;
  }

  px[4]=px[0];
  py[4]=py[0];

  TPolyLine* pl = new TPolyLine(5,px,py);

  pl->SetLineWidth(2);
  pl->SetLineColor(kBlue);

  pl->Draw("same");

  // example track
  double ang = 10.*TMath::DegToRad();

  double x1 = -20.*TMath::Sin(ang);
  double y1 = -20.*TMath::Cos(ang);

  double x2 =  20.*TMath::Sin(ang);
  double y2 =  20.*TMath::Cos(ang);

  TLine* tr = new TLine(x1,y1,x2,y2);

  tr->SetLineColor(kRed);
  tr->SetLineWidth(3);

  tr->Draw("same");
}
