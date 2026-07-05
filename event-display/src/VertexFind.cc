#include "VertexFind.hh"

VertexFind::VertexFind(){
}

TVector3
VertexFind::VertexPointHelix(const Double_t par1[7], const Double_t par2[7], Double_t& dist)
{
  //[0 ~ 2] : center position, [3] : radius, [4] : slope
  //helix function 1
  //x = [0] + [3]*sin(alpha)
  //y = [1] + alpha*[4]
  //z = [2] + [3]*cos(alpha)

  //helix funciton 2
  //x = [5] + [8]*sin(alpha)
  //y = [6] + alpha*[9]
  //z = [7] + [8]*cos(alpha)

  //originally 10.0
  double start_alpha1 = par1[5];
  double end_alpha1 = par1[6];

  double start_alpha2 = par2[5];
  double end_alpha2 = par2[6];
    
  static TF2 fvert_helix("fvert_helix","pow(([0]+[3]*sin(x))-([5]+[8]*sin(y)),2)+pow(([1]+x*[4])-([6]+y*[9]),2)+pow(([2]+[3]*cos(x))-([7]+[8]*cos(y)),2)",start_alpha1,end_alpha1,start_alpha2,end_alpha2);

  

  fvert_helix.SetParameter(0, par1[0]);
  fvert_helix.SetParameter(1, par1[1]);
  fvert_helix.SetParameter(2, par1[2]);
  fvert_helix.SetParameter(3, par1[3]);
  fvert_helix.SetParameter(4, par1[4]);
  fvert_helix.SetParameter(5, par2[0]);
  fvert_helix.SetParameter(6, par2[1]);
  fvert_helix.SetParameter(7, par2[2]);
  fvert_helix.SetParameter(8, par2[3]);
  fvert_helix.SetParameter(9, par2[4]);


  
  Double_t alpha1, alpha2;
  fvert_helix.GetMinimumXY(alpha1, alpha2);


  Double_t x1 = par1[0] + par1[3]*sin(alpha1);
  Double_t x2 = par2[0] + par2[3]*sin(alpha2);

  Double_t y1 = par1[1] + alpha1*par1[4];
  Double_t y2 = par2[1] + alpha2*par2[4];

  Double_t z1 = par1[2] + par1[3]*cos(alpha1);
  Double_t z2 = par2[2] + par2[3]*cos(alpha2);

  Double_t vx = (x1 + x2)/2.;
  Double_t vy = (y1 + y2)/2.;
  Double_t vz = (z1 + z2)/2.;

  TVector3 vertex;
  vertex.SetX(vx);
  vertex.SetY(vy);
  vertex.SetZ(vz);

  dist = TMath::Sqrt(fvert_helix.GetMinimum());

  return vertex;
}

//get momentum direction
TVector3 VertexFind::GetHelixTangentDirection(const TVector3& point, double par[7]) {

  TVector3 center;
  center.SetX(par[0]);
  center.SetY(par[1]);
  center.SetZ(par[2]);

  double radius = par[3];
  double slope = par[4];
  
  
  // Compute alpha from the given point
  double dx = point.X() - center.X();
  double dz = point.Z() - center.Z();
  double alpha = atan2(dx, dz); // Use atan2 to calculate the angle alpha

  // Calculate derivatives with respect to alpha (direction vector components)

  double dx_dalpha = radius * cos(alpha);
  double dy_dalpha = slope;
  double dz_dalpha = -radius * sin(alpha);
  
  // Create tangent direction vector from the derivatives

  if(dz_dalpha<=0){
    dx_dalpha*=-1.;
    dy_dalpha*=-1.;
    dz_dalpha*=-1.;
  }
  TVector3 tangentDirection(dx_dalpha, dy_dalpha, dz_dalpha);

  return tangentDirection;
}
