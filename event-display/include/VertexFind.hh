#ifndef VERTEXFIND_HH
#define VERTEXFIND_HH

#include <iostream>

#include <TVector3.h>
#include <TMath.h>
#include <TF2.h>


class VertexFind{
public:
  VertexFind();
  TVector3 VertexPointHelix(const Double_t helix_par1[7], const Double_t helix_par2[7], Double_t &dist);
  TVector3 GetHelixTangentDirection(const TVector3 &vertex_point, Double_t helix_par[7]);
};

#endif
