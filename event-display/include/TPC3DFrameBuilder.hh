#ifndef TPC_3D_FRAME_BUILDER_HH
#define TPC_3D_FRAME_BUILDER_HH

#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGeoShape.h>

class TPC3DFrameBuilder {
public :
  TPC3DFrameBuilder();
  TGeoVolume* TPC3DGeometry();

};

#endif
