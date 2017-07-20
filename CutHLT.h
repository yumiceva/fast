
#include "ROOT/TDataFrame.hxx"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"

auto cutHLT = [](ULong64_t t) { return t >> 19 & 1 || t >> 20 & 1; };

