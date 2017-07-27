#include "ROOT/TDataFrame.hxx"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"

using FourVector = ROOT::Math::XYZTVector;
using FourVectors = std::vector<FourVector>;
using CylFourVector = ROOT::Math::RhoEtaPhiVector;

//bool operator<( const FourVector& tlv1, const FourVector& tlv2 ) { return tlv1.Pt() < tlv2.Pt(); }

auto calcM3 = [](const FourVectors &list) {
  // Compute M3
  float M3 = -1;
  float sumpt = -1;
  int icomb = 0;
  do {
    FourVector tlv_M3 = list[0] + list[1] + list[2];
    if ( tlv_M3.Pt() >= sumpt )
      {
        sumpt = tlv_M3.Pt();
        M3 = tlv_M3.M();
      }
    icomb++;

  } while ( next_permutation(list.begin(), list.end() ) );
  
  return M3;

};
