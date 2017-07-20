#include "ROOT/TDataFrame.hxx"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"

using FourVector = ROOT::Math::XYZTVector;
using FourVectors = std::vector<FourVector>;
using CylFourVector = ROOT::Math::RhoEtaPhiVector;

auto getMphotons = [](const vector<float> phoEt, const vector<float> phoEta, const vector<float> phoSCEta,
                      const vector<float> phoPhi, const vector<unsigned short> phoIDbit, const vector<int> phohasPixelSeed) {
  FourVectors list;
  for (unsigned int t=0; t< phoEt.size(); ++t ) {
    float absSCEta = fabs( phoSCEta[t] );
    bool passEtaOverlap = (absSCEta < 1.4442) || (absSCEta > 1.566);
    bool passMediumPhotonID = phoIDbit[t] >> 1 & 1;
    if ( phoEt[t] > 15.0 &&
         fabs( phoEta[t] ) < 2.5 &&
         passEtaOverlap &&
         passMediumPhotonID &&
         !phohasPixelSeed[t] ) {
      CylFourVector cp4( phoEt[t], phoEta[t], phoPhi[t]);
      double M = 0;
      double E = sqrt(cp4.R() * cp4.R() + M * M);
      FourVector p4( cp4.X(), cp4.Y(), cp4.Z(), E );
      list.emplace_back( p4 );
    }
  }
  return list;
};

