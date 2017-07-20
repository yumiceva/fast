#include "ROOT/TDataFrame.hxx"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"

using FourVector = ROOT::Math::XYZTVector;
using FourVectors = std::vector<FourVector>;
using CylFourVector = ROOT::Math::RhoEtaPhiVector;

auto getVetoElectrons = [](const vector<float> elePt, const vector<float> eleEta, const vector<float> eleSCEta, const vector<float> elePhi, const vector<unsigned short> eleIDbit,
                           const vector<float> eleD0, const vector<float> eleDz) {
  FourVectors list;
  float cutD0 = 0.05; //barrel
  float cutDz = 0.1; // barrel
  for (unsigned int t=0; t< elePt.size(); ++t ) {
    bool cutbasedID_veto = ( eleIDbit[t] >> 0 & 1);
    bool passEtaEBEEGap = (fabs(eleSCEta[t]) < 1.4442) || (fabs(eleSCEta[t]) > 1.566);

    if ( fabs(eleEta[t]) > 1.479 ) {
      // endcap
      cutD0 = 0.1;
      cutDz = 0.2;
    }
    if ( elePt[t] > 15.0 &&
         fabs(eleEta[t]) < 2.5 &&
         cutbasedID_veto &&
         passEtaEBEEGap &&
         eleD0[t] < cutD0 &&
         eleDz[t] < cutDz ) {

      CylFourVector cp4( elePt[t], eleEta[t], elePhi[t]);
      double M = 0;
      double E = sqrt(cp4.R() * cp4.R() + M * M);
      FourVector p4( cp4.X(), cp4.Y(), cp4.Z(), E );
      list.emplace_back( p4 );
    }
  }
  return list;

};

