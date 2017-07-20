#include "ROOT/TDataFrame.hxx"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"

using FourVector = ROOT::Math::XYZTVector;
using FourVectors = std::vector<FourVector>;
using CylFourVector = ROOT::Math::RhoEtaPhiVector;


auto getTightJets = [](const vector<int> jetID, const vector<float> jetPt, const vector<float> jetEta, const vector<float> jetPhi, const vector<float> jetEn) {
  FourVectors list;
  for (unsigned int t=0; t< jetPt.size(); ++t ) {
    if ( (jetID[t] >>2 & 1) &&
         jetPt[t] > 30.0 &&
         fabs(jetEta[t]) < 2.4 ) {
      CylFourVector cp4( jetPt[t], jetEta[t], jetPhi[t]);
      FourVector p4( cp4.X(), cp4.Y(), cp4.Z(), jetEn[t] );
      list.emplace_back( p4 );
    }
  }
  return list;
};

auto getMbtags = [](const vector<float> jetCSV2BJetTags) {
  int btags = 0;
  for (unsigned int t=0; t< jetCSV2BJetTags.size(); ++t ) {

    if ( jetCSV2BJetTags[t] > 0.8484 ) btags++;
  }
  return btags;
};

