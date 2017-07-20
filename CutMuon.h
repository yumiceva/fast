#include "ROOT/TDataFrame.hxx"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"

using FourVector = ROOT::Math::XYZTVector;
using FourVectors = std::vector<FourVector>;
using CylFourVector = ROOT::Math::RhoEtaPhiVector;


auto getLooseMuons = [](const vector<Int_t> muType, const vector<float> muPt, const vector<float> muEta, const vector<float> muPhi, const vector<float> relIso) {
  FourVectors list;
  for (unsigned int t=0; t< muType.size(); ++t ) {
    if ( (muType[t] & (1<<1) || muType[t] & (1<<2) ) &&
         muType[t] & (1<<5) &&
         muPt[t] > 15 &&
         fabs(muEta[t]) < 2.4 &&
         relIso[t] < 0.25 ) {
      CylFourVector cp4( muPt[t], muEta[t], muPhi[t]);
      double M = 0;
      double E = sqrt(cp4.R() * cp4.R() + M * M);
      FourVector p4( cp4.X(), cp4.Y(), cp4.Z(), E );
      list.emplace_back( p4 );
    }
  }
  return list;
};

auto getTightMuons = [](const vector<Int_t> muType, const vector<float> muPt, const vector<float> muEta, const vector<float> muPhi, const vector<float> relIso,
                        const vector<float> muChi2NDF, const vector<int> muTrkLayers, const vector<int> muMuonHits, const vector<float> muD0,
                        const vector<float> muDz, const vector<int> muPixelHits, const vector<int> muStations) {
  FourVectors list;
  for (unsigned int t=0; t< muType.size(); ++t ) {
    if ( muType[t] & (1<<1) && muType[t] & (1<<2) &&
         muType[t] & (1<<5) &&
         muPt[t] > 30 &&
         fabs(muEta[t]) < 2.4 &&
         relIso[t] < 0.15 &&
         muChi2NDF[t] < 10 &&
         muTrkLayers[t] > 5 &&
         muMuonHits[t] > 0 &&
         muD0[t] < 0.2 &&
         fabs( muDz[t] ) < 0.5 &&
         muPixelHits[t] > 0 &&
         muStations[t] > 1) {

      CylFourVector cp4( muPt[t], muEta[t], muPhi[t]);
      double M = 0;
      double E = sqrt(cp4.R() * cp4.R() + M * M);
      FourVector p4( cp4.X(), cp4.Y(), cp4.Z(), E );
      list.emplace_back( p4 );
    }
  }
  return list;
};

auto calcRelIso = [](const vector<float> muPt, const vector<float> muPFChIso, const vector<float> muPFNeuIso, const vector<float> muPFPhoIso, const vector<float> muPFPUIso) {

  vector<float> relIsoVec;
  for (unsigned int t=0; t< muPt.size(); ++t ) {
    float relIso = ( muPFChIso[t] + fmax(0.0, muPFNeuIso[t] + muPFPhoIso[t] -0.5*muPFPUIso[t] ) ) / muPt[t];
    relIsoVec.emplace_back( relIso );
  }
  return relIsoVec;
};

