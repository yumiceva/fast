#include "ROOT/TDataFrame.hxx"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"

using FourVector = ROOT::Math::XYZTVector;
using FourVectors = std::vector<FourVector>;
using CylFourVector = ROOT::Math::RhoEtaPhiVector;


auto listDeltaR1 = [](const FourVectors &list, const FourVectors &list2) {
  FourVectors outlist;
  if ( list2.size() == 0 ) return list;
  for (auto &t : list) {
    bool keep = false;
    for (auto &t2 : list2) {
      float deltaR = sqrt( pow(t.Eta() - t2.Eta(),2) + pow(t.Phi() - t2.Phi(),2) );
      if ( deltaR >= 0.1 ) keep = true;
    }
    if (keep) outlist.emplace_back( t );
  }
  return outlist;
};

auto listDeltaR4 = [](const FourVectors &list, const FourVectors &list2) {
  FourVectors outlist;
  if ( list2.size() == 0 ) return list;
  for (auto &t : list) {
    bool keep = false;
    for (auto &t2 : list2) {
      float deltaR = sqrt( pow(t.Eta() - t2.Eta(),2) + pow(t.Phi() - t2.Phi(),2) );
      if ( deltaR >= 0.4 ) keep = true;
    }
    if (keep) outlist.emplace_back( t );
  }
  return outlist;
};

auto getminDeltaR = [](const FourVectors &list, const FourVectors &list2) {
  // minimum deltaR between a list of objects and the leading object of list 2
  float deltaR = 1000.;
  FourVector p4 = list2[0];
  for (auto &t : list) {
    float tmp = sqrt( pow(t.Eta() - p4.Eta(),2) + pow(t.Phi() - p4.Phi(),2) );
    if ( tmp <= deltaR ) deltaR = tmp;
  }
  return deltaR;
};

