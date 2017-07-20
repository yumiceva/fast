//////////////////////////////////////////////////////////
//
//
// Francisco Yumiceva
// Florida Institute of Technology
/////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <iterator>

#include "parser.h"
#include "CutHLT.h"
#include "CutMuon.h"
#include "CutElectron.h"
#include "CutJet.h"
#include "CutDeltaR.h"
#include "CutPhoton.h"
#include "PUReweight.h"

#include <boost/date_time/posix_time/posix_time.hpp>

#include "ROOT/TDataFrame.hxx"
#include "TChain.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"

using namespace std;
using namespace boost::posix_time;
using FourVector = ROOT::Math::XYZTVector;
using FourVectors = std::vector<FourVector>;
using CylFourVector = ROOT::Math::RhoEtaPhiVector;


int main(int ac, char** av)
{

  // Time Stamp --------------------------------------------------------------------------
  //get the current time from the clock -- one second resolution
  ptime now = second_clock::local_time();
  cout << "========= BEGIN JOB at " << to_simple_string( now ) << endl;
  // end Time Stamp ----------------------------------------------------------------------

  // Command line parser
  parser parser( ac, av);
  if ( parser.Exit() ) return 1;

  auto treeName = "ggNtuplizer/EventTree";
  //auto fileName = "/uscms_data/d2/dnoonan/13TeV_TTGamma/TTGamma_SingleLeptFromTbar_100k.root";
  auto inputfilenames = parser.InputFiles();
  auto fileName = inputfilenames[0];
  //auto outputfileName = "histos_100k.root";
  //auto fileName = "root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma_SingleLeptFromTbar.root";
  //auto outputfileName = "histos_ttgamma.root";

  if ( parser.MT() )
    ROOT::EnableImplicitMT();

  TChain chain(treeName);
  chain.Add(fileName.c_str());

  ROOT::Experimental::TDataFrame dataFrame( chain );

  // HLT
  //auto cutHLT = [](ULong64_t t) { return t >> 19 & 1 || t >> 20 & 1; };
  auto dt_HLT = dataFrame.Filter(cutHLT, {"HLTEleMuX"} );
  auto entries_HLT = dt_HLT.Count();
  cout << "HLT = " << *entries_HLT << endl;

  // PV
  auto dt_PV = dt_HLT.Filter( "isPVGood" );
  auto entries_PV = dt_PV.Count();
  cout << "PV = " << *entries_PV << endl;

  // Tight muons
  auto dt_tight_muons = dt_PV
    .Define( "relIso", calcRelIso, {"muPt","muPFChIso","muPFNeuIso","muPFPhoIso","muPFPUIso"})
    .Define( "tightmuons", getTightMuons, {"muType", "muPt","muEta", "muPhi","relIso","muChi2NDF","muTrkLayers","muMuonHits","muD0","muDz","muPixelHits","muStations"} );

  auto dt_one_T_muon = dt_tight_muons.Filter( [](const FourVectors &list) { return list.size() == 1;}, {"tightmuons"} );
  auto entries_1TMuon = dt_one_T_muon.Count();
  cout << "one tight muon = " << *entries_1TMuon << endl;
  // Veto loose muons
  auto dt_loose_muons = dt_one_T_muon.Define( "loosemuons", getLooseMuons, {"muType", "muPt","muEta", "muPhi","relIso"} );
  auto dt_veto_L_muons = dt_loose_muons.Filter( [](const FourVectors &list) { return list.size() == 1;}, {"loosemuons"} );
  auto entries_vetoLmuons = dt_veto_L_muons.Count();
  cout << "veto loose muons = " << *entries_vetoLmuons << endl;

  // Veto electrons
  auto dt_electrons = dt_veto_L_muons.Define( "vetoelectrons", getVetoElectrons, {"elePt","eleEta","eleSCEta","elePhi","eleIDbit","eleD0","eleDz"} );
  auto dt_veto_electrons = dt_electrons.Filter( [](const FourVectors &list) { return list.size() == 0;}, {"vetoelectrons"} );
  auto entries_vetoElectrons = dt_veto_electrons.Count();
  cout << "veto electrons = " << *entries_vetoElectrons << endl;
  // Jet selection
  auto dt_jets_raw = dt_veto_electrons.Define( "jets_raw", getTightJets, {"jetID","jetPt","jetEta","jetPhi","jetEn"} );
  //auto tmp = dt_jets_raw.Filter( [](const FourVectors &list) { return list.size() >= 3;}, {"jets_raw"} ).Count();
  //cout << "gte3j no deltaR = " << *tmp << endl;

  // DeltaR(jet, muon) cut
  auto dt_jets_no_muons = dt_jets_raw.Define( "jets_no_muons", listDeltaR4, {"jets_raw","tightmuons"} ); //drop jets
  // get medium photons
  auto dt_jets_gamma = dt_jets_no_muons.Define( "photons", getMphotons, {"phoEt","phoEta","phoSCEta","phoPhi","phoIDbit","phohasPixelSeed"} );
  // DeltaR(jet, gamma) cut
  auto dt_jets = dt_jets_gamma.Define( "jets", listDeltaR1, {"jets_no_muons","photons"} ); //drop jets

  auto dt_gte3j = dt_jets.Filter( [](const FourVectors &list) { return list.size() >= 3;}, {"jets"} );
  auto entries_gte3j = dt_gte3j.Count();
  cout << "gte3j = " << *entries_gte3j << endl;

  // Write a tree
  //cout << "Writing the skimmed tree ... ";
  //dt_gte3j.Snapshot(treeName,"tiny.root");
  //cout << "done" << endl;
  
  auto dt_gte3jgte1b = dt_gte3j.Define("btags", getMbtags, {"jetCSV2BJetTags"} )
    .Filter([](int btags) { return btags >= 1;}, {"btags"});
  auto entries_gte3jgte1b = dt_gte3jgte1b.Count();
  cout << "gte3jgte1b = " << *entries_gte3jgte1b <<endl;

  
  // HISTOGRAMS
  auto h_pv = dt_gte3jgte1b.Histo1D(TH1D("PV","Vertices",50,0,50), "nVtx");
  
  auto outputfileName = parser.OutputPath() + "histos_"+ parser.JobName() + ".root";
  cout << "Writing histograms to " << outputfileName;
  TFile *fFile = TFile::Open( outputfileName.c_str(), "RECREATE");
  h_pv->Write();
  cout << " done." << endl;


  cout << "========== JOB END at " << to_simple_string( second_clock::local_time() ) << endl;
  return 0;
}
