//////////////////////////////////////////////////////////
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
#include "CutM3.h"
#include "CutPhoton.h"
#include "Histos.h"
#include "PUReweight.h"

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>

#include "ROOT/TDataFrame.hxx"
#include "TChain.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"

using namespace std;
using namespace boost::posix_time;
using boost::format;
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
  auto inputfilenames = parser.InputFiles();
  auto fileName = inputfilenames[0];


  if ( parser.MT() )
    ROOT::EnableImplicitMT(8);

  // Load tree, for the moment one file per job
  TChain chain(treeName);
  chain.Add(fileName.c_str());

  // Construct a data frame
  ROOT::Experimental::TDataFrame dataFrame( chain );

  // event counter
  vector<string> counterlabels = {"Begin",
                        "HLT",
                        "PV",
                        "TightMuon",
                        "vetoLooseMuon",
                        "vetoElectron",
                        "nodeltaR_gte3j",
                        "gte3j"

  };
  map< string, float > counter;
  for ( vector<string>::iterator ii = counterlabels.begin(); ii != counterlabels.end(); ++ii) {
    counter[*ii] = 0.;
  }
  // Total entries
  counter["Begin"] = chain.GetEntries();

  // HLT __________________________________________________
  auto dt_HLT = dataFrame.Filter(cutHLT, {"HLTEleMuX"} );
  counter["HLT"] = *(dt_HLT.Count());
  cout << format("%-25s %12.1f\n") % "HLT" % counter["HLT"];
  

  // PV __________________________________________________
  auto dt_PV = dt_HLT.Filter( "isPVGood" );
  counter["PV"] = *(dt_PV.Count());
  cout << format("%-25s %12.1f\n") % "PV" % counter["PV"]; 


  // Tight muons _________________________________________
  auto dt_tight_muons = dt_PV
    .Define( "relIso", calcRelIso, {"muPt","muPFChIso","muPFNeuIso","muPFPhoIso","muPFPUIso"})
    .Define( "tightmuons", getTightMuons, {"muType", "muPt","muEta", "muPhi","relIso","muChi2NDF","muTrkLayers","muMuonHits","muD0","muDz","muPixelHits","muStations"} );

  auto dt_one_T_muon = dt_tight_muons.Filter( [](const FourVectors &list) { return list.size() == 1;}, {"tightmuons"} );
  counter["TightMuon"] = *(dt_one_T_muon.Count());
  cout << format("%-25s %12.1f\n") % "one tight muon" % counter["TightMuon"];


  // Veto loose muons ____________________________________
  auto dt_loose_muons = dt_one_T_muon.Define( "loosemuons", getLooseMuons, {"muType", "muPt","muEta", "muPhi","relIso"} );
  auto dt_veto_L_muons = dt_loose_muons.Filter( [](const FourVectors &list) { return list.size() == 1;}, {"loosemuons"} );
  counter["vetoLooseMuon"] = *(dt_veto_L_muons.Count());
  cout << format("%-25s %12.1f\n") % "veto loose muons" % counter["vetoLooseMuon"];


  // Veto electrons ______________________________________
  auto dt_electrons = dt_veto_L_muons.Define( "vetoelectrons", getVetoElectrons, {"elePt","eleEta","eleSCEta","elePhi","eleIDbit","eleD0","eleDz"} );
  auto dt_veto_electrons = dt_electrons.Filter( [](const FourVectors &list) { return list.size() == 0;}, {"vetoelectrons"} );
  counter["vetoElectron"] = *(dt_veto_electrons.Count());
  cout << format("%-25s %12.1f\n") % "veto electrons" % counter["vetoElectron"];


  // Jet selection ________________________________________
  auto dt_jets_raw = dt_veto_electrons.Define( "jets_raw", getTightJets, {"jetID","jetPt","jetEta","jetPhi","jetEn"} );

  auto dt_skim = dt_jets_raw.Filter( [](const FourVectors &list) { return list.size() >= 3;}, {"jets_raw"} );
  counter["nodeltaR_gte3j"] = *(dt_skim.Count());
  cout << format("%-25s %12.1f\n") % "no deltaR gte3j" % counter["nodeltaR_gte3j"];


  // Write a tree _________________________________________
  if ( parser.SkimName() != "" && false) {

    if ( parser.MT() )
      ROOT::DisableImplicitMT();

    cout << endl << "Writing a skimmed tree ... " << endl;
    auto skimfilename = parser.SkimName();
    dt_skim.Snapshot(treeName, skimfilename.c_str() );
    cout << "done." << endl;

    // write cutflow histogram in the skimmed file
    TH1F *hcutflow = new TH1F("cutflow","cutflow skim", int(counterlabels.size()), 0, int(counterlabels.size()) );
    int ibin = 1;
    for ( vector<string>::iterator ii = counterlabels.begin(); ii != counterlabels.end(); ++ii) {
      hcutflow->SetBinContent( ibin, counter[*ii] );
      hcutflow->GetXaxis()->SetBinLabel(ibin, (*ii).c_str() );
      ibin++;
    }
    TFile *skimFile = TFile::Open( skimfilename.c_str(), "UPDATE" );
    hcutflow->Write();
    skimFile->Close();
    cout << "Cutflow histogram saved." << endl;
    // stop here
    cout << "========== JOB END at " << to_simple_string( second_clock::local_time() ) << endl;
    return 0;
  }
  
  // DeltaR(jet, muon) cut __________________________
  auto dt_jets_no_muons = dt_jets_raw.Define( "jets_no_muons", listDeltaR4, {"jets_raw","tightmuons"} ); //drop jets

  // get medium photons _____________________________
  auto dt_jets_gamma = dt_jets_no_muons.Define( "photons_raw", getMphotons, {"phoEt","phoEta","phoSCEta","phoPhi","phoIDbit","phohasPixelSeed"} );

  // DeltaR(jet, gamma) cut _________________________
  auto dt_jets = dt_jets_gamma.Define( "jets", listDeltaR1, {"jets_no_muons","photons_raw"} ); //drop jets

  auto dt_gte3j = dt_jets.Filter( [](const FourVectors &list) { return list.size() >= 3;}, {"jets"} );
  counter["gte3j"] = *(dt_gte3j.Count());
  cout << format("%-25s %12.1f\n") % "gte3j" % counter["gte3j"];

  // One b tag ______________________________________
  auto dt_gte3jgte1b = dt_gte3j.Define("btags", getMbtags, {"jetCSV2BJetTags"} )
    .Filter([](int btags) { return btags >= 1;}, {"btags"});
  counter["gte3jgte1b"] = *(dt_gte3jgte1b.Count());
  cout << format("%-25s %12.1f\n") % "gte3jgte1b" % counter["gte3jgte1b"];

  // DeltaR( gamma, muon ) cut ______________________
  //auto dt_pho_nomu = dt_gte3jgte1b.Define("photons_nomuons", listDeltaR4, {"photons_raw","tightmuons"} ); // drop photons
  //auto dt_pho_nomu_noj = dt_pho_nomu.Define("photons_nomunoj", listDeltaR4, {"photons_nomuons","jets"} ); // drop photons
  auto dt_pho_nomu_noj = dt_gte3jgte1b.Define("photons_nomuons", listDeltaR4, {"photons_raw","tightmuons"} )
    .Define("photons_nomunoj", listDeltaR4, {"photons_nomuons","jets"} ); // drop photons


  if ( parser.SkimName() != "") {

    cout << endl << "Writing a mini skimmed tree ... " << endl;
    auto skimfilename = parser.SkimName();
    dt_skim.Snapshot(treeName, skimfilename.c_str(), {"nVtx","tightmuons","jets","jetCSV2BJetTags","photons_nomunoj"} );
    cout << "done." << endl;

  }

  //counter["pho_nomu"] = *(dt_pho_nomu.Count());
  //cout << format("%-25s %12.1f\n") % "photons no muons" % counter["pho_nomu"];
  auto dt_deltaR_phomu = dt_pho_nomu_noj.Define("deltaRphomu", getminDeltaR, {"photons_nomunoj","tightmuons"} ); // column of deltaR phomu
  auto dt_deltaR_phoj = dt_pho_nomu_noj.Define("deltaRphoj", getminDeltaR, {"photons_nomunoj","jets"} ); // column of deltaR phoj


  // Calculate M3
  //auto dt_M3 = dt_pho_nomu_noj.Define("M3", calcM3, {"jets"} ); // column of M3

  // HISTOGRAMS _____________________________________
  cout << "Creating histograms ..." << endl;
  //vector< ROOT::Experimental::TDF::TResultProxy<TH1D> > list_histos;
  //auto hh_pv_noPU = dt_gte3jgte1b.Histo1D( TH1D{"h","h",10,0,10}, string("nVtx") );//TH1D(h_pv_noPU), "nVtx" );
  //list_histos.emplace_back( hh_pv_noPU );

  vector< TH1D > list_histos;

  dt_gte3jgte1b.Foreach([](int b1) { h_pv_noPU.Fill(b1); }, {"nVtx"} );
  list_histos.emplace_back( h_pv_noPU );

  auto hh_mu1_pt = dt_gte3jgte1b.Define("themuPt", [](const FourVectors &list) { return list[0].Pt(); },{"tightmuons"} );
  hh_mu1_pt.Foreach([](double b1) { h_mu1_pt.Fill(b1); }, {"themuPt"} );
  list_histos.emplace_back( h_mu1_pt );

  /*
  auto hh_mu1_pt = dt_gte3jgte1b.Define("themuPt", [](const FourVectors &list) { return list[0].Pt(); },{"tightmuons"} )
    .Histo1D( TH1D(h_mu1_pt) , "themuPt");
  list_histos.emplace_back( hh_mu1_pt );
  auto hh_deltaR_phomu = dt_deltaR_phomu.Histo1D( TH1D(h_deltaR_phomu), "deltaRphomu" );
  list_histos.emplace_back( hh_deltaR_phomu );
  auto hh_deltaR_phoj = dt_deltaR_phoj.Histo1D( TH1D(h_deltaR_phoj), "deltaRphoj" );
  list_histos.emplace_back( hh_deltaR_phoj );
  auto hh_M3 = dt_m3.Histo1D( TH1D(h_M3), "M3" );
  */

  // Write histograms _______________________________
  auto outputfileName = parser.OutputPath() + "histos_"+ parser.JobName() + ".root";
  cout << "Writing histograms to " << "\033[1m" << outputfileName << "\033[0m";
  TFile *fFile = TFile::Open( outputfileName.c_str(), "RECREATE");
  for ( auto &t : list_histos ) t.Write();
  cout << " done." << endl;


  cout << "========== JOB END at " << to_simple_string( second_clock::local_time() ) << endl;
  return 0;
}
