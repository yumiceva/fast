#include "TCanvas.h"
#include "TChain.h"
#include "TH1F.h"
#include "TTree.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"
#include "ROOT/TDataFrame.hxx"

#include "PUReweight.h"

using FourVector = ROOT::Math::XYZTVector;
using FourVectors = std::vector<FourVector>;
using CylFourVector = ROOT::Math::RhoEtaPhiVector;


auto cutHLT = [](ULong64_t t) { return t >> 19 & 1 || t >> 20 & 1; };
auto cutPV  = [](Bool_t t) { return t; };
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

  float deltaR = 1000.;
  FourVector p4 = list2[0];
  for (auto &t : list) {
    float tmp = sqrt( pow(t.Eta() - p4.Eta(),2) + pow(t.Phi() - p4.Phi(),2) );
    if ( tmp <= deltaR ) deltaR = tmp;
  }
  return deltaR;
};

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

void tdf()
{
  auto treeName = "ggNtuplizer/EventTree";
  auto fileName = "/uscms_data/d2/dnoonan/13TeV_TTGamma/TTGamma_SingleLeptFromTbar_100k.root";
  auto outputfileName = "histos_100k.root";
  //auto fileName = "root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma_SingleLeptFromTbar.root";
  //auto outputfileName = "histos_ttgamma.root";

  //ROOT::EnableImplicitMT();

  TChain chain(treeName);
  chain.Add(fileName);

  ROOT::Experimental::TDataFrame dataFrame( chain );
  
  // PU reweight
  //==============
  vector<string> puFileNames;
  puFileNames.emplace_back( fileName );
  PUReweight *PUweighter = new PUReweight(1, puFileNames, "Data_2016BCDGH_Pileup.root");

  auto getPUweight = [PUweighter] ( int nPUInfo, vector<int> puBX, vector<float> puTrue ) {
    float w = PUweighter->getWeight( nPUInfo, &puBX, &puTrue );
    if ( w <= 0 ) cout << "PU w = " << w <<endl;
    return w;
  };
  //==============

  // HLT
  auto dt_HLT = dataFrame.Filter(cutHLT, {"HLTEleMuX"} );
  auto entries_HLT = dt_HLT.Count();
  cout << "HLT = " << *entries_HLT << endl;
  // PV
  auto dt_PV = dt_HLT.Filter( cutPV, {"isPVGood"} );
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

  auto dt_deltaR_muj = dt_gte3jgte1b.Define("deltaRmuj", getminDeltaR, {"jets","tightmuons"} );
  auto h_deltaR = dt_deltaR_muj.Histo1D(TH1D("deltaR","min deltaR(j,mu)",50,0,5), "deltaRmuj" );
  
  // PU weight
  auto dt_gte3jgte1b_PU = dt_gte3jgte1b.Define( "PUweight", getPUweight, {"nPUInfo","puBX", "puTrue"} );
 
  auto h_pvPU = dt_gte3jgte1b_PU.Define("nVtxPU", [](float weight, int nVtx){ return nVtx*weight;}, {"PUweight","nVtx"} )
    .Histo1D(TH1D("PVPU","Vertices PU",50,0,50), "nVtxPU");
  
  //TCanvas c1pu;
  //h_pvPU->Draw();
  //c1pu.Print("pvpu.png");

  auto h_pv = dt_gte3jgte1b.Histo1D(TH1D("PV","Vertices",50,0,50), "nVtx");
  //TCanvas c1;
  //h_pv->Draw();
  //c1.Print("pv.png");
  
  auto h_mu_pt = dt_gte3jgte1b.Define("themuPt", [](const FourVectors &list) { return list[0].Pt(); },{"tightmuons"} )
    .Histo1D(TH1D("muPt","p_{T}",100,30,400), "themuPt");
  //TCanvas c2;
  //h_mu_pt->Draw();
  //c2.Print("mu_pt.png");
  cout << "entries in mu_pt " << h_mu_pt->GetEntries() << endl;

  TFile *fFile = TFile::Open( outputfileName, "RECREATE");
  cout << "Output root file: " << outputfileName << endl;
  
  h_deltaR->Write();
  h_pv->Write();
  h_pvPU->Write();
  h_mu_pt->Write();

  fFile->Close();
  //auto Muons = 

}

