#include "TH1D.h"

// PV
TH1D h_pv_noPU = TH1D("pv_noPU","Vertices", 50, 0, 50);
TH1D h_pv = TH1D("pv","Vertices", 50,0, 50);

// Muon
TH1D h_mu1_pt = TH1D("h_mu1_pt","muon p_{T} [GeV]", 150, 30, 300);
TH1D h_mu1_eta = TH1D("h_mu1_eta","muon eta", 100, -2.4, 2.4);

// Electron

// Jet
TH1D h_jet1_pt = TH1D("h_jet1_pt","leading jet p_{T} [GeV]", 150, 30, 500);

// DeltaR
TH1D h_deltaR_phomu = TH1D("h_deltaR_phomu","{#Delta}R(#gamma,#mu)",50, 0, 7);
