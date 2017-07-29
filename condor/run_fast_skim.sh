#!/bin/bash

job=$1

if [ -z ${_CONDOR_SCRATCH_DIR} ] ; then 
	echo "Running Interactively"
else
	echo "Running In Batch"
	cd ${_CONDOR_SCRATCH_DIR}
	echo ${_CONDOR_SCRATCH_DIR}
	echo "tar -xzf CMSSW.tgz"
	tar -xzf CMSSW.tgz
	cd CMSSW_9_3_ROOT6_X_2017-07-27-2300/src/
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	cd  fast/
fi

eval `scramv1 runtime -sh`



outputdir="root://cmseos.fnal.gov//store/user/yumiceva/ttgamma13TeV/skims"


files=("TTGamma_SingleLeptFromTbar_" \
"TTGamma_SingleLeptFromT_" \
"TTGamma_Dilepton_" \
"TTGamma_Hadronic_" \
"TTbar_" \
"W1jets_" \
"W2jets_" \
"W3jets_" \
"W4jets_" \
"DYjets_" \
"ST_s-channel_" \
"ST_t-channel_" \
"ST_tbar-channel_" \
"ST_tW-channel_" \
"ST_tbarW-channel_" \
"TTW_" \
"TTZ_" \
"WGamma_" \
"ZGamma_" \
"QCD_20to30EM_" \
"QCD_30to50EM_" \
"QCD_50to80EM_" \
"QCD_80to120EM_" \
"QCD_120to170EM_" \
"QCD_170to300EM_" \
"QCD_300toInfEM_" \
"QCD_20to30Mu_" \
"QCD_30to50Mu_" \
"QCD_50to80Mu_" \
"QCD_80to120Mu_" \
"QCD_120to170Mu_" \
"QCD_170to300Mu_" \
"QCD_300to470Mu_" \
"QCD_470to600Mu_" \
"QCD_600to800Mu_" \
"QCD_800to1000Mu_" \
"QCD_1000toInfMu_" \
"Data_SingleMu_b_" \
"Data_SingleMu_c_" \
"Data_SingleMu_d_" \
"Data_SingleMu_e_" \
"Data_SingleMu_f_" \
"Data_SingleMu_g_" \
"Data_SingleMu_h_")


inputfiles=("root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma_SingleLeptFromTbar.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma_semileptfromT.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma_dilept.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ttgamma_hadronic.root" \
"root://cmseos.fnal.gov//store/user/yumiceva/ntuples_2016/TTbar.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/W1jets.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/W2jets.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/W3jets.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/W4jets.root" \
"root://cmseos.fnal.gov//store/user/yumiceva/ntuples_2016/ntuple_DYjets.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ST_s.root" \
"root://cmseos.fnal.gov//store/user/yumiceva/ntuples_2016/ST_t_top.root" \
"root://cmseos.fnal.gov//store/user/yumiceva/ntuples_2016/ST_t_bar.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ST_tW_top.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ST_t_tWbar.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/TTWJets.root" \
"root://cmseos.fnal.gov//store/user/yumiceva/ntuples_2016/TTZ.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/WGToLNuG_aMCatNLO_summer16.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/ZGTo2LG_aMCatNLO_summer16.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_20to30EM.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_30to50EM.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_50to80EM.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_80to120EM.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_120to170EM.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_170to300EM.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_300toInfEM.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_20to30Mu.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_30to50Mu.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_50to80Mu.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_80to120Mu.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_120to170Mu.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_170to300Mu.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_300to470Mu.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_470to600Mu.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_600to800Mu.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_800to1000Mu.root" \
"root://cmseos.fnal.gov//store/user/troy2012/ntuples_2016/QCD_1000toInfMu.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016B_FebReminiAOD.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016C_FebReminiAOD.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016D_FebReminiAOD.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016E_FebReminiAOD.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016F_FebReminiAOD1.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016F_FebReminiAOD2.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016G_FebReminiAOD.root" \
"root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016H_FebReminiAODv2.root root://cmseos.fnal.gov//store/user/dnoonan/13TeV_ggNTuples/job_SingleMu_Run2016H_FebReminiAODv3.root")


echo "./fast -i ${inputfiles[job]} -j ${files[job]} -s skim_${files[job]}.root"

./fast -i ${inputfiles[job]} -j ${files[job]} -s skim_${files[job]}.root

echo "xrdcp -f skim_${files[job]}.root ${outputdir}/."

xrdcp -f skim_${files[job]}.root ${outputdir}/.
