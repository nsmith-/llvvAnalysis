llvvAnalysis
==============

General:

	git clone git@github.com:rjwang/llvvAnalysis.git


Prepare:

	setenv SCRAM_ARCH slc6_amd64_gcc491
	scramv1 project CMSSW CMSSW_7_4_14
	cd CMSSW_7_4_14/src
	cmsenv
	git cms-merge-topic ikrav:egm_id_7.4.12_v1

	## set up the RooStats-based statistics tools for Higgs PAG
	## https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit
	git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
	cd HiggsAnalysis/CombinedLimit
	git checkout 74x-root6
	git fetch origin
	git checkout v6.0.0
	scramv1 b vclean; scramv1 b # always make a clean build, as scram doesn't always see updates to src/LinkDef.h

MiniAOD PHYS14:

	source /afs/cern.ch/cms/cmsset_default.csh
	voms-proxy-init
	xrdcp root://xrootd.unl.edu//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root /tmp/`whoami`/
