import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import re

options = VarParsing('analysis')
options.register(
    "jobNumber",
    -1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Job index"
)
options.parseArguments()

process = cms.Process('magic')

inputFile = options.inputFiles[0]
match = re.match(".*((?:Data|MC)13TeV[^/]*)", inputFile)
if match:
    outputFile = match.group(1)
    if options.jobNumber > -1:
        outputFile += '_%d' % options.jobNumber
    outputFile += '.root'
else:
    outputFile = inputFile.replace('.root','_out.root')

process.config = cms.PSet(
    input = cms.string(inputFile),
    output = cms.string(outputFile),
    dirName = cms.string("mainAnalyzer/data"),
    runSystematics = cms.bool(True),
    runOptimization = cms.bool(True),
    evStart = cms.int32(0),
    evEnd = cms.int32(-1),
    BtagSF = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/CSVv2.csv"),
    jesUncFileName = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/Summer15_25nsV6_MC_Uncertainty_AK4PFchs.txt"),
    pdfInput = cms.string("root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/monoZ/llvvNtuple_25ns_PDFs/"),
    PU_Central = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/PU_Central.root"),
    PU_Up = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/PU_minBiasUP.root"),
    PU_Down = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/PU_minBiasDOWN.root"),
)
