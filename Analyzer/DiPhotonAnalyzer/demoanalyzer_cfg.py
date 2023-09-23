# Parts of this code are adapted from: https://github.com/cms-opendata-analyses/HiggsExample20112012/blob/master/Level4/demoanalyzer_cfg_level4data.py

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

# define JSON file for 2012 data
goodJSON = '/home/cms-opendata/CMSSW_5_3_32/src/Demo/DemoAnalyzer/Datasets/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'

myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')


files2012data = FileUtils.loadListFromFile ('/home/cms-opendata/CMSSW_5_3_32/src/Demo/DemoAnalyzer/Datasets/CMS_Run2012A_Photon_AOD_22Jan2013-v1_20000_file_index.txt')
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*files2012data    
    )   
)

# apply JSON file
#   (needs to be placed *after* the process.source input file definition!)
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)

process.demo = cms.EDAnalyzer('DemoAnalyzer'
)

process.TFileService = cms.Service("TFileService",
       fileName = cms.string('T5_test.root')
                                   )

process.p = cms.Path(process.demo)
