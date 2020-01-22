import os

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'deepJetTree.root'
options.inputFiles = '/store/mc/RunIIFall17MiniAODv2/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-7_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/60000/20AD20B4-B563-E811-86C5-1866DAEA6BC4.root'
options.maxEvents = -1

options.parseArguments()

#####################
### setup process ###
#####################

process = cms.Process("MuonNtuple")

process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
)

#################
### GlobalTag ###
#################
envvar = 'mcgt'
from Configuration.AlCa.GlobalTag import GlobalTag
GT = {'mcgt': '94X_mc2017_realistic_v17', 'datagt': '94X_dataRun2_v11'}
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
# https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
process.GlobalTag = GlobalTag(process.GlobalTag, GT[envvar], '')

#############################
### Setup rest of running ###
#############################
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(options.outputFile),
)

process.schedule = cms.Schedule()

process.deepJetTree = cms.EDAnalyzer('DeepJetTree',
    jetSrc = cms.InputTag('slimmedJets'),
    genSrc  = cms.InputTag('prunedGenParticles'),
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
)
process.deepJetTreePath = cms.Path()
process.deepJetTreePath += process.deepJetTree
process.schedule.append(process.deepJetTreePath)
