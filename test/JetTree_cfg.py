import os

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'miniTree.root'
#options.inputFiles= '/store/mc/RunIIFall17MiniAOD/WZ_TuneCP5_13TeV-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/60000/1C7C7FB8-14E6-E711-90A9-0025905A60FE.root' # WZ
options.inputFiles = '/store/mc/RunIIFall17MiniAODv2/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-5_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/5E689511-6A43-E811-8ABC-00269E95B17C.root'
#options.inputFiles = '/store/data/Run2017F/SingleMuon/MINIAOD/17Nov2017-v1/010000/183BE17C-A5EC-E711-8F69-0025905A607E.root' # ReReco
options.maxEvents = -1
options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Events to skip")
options.register('reportEvery', 100, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Report every")
options.register('isMC', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Sample is MC")
options.register('isREMINIAOD', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Sample is ReMiniAOD")
options.register('reHLT', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Sample is reHLT")
options.register('runMetFilter', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Run the recommended MET filters")
options.register('crab', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Make changes needed for crab")

options.parseArguments()

#####################
### setup process ###
#####################

process = cms.Process("MiniNtuple")

process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi") # Needed by EGamma energy correction

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
)

process.RandomNumberGeneratorService = cms.Service(
    "RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    ),
    calibratedPatPhotons = cms.PSet(
        initialSeed = cms.untracked.uint32(2),
        engineName = cms.untracked.string('TRandom3')
    ),
    mRoch = cms.PSet(
        initialSeed = cms.untracked.uint32(3),
        engineName = cms.untracked.string('TRandom3')
    ),
)

#################
### GlobalTag ###
#################
envvar = 'mcgt' if options.isMC else 'datagt'
from Configuration.AlCa.GlobalTag import GlobalTag
#GT = {'mcgt': 'auto:run2_mc', 'datagt': 'auto:run2_data'}
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
# https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
GT = {'mcgt': '94X_mc2017_realistic_v10', 'datagt': '94X_dataRun2_ReReco_EOY17_v2'}
process.GlobalTag = GlobalTag(process.GlobalTag, GT[envvar], '')

##################
### JEC source ###
##################
# this is if we need to override the jec in global tag
#sqfile = 'DevTools/Ntuplizer/data/{0}.db'.format('Summer16_23Sep2016V3_MC' if options.isMC else 'Summer16_23Sep2016AllV3_DATA')
#if options.crab: sqfile = 'src/{0}'.format(sqfile) # uncomment to submit to crab
#tag = 'JetCorrectorParametersCollection_Summer16_23Sep2016AllV3_DATA_AK4PFchs'
#if options.isMC: tag = 'JetCorrectorParametersCollection_Summer16_23Sep2016V3_MC_AK4PFchs' # MoriondMC
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
#process.jec = cms.ESSource("PoolDBESSource",
#    DBParameters = cms.PSet(
#        messageLevel = cms.untracked.int32(0)
#    ),
#    timetype = cms.string('runnumber'),
#    toGet = cms.VPSet(
#        cms.PSet(
#            record = cms.string('JetCorrectionsRecord'),
#            tag    = cms.string(tag),
#            label  = cms.untracked.string('AK4PFchs')
#        ),
#    ), 
#    connect = cms.string('sqlite:{0}'.format(sqfile)),
#)
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

#############################
### Setup rest of running ###
#############################
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    skipEvents = cms.untracked.uint32(options.skipEvents),
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(options.outputFile),
)

process.schedule = cms.Schedule()

###########################
### Profiling utilities ###
###########################

#process.ProfilerService = cms.Service (
#      "ProfilerService",
#       firstEvent = cms.untracked.int32(2),
#       lastEvent = cms.untracked.int32(500),
#       paths = cms.untracked.vstring('schedule') 
#)

#process.SimpleMemoryCheck = cms.Service(
#    "SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)

### To use IgProf's neat memory profiling tools, uncomment the following 
### lines then run this cfg with igprof like so:
###      $ igprof -d -mp -z -o igprof.mp.gz cmsRun ... 
### this will create a memory profile every 250 events so you can track use
### Turn the profile into text with
###      $ igprof-analyse -d -v -g -r MEM_LIVE igprof.yourOutputFile.gz > igreport_live.res
### To do a performance profile instead of a memory profile, change -mp to -pp
### in the first command and remove  -r MEM_LIVE from the second
### For interpretation of the output, see http://igprof.org/text-output-format.html

#from IgTools.IgProf.IgProfTrigger import igprof
#process.load("IgTools.IgProf.IgProfTrigger")
#process.igprofPath = cms.Path(process.igprof)
#process.igprof.reportEventInterval     = cms.untracked.int32(250)
#process.igprof.reportToFileAtBeginJob  = cms.untracked.string("|gzip -c>igprof.begin-job.gz")
#process.igprof.reportToFileAtEvent = cms.untracked.string("|gzip -c>igprof.%I.%E.%L.%R.event.gz")
#process.schedule.append(process.igprofPath)

# first create collections to analyze
collections = {
    'genParticles' : 'prunedGenParticles',
    'electrons'    : 'slimmedElectrons',
    'muons'        : 'slimmedMuons',
    'taus'         : 'slimmedTaus',
    'photons'      : 'slimmedPhotons',
    'jets'         : 'slimmedJets',
    'pfmet'        : 'slimmedMETs',
    'rho'          : 'fixedGridRhoFastjetAll',
    'vertices'     : 'offlineSlimmedPrimaryVertices',
    'packed'       : 'packedPFCandidates',
}

# the selections for each object (to be included in ntuple)
# will always be the last thing done to the collection, so can use embedded things from previous steps
selections = {
    'electrons'   : 'pt>7 && abs(eta)<2.5',
    'muons'       : 'pt>4 && abs(eta)<2.4',
    'taus'        : 'pt>17 && abs(eta)<2.3',
    'photons'     : 'pt>10 && abs(eta)<3.0',
    'jets'        : 'pt>15 && abs(eta)<2.5',
}
if options.isMC:
    selections['genParticles'] = 'pt>4'

# requirements to store events
minCounts = {
    #'electrons' : 0,
    #'muons'     : 0,
    #'taus'      : 0,
    #'photons'   : 0,
    'jets'      : 1,
}

# maximum candidates to store
# selects the first n in the collection
# patobjects are pt ordered
# vertices has pv first
maxCounts = {
    'vertices': 1,
}

# selection for cleaning (objects should match final selection)
# just do at analysis level
cleaning = {
}


# now do any customization/cleaning
print 'Customizing jets'
from DevTools.Ntuplizer.customizeJets import customizeJets
collections = customizeJets(
    process,
    collections,
    isMC=bool(options.isMC),
    isREMINIAOD=bool(options.isREMINIAOD),
    reHLT=bool(options.reHLT),
)


# select desired objects
print 'Selecting objects'
from DevTools.Ntuplizer.objectTools import objectSelector, objectCleaner
for coll in selections:
    collections[coll] = objectSelector(process,coll,collections[coll],selections[coll])

# add the analyzer
process.load("DevTools.Ntuplizer.JetTree_cfi")

process.miniTree.isData = not options.isMC
process.miniTree.filterResults = cms.InputTag('TriggerResults', '', 'PAT') if options.isMC else cms.InputTag('TriggerResults', '', 'RECO')
#process.miniTree.filterResults = cms.InputTag('TriggerResults', '', 'PAT')
process.miniTree.vertexCollections.vertices.collection = collections['vertices']
if options.isMC:
    from DevTools.Ntuplizer.branchTemplates import genParticleBranches 
    process.miniTree.collections.genParticles = cms.PSet(
        collection = cms.InputTag(collections['genParticles']),
        branches = genParticleBranches,
    )
process.miniTree.collections.jets.collection = collections['jets']
process.miniTree.rho = collections['rho']
for coll, count in minCounts.iteritems():
    if  process.miniTree.vertexCollections.hasParameter(coll):
        getattr(process.miniTree.vertexCollections,coll).minCount = cms.int32(count)
    if process.miniTree.collections.hasParameter(coll):
        getattr(process.miniTree.collections,coll).minCount = cms.int32(count)
    else:
        print 'Unrecognized collection {0}'.format(coll)
for coll, count in maxCounts.iteritems():
    if process.miniTree.vertexCollections.hasParameter(coll):
        getattr(process.miniTree.vertexCollections,coll).maxCount = cms.int32(count)
    elif process.miniTree.collections.hasParameter(coll):
        getattr(process.miniTree.collections,coll).maxCount = cms.int32(count)
    else:
        print 'Unrecognized collection {0}'.format(coll)

process.miniTreePath = cms.Path()
process.miniTreePath += process.miniTree
process.schedule.append(process.miniTreePath)
