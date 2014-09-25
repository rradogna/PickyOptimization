import FWCore.ParameterSet.Config as cms

process = cms.Process("pickyANA")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('PRE_STA71_V3::All')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/r/rradogna/TuneP/CMSSW_7_1_0_pre7/src/PikyOpt/reco_RAW2DIGI_RECO.root'
    )
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("testPicky.root"),
      closeFileFast = cms.untracked.bool(True)
  )
          
process.analizePicky = cms.EDAnalyzer('PickyAnalyzer',
#			       	RootfileName = cms.untracked.string('testPicky.root'),
				DataType = cms.untracked.string('SimData'),
				MuonCollectionLabel = cms.untracked.InputTag('muons','','RECO'),
)


  
#process.MessageLogger = cms.Service("MessageLogger",
#    cout = cms.untracked.PSet(
#        threshold = cms.untracked.string('INFO'),
#        noLineBreaks = cms.untracked.bool(True)
#    ),
#    destinations = cms.untracked.vstring('cout')
#)  

process.p = cms.Path(process.analizePicky)

#process.this_is_the_end = cms.EndPath(process.out)
