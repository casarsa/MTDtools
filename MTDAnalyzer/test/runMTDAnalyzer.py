import FWCore.ParameterSet.Config as cms

process = cms.Process("MTDAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(100),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:step3_0.root'
    )
)

process.MTDAnalyzer = cms.EDAnalyzer('MTDAnalyzer',
                                     BTLMinimumEnergy = cms.double(1.)  # [MeV]
                                     )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('MTDAnalyzer_histo.root')
                                   )

process.p = cms.Path(process.MTDAnalyzer)
