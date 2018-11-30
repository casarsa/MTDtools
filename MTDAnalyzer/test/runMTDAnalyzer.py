import FWCore.ParameterSet.Config as cms

process = cms.Process("MTDAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(100),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:step3.root',
        #'file:/eos/infnts/cms/store/user/casarsa/MTD/SingleMuPt10_pythia8_2023D24_noPU/step3.root'
    )
)

process.MTDAnalyzer = cms.EDAnalyzer('MTDAnalyzer',
                                     BTLIntegrationWindow = cms.double(25.), # [ns]
                                     BTLMinimumEnergy     = cms.double(2.)   # [MeV]
                                     )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('MTDAnalyzer_histo.root')
                                   )

process.p = cms.Path(process.MTDAnalyzer)
