import FWCore.ParameterSet.Config as cms

process = cms.Process("MTDAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.Geometry.GeometryExtended2023D35_cff")

process.load("Geometry.MTDNumberingBuilder.mtdNumberingGeometry_cfi")

process.load("Geometry.MTDNumberingBuilder.mtdTopology_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdGeometry_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdParameters_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(100),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:step3.root',
    )
)


process.mtdGeometry = cms.ESProducer("MTDDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(True)
)


process.MTDAnalyzer = cms.EDAnalyzer('MTDAnalyzer',
                                     BTLIntegrationWindow = cms.double(25.), # [ns]
                                     BTLMinimumEnergy     = cms.double(2.)   # [MeV]
                                     )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('MTDAnalyzer_histo.root')
                                   )

process.p = cms.Path(process.MTDAnalyzer)
