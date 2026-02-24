
createAnnotationsFromPixelClassifier("SpleenRegions", 1000.0, 100.0, "DELETE_EXISTING")

createAnnotationsFromPixelClassifier("SmallVessel", 10.0, 0.0, "SPLIT")

addShapeMeasurements("AREA", "LENGTH", "CIRCULARITY", "SOLIDITY", "MAX_DIAMETER", "MIN_DIAMETER", "NUCLEUS_CELL_RATIO")

qupath.ext.instanseg.core.InstanSeg.builder()
    .modelPath("/home/smith6jt/QuPath/v0.6/downloaded/fluorescence_nuclei_and_cells-0.1.1")
    .device("gpu0")
    .inputChannels([ColorTransforms.createChannelExtractor("DAPI"), ColorTransforms.createChannelExtractor("CD45RO"), ColorTransforms.createChannelExtractor("Ki67"), ColorTransforms.createChannelExtractor("FOXP3"), ColorTransforms.createChannelExtractor("CD38"), ColorTransforms.createChannelExtractor("CD20"), ColorTransforms.createChannelExtractor("CD4"), ColorTransforms.createChannelExtractor("CD44"), ColorTransforms.createChannelExtractor("CD31"), ColorTransforms.createChannelExtractor("CD11c"), ColorTransforms.createChannelExtractor("CD34"), ColorTransforms.createChannelExtractor("CD107a"), ColorTransforms.createChannelExtractor("PDL1"), ColorTransforms.createChannelExtractor("CD163"), ColorTransforms.createChannelExtractor("HLA-DR"), ColorTransforms.createChannelExtractor("CD68"), ColorTransforms.createChannelExtractor("CD8"), ColorTransforms.createChannelExtractor("CD21"), ColorTransforms.createChannelExtractor("CD66"), ColorTransforms.createChannelExtractor("CD141"), ColorTransforms.createChannelExtractor("CD57"), ColorTransforms.createChannelExtractor("CD3e"), ColorTransforms.createChannelExtractor("HLA-A"), ColorTransforms.createChannelExtractor("PD-1"), ColorTransforms.createChannelExtractor("CD45"), ColorTransforms.createChannelExtractor("Podoplanin")])
    .outputChannels()
    .tileDims(512)
    .interTilePadding(32)
    .nThreads(24)
    .makeMeasurements(true)
    .randomColors(false)
    .outputType("Default")
    .build()
    .detectObjects()
detectionCentroidDistances(false)
detectionToAnnotationDistancesSigned(false)
