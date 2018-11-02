library(RUnit)
library(trenaSGM)
library(motifStack)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_build.trem2.2200bp.motifDB.model()
   test_build.trem2.2200bp.TFClass.model()
   test_build.trem2.2200bp.TFClass.and.MotifDb.model()

   test_IKZF1()
   test_IRF5()
   test_BHLHE41()

   test_build.trem2.enhancers.TFClass.and.MotifDb.model()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   build.spec <- list(title="trem2.rmm.2000up.200down",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, "sapiens", "jaspar2018"),
                      motifDiscovery="matchPWM",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=TRUE)

   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)

   checkTrue("RegionsMotifMatchingModelBuilder" %in% is(builder))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.2200bp.motifDB.model <- function()
{
   printf("--- test_build.trem2.2200bp.motifDb.model")

   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   build.spec <- list(title="trem2.rmm.2000up.200down",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, "sapiens", "jaspar2018"),
                      matchThreshold=90,
                      motifDiscovery="matchPWM",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=FALSE)

   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regRegions <- x$regulatoryRegions
   tbl.model <- x$model
   expected.tfs <- sort(c("SPI1", "STAT3", "RUNX1", "NFATC2", "FLI1"))
   checkEquals(sort(tbl.model$gene), expected.tfs)
   checkTrue(all(expected.tfs %in% tbl.regRegions$geneSymbol))
   checkTrue(all(tbl.regRegions$chrom == chromosome))
   checkTrue(all(tbl.regRegions$fmotifStart >= min(tbl.regions$start)))
   checkTrue(all(tbl.regRegions$motifEnd <= max(tbl.regions$end)))

} # test_build.trem2.2200bp.motifDB.model
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.2200bp.TFClass.model <- function()
{
   printf("--- test_build.trem2.2200bp.TFClass.model")

   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   build.spec <- list(title="trem2.rmm.2000up.200down",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, "sapiens", "jaspar2018"),
                      matchThreshold=90,
                      motifDiscovery="matchPWM",
                      tfPool=allKnownTFs(),
                      tfMapping="TFClass",
                      tfPrefilterCorrelation=0.4,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=FALSE)


   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regRegions <- x$regulatoryRegions
   tbl.model <- x$model
   expected.top.tfs <- sort(c("LYL1", "TFEC", "NFATC2", "TAL1", "ELK3"))
   checkEquals(sort(tbl.model$gene[1:5]), expected.top.tfs)
   checkTrue(all(expected.top.tfs %in% tbl.regRegions$geneSymbol))
   checkTrue(all(tbl.regRegions$chrom == chromosome))
   checkTrue(all(tbl.regRegions$fmotifStart >= min(tbl.regions$start)))
   checkTrue(all(tbl.regRegions$motifEnd <= max(tbl.regions$end)))

} # test_build.trem2.2200bp.TFClass.model
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.2200bp.TFClass.and.MotifDb.model <- function()
{
   printf("--- test_build.trem2.2200bp.TFClass.and.MotifDbmodel")

   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   build.spec <- list(title="trem2.rmm.2000up.200down",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, "sapiens", "jaspar2018"),
                      matchThreshold=90,
                      motifDiscovery="matchPWM",
                      tfPool=allKnownTFs(),
                      tfMapping=c("MotifDb", "TFClass"),
                      tfPrefilterCorrelation=0.4,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=FALSE)


   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regRegions <- x$regulatoryRegions
   tbl.model <- x$model
   expected.top.tfs <- sort(c("LYL1", "SPI1", "TFEC", "NFATC2", "TAL1"))
   checkEquals(sort(tbl.model$gene[1:5]), expected.top.tfs)
   checkTrue(all(expected.top.tfs %in% tbl.regRegions$geneSymbol))
   checkTrue(all(tbl.regRegions$chrom == chromosome))
   checkTrue(all(tbl.regRegions$fmotifStart >= min(tbl.regions$start)))
   checkTrue(all(tbl.regRegions$motifEnd <= max(tbl.regions$end)))

} # test_build.trem2.2200bp.TFClass.and.MotifDb.model
#------------------------------------------------------------------------------------------------------------------------
# with my historical preference for jaspar2018, the crucial tf of TREM2, IKZF1, is missed: only HOCOMOCO
# and SwissRegulon report a motif for it.  this test demonstrates the solution, and thus the importance of
# motif sources other than MotifDb
test_IKZF1 <- function()
{
   printf("--- test_IKZF1")
   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   build.spec <- list(title="trem2.2000up.200down-find.ikzf1",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, andStrings="sapiens", orString=c("jaspar2018", "hocomoco")),
                      matchThreshold=90,
                      motifDiscovery="matchPWM",
                      tfPool=allKnownTFs(),
                      tfMapping=c("MotifDb", "TFClass"),
                      tfPrefilterCorrelation=0.4,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=FALSE)

   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regRegions <- x$regulatoryRegions
   tbl.model <- x$model
   expected.top.tfs <- sort(c("IKZF1", "LYL1", "SPI1", "TFEC", "NFATC2", "TAL1"))
   checkEquals(sort(tbl.model$gene[1:6]), expected.top.tfs)
   checkTrue(all(expected.top.tfs %in% tbl.regRegions$geneSymbol))
   checkTrue(all(tbl.regRegions$chrom == chromosome))
   checkTrue(all(tbl.regRegions$fmotifStart >= min(tbl.regions$start)))
   checkTrue(all(tbl.regRegions$motifEnd <= max(tbl.regions$end)))

} # test_IKZF1
#------------------------------------------------------------------------------------------------------------------------
# irf5 matches (about 4kb upstream of tss) only two older motifs, and then only at 82%
# show that here
test_IRF5 <- function()
{
   printf("--- test_IRF5")
   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE)

   build.spec <- list(title="trem2.2000up.200down-find.irf5",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, andStrings=c("sapiens", "irf5")),
                      matchThreshold=80,
                      motifDiscovery="matchPWM",
                      tfPool=allKnownTFs(),
                      tfMapping=c("MotifDb", "TFClass"),
                      tfPrefilterCorrelation=0.4,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=FALSE)

   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   suppressWarnings(x <- build(builder))
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regRegions <- x$regulatoryRegions
   tbl.model <- x$model
   checkEquals(tbl.model$gene[1], "IRF5")    # just one tf, IRF5
   checkTrue(all(tbl.regRegions$chrom == chromosome))
   checkTrue(all(tbl.regRegions$fmotifStart >= min(tbl.regions$start)))
   checkTrue(all(tbl.regRegions$motifEnd <= max(tbl.regions$end)))

} # test_IRF5
#------------------------------------------------------------------------------------------------------------------------
# fails to make it into a standard regions model unless pwm match threshold is dropped to 84.
# however, looking at the sequence of the 84.4% match,  chr6:41158763-41158772, it looks like
# a strong pwm/sequence match to me, perfect until base 9/10, where [AC] is actually a T
#  showGeomicRegion(igv, "chr6:41,158,763-41,158,772")
#  x <- query(MotifDb, andStrings=c("sapiens", "jaspar", "bhlhe41"))
#  motifStack(lapply(names(x), function(mName) new("pfm", x[[mName]], name=mName)))
test_BHLHE41 <- function()
{
   printf("--- test_BHLHE41")
   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE)

   build.spec <- list(title="trem2.2000up.200down-find.bhlhe41",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, andStrings=c("sapiens", "jaspar2018", "bhlhe41")),
                      matchThreshold=80,
                      motifDiscovery="matchPWM",
                      tfPool=allKnownTFs(),
                      tfMapping=c("MotifDb", "TFClass"),
                      tfPrefilterCorrelation=0.4,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=FALSE)

   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   suppressWarnings(x <- build(builder))
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regRegions <- x$regulatoryRegions
   tbl.model <- x$model
   checkEquals(tbl.model$gene[1], "BHLHE41")    # just one tf, BHLHE41
   checkTrue(all(tbl.regRegions$chrom == chromosome))
   checkTrue(all(tbl.regRegions$fmotifStart >= min(tbl.regions$start)))
   checkTrue(all(tbl.regRegions$motifEnd <= max(tbl.regions$end)))

   checkTrue(all(tbl.regRegions$motifRelativeScore > 0.81))
   checkTrue(all(tbl.regRegions$motifRelativeScore < 0.85))
   checkTrue(tbl.model$pearsonCoeff[1] > 0.6)

} # test_BHLHE41
#------------------------------------------------------------------------------------------------------------------------
# 105 seconds to run
# how closely can we recreate cory's footprint/fimo database TREM2 model?
test_build.trem2.enhancers.TFClass.and.MotifDb.model <- function()
{
   printf("--- test_build.trem2.enhancers.TFClass.and.MotifDb.model")
   load(system.file(package="trenaSGM", "extdata", "enhancers.v47.TREM2.RData"))

   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186

   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)
   tbl.regions$width <- with(tbl.regions, 1 + end - start)
   tbl.regionsCombined <- rbind(tbl.regions, tbl.enhancers)

   build.spec <- list(title="trem2.rmm.2000up.200down",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regionsCombined,
                      matrix=mtx,
                      pfms=query(MotifDb, andStrings="sapiens", orString=c("jaspar2018", "hocomoco")),
                      matchThreshold=90,
                      tfPool=allKnownTFs(),
                      motifDiscovery="matchPWM",
                      tfPool=allKnownTFs(),
                      tfMapping=c("MotifDb", "TFClass"),
                      tfPrefilterCorrelation=0.4,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=FALSE)

   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regRegions <- x$regulatoryRegions
   tbl.model <- x$model

   corys.tfs <-  c("IRF5", "IKZF1", "LYL1", "SPI1", "CEBPA", "TFEC", "BHLHE41", "IRF8", "TAL1",
                   "BACH1", "NR6A1", "FOXP1", "RUNX1", "NR1I3", "MAFB", "HHEX", "ELK3", "POU2F2",
                   "MAF", "ZBTB18")
   checkEquals(length(intersect(tbl.model$gene[1:10], corys.tfs[1:10])), 8)
   checkEquals(setdiff(corys.tfs[1:10], tbl.model$gene[1:10]), c("IRF5", "BHLHE41"))

} # test_build.trem2.enhancers.TFClass.and.MotifDb.model
#------------------------------------------------------------------------------------------------------------------------
test_build.bogusTargetGene <- function()
{
   printf("--- test_build.bogusTargetGene")

   genome <- "hg38"
   targetGene <- "bogus"
   chromosome <- "chr6"
   tss <- 41163186

   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)
   tbl.regions$width <- with(tbl.regions, 1 + end - start)

   build.spec <- list(title="trem2.rmm.2000up.200down",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, andStrings="sapiens", orString=c("jaspar2018", "hocomoco")),
                      matchThreshold=90,
                      tfPool=allKnownTFs(),
                      motifDiscovery="matchPWM",
                      tfPool=allKnownTFs(),
                      tfMapping=c("MotifDb", "TFClass"),
                      tfPrefilterCorrelation=0.4,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=FALSE)

   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   checkEquals(nrow(x$model), 0)
   checkEquals(nrow(x$regulatoryRegions), 0)

} # test_bogusTargetGene
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
