library(RUnit)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_build.trem2.noDNA.allTFs()

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
                      type="noDNA_tfsSupplied",
                      tfs=c("HLF", "STAT4", "SATB2", "SATB1", "TSHZ3", "TSHZ2", "FOXP2"),
                      matrix=mtx,
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- NoDnaModelBuilder(genome, targetGene, build.spec, quiet=TRUE)

   checkTrue("NoDnaModelBuilder" %in% is(builder))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.noDNA.allTFs <- function()
{
   printf("--- test_build.trem2.noDNA.allTFs")

   genome <- "hg38"
   targetGene <- "TREM2"

   candidate.tfs <-   c("HLF", "STAT4", "SATB2", "SATB1", "TSHZ3", "TSHZ2", "FOXP2",
                        "FOXP1", "LHX6", "BACH1", "SOX12", "FOXD4L1", "NFE2L2", "ZHX3",
                        "ZBTB16", "ZHX1", "TAF1", "STAT6", "POU4F1", "FOXD1", "ATF2",
                        "BCL6B", "STAT5B", "NR5A2", "FOXE3", "STAT3", "ATF7", "STAT2")


   build.spec <- list(title="trem2.noDNA.allTFs",
                      type="noDNA_tfsSupplied",
                      matrix=mtx,
                      tfs=candidate.tfs,
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   checkEquals(x$regulatoryRegions, data.frame())
   tbl.model <- x$model
   expected.tfs <- sort(c("BACH1", "FOXP1", "STAT3"))
   checkTrue(all(expected.tfs %in% tbl.model$gene))
   checkTrue(all(expected.tfs %in% candidate.tfs))
   checkTrue(all(is.na(tbl.model$bindingSites)))

} # test_build.trem2.noDNA.allTFs
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
