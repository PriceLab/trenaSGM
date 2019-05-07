library(RUnit)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.new.RData"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_build.trem2.noDNA.13.known.TFs()
   test_build.trem2.noDNA.all.known.TFs()
   test_build.trem2.noDNA.bogus.targetGene()

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
                      type="noDNA.tfsSupplied",
                      tfPool=c("HLF", "STAT4", "SATB2", "SATB1", "TSHZ3", "TSHZ2", "FOXP2"),
                      matrix=mtx,
                      tfPrefilterCorrelation=0.4,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=TRUE)

   builder <- NoDnaModelBuilder(genome, targetGene, build.spec, quiet=TRUE)

   checkTrue("NoDnaModelBuilder" %in% is(builder))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.noDNA.13.known.TFs <- function()
{
   printf("--- test_build.trem2.noDNA.13.known.TFs")

   genome <- "hg38"
   targetGene <- "TREM2"

   candidate.tfs <- c("IRF5", "IKZF1", "LYL1", "SPI1", "CEBPA", "TFEC",
                      "BHLHE41", "IRF8", "TAL1","ELK3", "POU2F2", "MAFB",
                      "ZBTB18", "bogus")

   build.spec <- list(title="trem2.noDNA.allTFs",
                      type="noDNA.tfsSupplied",
                      matrix=mtx,
                      tfPool=allKnownTFs(),
                      tfs=candidate.tfs,
                      tfPrefilterCorrelation=0.2,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=TRUE)

   builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   checkEquals(x$regulatoryRegions, data.frame())
   tbl.model <- x$model
      # with a relaxed tfPrefilterCorrelation, and the hand-picked TFs listed above
      # all but "bogus" make the cut
   checkEquals(setdiff(candidate.tfs, tbl.model$gene), "bogus")
      # noDNA implies no bindingSites
   checkTrue(all(is.na(tbl.model$bindingSites)))

      # all pearsonCoeff above prefilter threshold?
   checkTrue(all(abs(tbl.model$pearsonCoeff) > 0.2))

      #--------------------------------------------------------------------------------
      # run again with a stricter tfPrefilterCorrelation
      #--------------------------------------------------------------------------------
   build.spec$tfPrefilterCorrelation <- 0.7
   builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   tbl.model <- x$model
   checkTrue(all(abs(tbl.model$pearsonCoeff) > 0.7))
   checkTrue(all(tbl.model$gene %in% candidate.tfs))

} # test_build.trem2.noDNA.13.known.TFS
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.noDNA.all.known.TFs <- function()
{
   printf("--- test_build.trem2.noDNA.all.known.TFs")

   genome <- "hg38"
   targetGene <- "TREM2"

   candidate.tfs <- allKnownTFs()

   build.spec <- list(title="trem2.noDNA.allTFs",
                      type="noDNA.tfsSupplied",
                      matrix=mtx,
                      tfPool=allKnownTFs(),
                      tfs=candidate.tfs,
                      tfPrefilterCorrelation=0.7,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="pearsonCoeff",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=TRUE)

   builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)
   tbl.model <- x$model
   checkEquals(dim(x$regulatoryRegions), c(0,0))
   checkTrue(all(tbl.model$peasonCoeff > 0.7))
      # the order
   checkEquals(tbl.model$gene, c("PLEK", "IRF5", "IKZF1", "LYL1", "SPI1", "TFEC"))

} # test_build.trem2.noDNA.all.known.TFs
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.noDNA.bogus.targetGene <- function()
{
   printf("--- test_build.trem2.noDNA.bogus.targetGene")

   genome <- "hg38"
   targetGene <- "bogusGene"

   candidate.tfs <- allKnownTFs()

   build.spec <- list(title="trem2.noDNA.allTFs",
                      type="noDNA.tfsSupplied",
                      matrix=mtx,
                      tfPool=allKnownTFs(),
                      tfs=candidate.tfs,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      tfPrefilterCorrelation=0.7,
                      orderModelByColumn="pearsonCoeff",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                      quiet=TRUE)

   checkException(builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE), silent=TRUE)

} # test_build.trem2.noDNA.bogus.targetGene
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
