library(RUnit)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_constructor_failure()
   test_.replaceGeneSymbolsWithEnsemblGeneIDs()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      tfPrefilterCorrelation=0.4,
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   modelBuilder <- ModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)

   checkEquals(is(modelBuilder), "ModelBuilder")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_constructor_failure <- function()
{
   printf("--- test_constructor_failure")

   genome <- "hg38"
   targetGene <- "bogus"
   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      tfPrefilterCorrelation=0.4,
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   checkException(modelBuilder <- ModelBuilder(genome, targetGene,  build.spec, quiet=TRUE), silent=TRUE)

} # test_constructor_failure
#------------------------------------------------------------------------------------------------------------------------
test_.replaceGeneSymbolsWithEnsemblGeneIDs <- function()
{
   printf("--- test_.replaceGeneSymbolsWithEnsemblGeneIDs")

   if(!exists("tbl.regulatoryRegions"))
      load(system.file(package="trenaSGM", "extdata", "tbl.regForSymbolEnsemblTranslationTest.RData"))

      # start with s small simple subet, with 2 NA geneSybmbols
      # this initially random subset has two NAs, quite a few mouse genes

   set.seed(17); indices <- sample(1:nrow(tbl.regulatoryRegions), 30)
   tbl.test <- tbl.regulatoryRegions[indices, c("motifName", "geneSymbol")]
   tbl.fixed <- trenaSGM:::.replaceGeneSymbolsWithEnsemblGeneIDs(tbl.test)
   checkEquals(nrow(tbl.test), 30)
   checkEquals(nrow(tbl.fixed), 28)
   checkEquals(length(grep("^NA", tbl.test$motifName)), 2)

   #geneSymbol.na <- which(is.na(tbl.test$geneSymbol))
   #fixed.na <- which(is.na(tbl.fixed$geneSymbol))
   #checkEquals(geneSymbol.na, fixed.na)
   successful <- length(grep("^ENSG", tbl.fixed$geneSymbol))
   checkEquals(successful, 22)   # 2 na, 6 unmapped gene symbols, 22 ensemblIDs

      # now a larger table
   set.seed(17); indices <- sample(1:nrow(tbl.regulatoryRegions), 800)[300:350]
   tbl.test <- tbl.regulatoryRegions[indices, c("motifName", "geneSymbol")]
   tbl.fixed <- trenaSGM:::.replaceGeneSymbolsWithEnsemblGeneIDs(tbl.test)
   #tbl.fixed <- trenaSGM:::.replaceGeneSymbolsWithEnsemblGeneIDs(tbl.test)

   tbl.2 <- tbl.regulatoryRegions[, c("motifName", "geneSymbol")]
   tbl.2f <- trenaSGM:::.replaceGeneSymbolsWithEnsemblGeneIDs(tbl.2)
   checkTrue(nrow(tbl.2f) > 900)
   checkEquals(length(grep("^NA$", tbl.2f$geneSymbol)), 0)
   checkTrue(length(grep("^ENSG0", tbl.2f$geneSymbol)) > 700)

} # test_.replaceGeneSymbolsWithEnsemblGeneIDs
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
