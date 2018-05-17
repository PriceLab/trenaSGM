library(RUnit)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()

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

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      chrom=chromosome,
                      start=start,
                      end=end,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   modelBuilder <- ModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)

   checkEquals(is(modelBuilder), "ModelBuilder")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
