library(RUnit)
library(trenaSGM)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_trem2_fpdb()
   test_summarizeModels()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   sgm <- trenaSGM("hg38", "TREM2")
   checkEquals(is(sgm), "trenaSGM")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_trem2_fpdb <- function()
{
   printf("--- test_trem2_fpdb")

   genome <- "hg38"
   targetGene <- "TREM2"

   sgm <- trenaSGM(genome, targetGene)

   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186

      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream

   build.spec <- list(title="fp.2000up.200down.04",
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
                      orderModelBy="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   build.spec.2 <- build.spec
   build.spec.2$title <- "fp.2000up.200down.02"
   build.spec.2$tfPrefilterCorrelation=0.2
   build.spec.2$orderModelBy="pcaMax"

   strategies <- list(one=build.spec, two=build.spec.2)

   models <- calculate(sgm, strategies)
   checkTrue(is.list(models))
   checkEquals(names(models), c("one", "two"))
   checkEquals(names(models$one), c("model", "regulatoryRegions"))
   checkEquals(names(models$two), c("model", "regulatoryRegions"))

      # the two different tfPrefilterCorrelation values produce models with very different row counts
   checkTrue(nrow(models$one$model) < 8)
   checkTrue(nrow(models$two$model) > 20)
   checkTrue(nrow(models$one$regulatoryRegions) < 35)
   checkTrue(nrow(models$two$regulatoryRegions) > 250)

      # check the sort order on the requested orderByColumns

   checkEquals(models$one$model$rfScore, sort(models$one$model$rfScore, decreasing=TRUE))
   checkEquals(models$two$model$pcaMax, sort(models$two$model$pcaMax, decreasing=TRUE))

} # test_trem2_fpdb
#------------------------------------------------------------------------------------------------------------------------
test_summarizeModels <- function()
{
   printf("--- test_summarizeModels")

   genome <- "hg38"
   targetGene <- "TREM2"

   sgm <- trenaSGM(genome, targetGene)

   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186

      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream

   build.spec <- list(title="fp.2000up.200down.04",
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
                      orderModelBy="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   build.spec.2 <- build.spec
   build.spec.2$title <- "fp.2000up.200down.02"
   build.spec.2$tfPrefilterCorrelation=0.2

   strategies <- list(one=build.spec, two=build.spec.2)
   x <- calculate(sgm, strategies)

   tbl.summary <- summarizeModels(sgm, orderBy="rfScore", maxTFpredictors=8)
   checkEquals(as.numeric(apply(tbl.summary, 1, function(row) sum(row[1:2], na.rm=TRUE))),
               tbl.summary[, "rank.sum"])
      # two tfs not observed in the 0.4 abs(correlation) prefilter:
   checkEquals(length(which(is.na(tbl.summary[, 1]))), 2)

} # test_summarizeModels
#------------------------------------------------------------------------------------------------------------------------
