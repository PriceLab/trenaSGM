library(RUnit)
library(trenaSGM)
library(MotifDb)
library(motifStack)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx")){
   filename <- system.file(package="trenaSGM", "extdata", "mayo.tcx.RData")
   if(!file.exists(filename))
      filename <- file.path("../extdata", "mayo.tcx.RData")
   stopifnot(file.exists(filename))
   load(filename)
   }


genome <- "hg38"
targetGene <- "TREM2"

if(!exists("sgm"))
   sgm <- trenaSGM(genome, targetGene)

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

   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186

      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   build.spec.2 <- build.spec
   build.spec.2$title <- "fp.2000up.200down.02"
   build.spec.2$tfPrefilterCorrelation=0.2
   build.spec.2$orderModelByColumn="pcaMax"
   build.spec.2$tfMapping <- c("MotifDb", "TFClass")

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

      # now summarize the two models.
   tbl.summary <- summarizeModels(sgm, orderBy="rfScore", maxTFpredictors=6)
   checkEquals(ncol(tbl.summary), 4)
   checkEquals(colnames(tbl.summary), c("one", "two", "rank.sum", "observed"))
   checkTrue(nrow(tbl.summary) > 5)   # expect 8, but stochasticity may arise

} # test_trem2_fpdb
#------------------------------------------------------------------------------------------------------------------------
test_summarizeModels <- function()
{
   printf("--- test_summarizeModels")

   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186

      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream

   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   spec.1 <- list(title="fp.2000up.200down.04",
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   spec.2 <- spec.1
   spec.2$title <- "fp.2000up.200down.02"
   spec.2$tfPrefilterCorrelation=0.2

   spec.3 <- list(title="trem2.rmm.2000up.200down",
                  type="regions.motifMatching",
                  tss=tss,
                  regions=tbl.regions,
                  matrix=mtx,
                  pfms=query(MotifDb, "sapiens", "jaspar2018"),
                  matchThreshold=90,
                  motifDiscovery="matchPWM",
                  tfMapping="MotifDB",
                  tfPrefilterCorrelation=0.4,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   candidate.tfs <- c("IRF5", "IKZF1", "LYL1", "SPI1", "CEBPA", "TFEC",
                      "BHLHE41", "IRF8", "TAL1","ELK3", "POU2F2", "MAFB",
                      "ZBTB18", "bogus")
   spec.4 <- list(title="trem2.noDNA.allTFs",
                  type="noDNA.tfsSupplied",
                  matrix=mtx,
                  tfs=candidate.tfs,
                  tfPrefilterCorrelation=0.4,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   strategies <- list(one=spec.1, two=spec.2, three=spec.3, four=spec.4)

   x <- calculate(sgm, strategies)

   tbl.summary <- summarizeModels(sgm, orderBy="rfScore", maxTFpredictors=8)
     # for a tf to be in the summary, it must appear in at least one model
     # thus there should be a non-NA entry in the first four columns of each row
     # (the fifth and sixth columns are counts, and will always be non-NA: so just
     #  check the first four columns)
   checkTrue(all(as.integer(apply(tbl.summary, 1, function(row) length(which(!is.na(row[1:4]))))) > 0))
      # two tfs not observed in the 0.4 abs(correlation) prefilter:
   checkEquals(ncol(tbl.summary), 6)
   checkEquals(colnames(tbl.summary), c("one", "two", "three", "four", "rank.sum", "observed"))
   checkTrue(nrow(tbl.summary) > 10)   # 17 on (23 may 2018)

} # test_summarizeModels
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()

