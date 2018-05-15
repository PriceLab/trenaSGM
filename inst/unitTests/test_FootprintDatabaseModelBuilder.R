library(RUnit)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_build()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   fp.specs <- list(title="fp.2000up.200down",
                    type="database.footprints",
                    matrix=mtx,
                    chrom="chr6",
                    tss=41163186,
                    upstream=2000,
                    downstream=200,
                    db.host="khaleesi.systemsbiology.net",
                    databases=list("brain_hint_20"),
                    motifDiscovery="builtinFimo",
                    tfMapping="MotifDB",
                    solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                    tfPrefilterCorrelation=0.2)

   fpBuilder <- FootprintDatabaseModelBuilder("hg38", "TREM2", fp.specs, quiet=TRUE)

   checkEquals(is(fpBuilder), "FootprintDatabaseModelBuilder")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_build <- function()
{
   printf("--- test_build")

   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream

   fp.specs <- list(title="fp.2000up.200down",
                    type="database.footprints",
                    chrom=chromosome,
                    start=start,
                    end=end,
                    tss=tss,
                    db.host="khaleesi.systemsbiology.net",
                    databases=list("brain_hint_20"),
                    motifDiscovery="builtinFimo",
                    tfMapping="MotifDB",
                    tfPrefilterCorrelation=0.4,
                    solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder("hg38", "TREM2", fp.specs, quiet=TRUE)
   x <- build(fpBuilder, mtx)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regions <- x$regulatoryRegions
   tbl.model <- x$model
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   expected.tfs <- c("IKZF1", "IRF8", "CEBPA", "TAL1",  "IRF2")
   checkEquals(tbl.model$gene, expected.tfs)
   checkTrue(all(expected.tfs %in% tbl.regions$geneSymbol))
   checkTrue(all(tbl.regions$chrom == chromosome))
   checkTrue(all(tbl.regions$fp_start >= start))
   checkTrue(all(tbl.regions$fp_end <= end))


} # test_build
#------------------------------------------------------------------------------------------------------------------------
