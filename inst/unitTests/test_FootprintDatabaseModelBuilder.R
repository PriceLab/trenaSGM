library(RUnit)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_build.small.fimo.motifDB.mapping.cor04()
   test_build.small.fimo.motifDB.mapping.cor.02()
   test_build.10kb.fimo.motifDB.mapping.cor04()
   test_build.10kb.fimo.tfclass.mapping.cor04()

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
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      #chrom=chromosome,
                      #start=start,
                      #end=end,
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

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)

   checkTrue("FootprintDatabaseModelBuilder" %in% is(fpBuilder))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_build.small.fimo.motifDB.mapping.cor04 <- function()
{
   printf("--- test_build.small.fimo.motifDB.mapping.cor04")

   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186

      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.2000up.200down",
                      type="database.footprints",
                      #chrom=chromosome,
                      #start=start,
                      #end=end,
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

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
   x <- build(fpBuilder)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regions <- x$regulatoryRegions
   tbl.model <- x$model
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   #browser()
   #xyz <- "test_FPDB, line 93"
   expected.tfs <- c("IKZF1", "CEBPA", "IRF8", "TAL1", "NR6A1", "IRF2")

   checkEquals(tbl.model$gene, expected.tfs)
   checkTrue(all(expected.tfs %in% tbl.regions$geneSymbol))
   checkTrue(all(tbl.regions$chrom == chromosome))
   checkTrue(all(tbl.regions$fp_start >= start))
   checkTrue(all(tbl.regions$fp_end <= end))

      #---------------------------------------------------------------------------
      # now split up the small 2200 bp promoter region into duplicate
      # chrom/start/stop rows.  the footprints and models return should be
      # identical with what was obtained above
      #---------------------------------------------------------------------------

   tbl.regions <- data.frame(chrom=rep(chromosome, 2),
                             start=rep(start, 2),
                             end=rep(end, 2),
                             stringsAsFactors=FALSE)
   build.spec$region <- tbl.regions
   fpBuilder.2 <- FootprintDatabaseModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
   x.2 <- build(fpBuilder.2)

   checkEquals(names(x.2), c("model", "regulatoryRegions"))
   tbl.regions.2 <- x.2$regulatoryRegions
   tbl.model.2 <- x.2$model
   #tbl.model.2 <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   expected.tfs <- c("IKZF1", "CEBPA", "IRF8", "TAL1", "NR6A1", "IRF2")

   #browser()
   checkEquals(tbl.model.2$gene, expected.tfs)
   checkTrue(all(expected.tfs %in% tbl.regions.2$geneSymbol))
   checkTrue(all(tbl.regions.2$chrom == chromosome))
   checkTrue(all(tbl.regions.2$fp_start >= start))
   checkTrue(all(tbl.regions.2$fp_end <= end))

   #browser()
   #xyz <- "test_FPDB, line 93"


      #---------------------------------------------------------------------------
      # now split up the small 2200 bp promoter region into two partly overlapping
      # chrom/start/stop rows.  the footprints and models return should be
      # identical with what was obtained above
      #---------------------------------------------------------------------------

   tbl.regions <- data.frame(chrom=rep(chromosome, 2),
                             start=c(41162986, 41163986),
                             end=c(41163886, 41165186),
                             stringsAsFactors=FALSE)
   build.spec$region <- tbl.regions
   fpBuilder.3 <- FootprintDatabaseModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
   x.3 <- build(fpBuilder.3)

   checkEquals(names(x.3), c("model", "regulatoryRegions"))
   tbl.regions.3 <- x.3$regulatoryRegions
   tbl.model.3 <- x.3$model
   expected.tfs <- c("IKZF1", "CEBPA", "IRF8", "TAL1", "NR6A1", "IRF2")

   checkEquals(tbl.model.3$gene, expected.tfs)
   checkTrue(all(expected.tfs %in% tbl.regions.3$geneSymbol))
   checkTrue(all(tbl.regions.3$chrom == chromosome))
   checkTrue(all(tbl.regions.3$fp_start >= start))
   checkTrue(all(tbl.regions.3$fp_end <= end))

} # test_build.small.fimo.motifDB.mapping.cor04
#------------------------------------------------------------------------------------------------------------------------
test_build.small.fimo.motifDB.mapping.cor.02 <- function()
{
   printf("--- test_build.small.fimo.motifDB.mapping.cor.02")

   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.2000up.200down",
                      type="database.footprints",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.2,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder("hg38", "TREM2", build.spec, quiet=TRUE)
   x <- build(fpBuilder)
   tbl.regions <- x$regulatoryRegions
   tbl.model <- x$model
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   checkTrue(nrow(tbl.model) > 20)
   top.tfs <- subset(tbl.model, rfScore >= 4)$gene
   checkTrue(all(top.tfs %in% tbl.regions$geneSymbol))

} # test_build.small.fimo.motifDB.mapping.cor02
#------------------------------------------------------------------------------------------------------------------------
test_build.10kb.fimo.motifDB.mapping.cor04 <- function()
{
   printf("--- test_build.10kb.fimo.motifDB.mapping.cor04")

   chromosome <- "chr6"
   upstream <- 5000
   downstream <- 5000
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.5kbup.5kbdown",
                      type="database.footprints",
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

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FootprintDatabaseModelBuilder("hg38", "TREM2", build.spec, quiet=TRUE)
   x <- build(fpBuilder)
   tbl.regions <- x$regulatoryRegions
   tbl.model <- x$model
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   checkTrue(nrow(tbl.model) > 10)
   top.tfs <- subset(tbl.model, rfScore >= 4)$gene
   checkTrue(all(top.tfs %in% tbl.regions$geneSymbol))

} # test_build.10kb.fimo.motifDB.mapping.cor04
#------------------------------------------------------------------------------------------------------------------------
test_build.10kb.fimo.tfclass.mapping.cor04 <- function()
{
   printf("--- test_build.10kb.fimo.tfclass.mapping.cor04")

   chromosome <- "chr6"
   upstream <- 5000
   downstream <- 5000
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.5kbup.5kbdown",
                      type="database.footprints",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfMapping="TFClass",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FootprintDatabaseModelBuilder("hg38", "TREM2", build.spec, quiet=TRUE)
   x <- build(fpBuilder)
   tbl.regions <- x$regulatoryRegions
   tbl.model <- x$model
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   checkTrue(nrow(tbl.model) > 10)
   top.10.tfs <- tbl.model$gene[1:10]   # already sorted by rfScore
   expected.tfs <- c("SPI1", "TFEC", "CEBPA", "BACH1", "TAL1", "FLI1", "ELK3", "KLF3", "FOXP1", "MESP1")
      # allow for some randomness in the result
   checkTrue(length(intersect(top.10.tfs, expected.tfs)) >= 8)
   checkTrue(all(top.10.tfs %in% tbl.regions$geneSymbol))

} # test_build.10kb.fimo.motifDB.mapping.cor04
#------------------------------------------------------------------------------------------------------------------------
