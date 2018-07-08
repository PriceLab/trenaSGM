library(RUnit)
library(trenaSGM)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
Sys.setlocale("LC_ALL", "C")
if(!exists("mtx")){
    load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
    all.genes <- rownames(mtx)
    suppressMessages(tbl.map <- select(org.Hs.eg.db, keys=all.genes, keytype="SYMBOL", columns="ENSEMBL"))
    dups <- which(duplicated(tbl.map$SYMBOL))
    if(length(dups) > 0)
       tbl.map <- tbl.map[-dups,]
    checkTrue(all(tbl.map$SYMBOL %in% rownames(mtx)))
    mtx.ensembl <- mtx[tbl.map$SYMBOL,]
    indices <- match(tbl.map$SYMBOL, rownames(mtx.ensembl))
    rownames(mtx.ensembl) <- tbl.map$ENSEMBL[indices]
    deleters <- which(is.na(rownames(mtx.ensembl)))
    if(length(deleters) > 0)
        mtx.ensembl <- mtx.ensembl[-deleters,]
    }

if(!exists("tbl.enhancers"))
  load(system.file(package="trenaSGM", "extdata", "enhancers.TREM2.RData"))

if(!exists("tbl.trena")){
   printf("loading cory's trem2 model for comparison")
   load(system.file(package="trenaSGM", "extdata", "ENSG00000095970.RData"))
   ensembl.ids <- tbl.trena$gene
   suppressWarnings(tbl.map <-  select(org.Hs.eg.db, keys=ensembl.ids, keytype="ENSEMBL", columns=c("SYMBOL", "ENSEMBL")))
   tbl.map <- tbl.map[-which(duplicated(tbl.map$ENSEMBL)),]   # two dups
   rownames(tbl.trena) <- tbl.map$SYMBOL
   }

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_build.small.fimo.motifDB.mapping.cor04()
   test_build.small.fimo.motifDB.mapping.cor.02()
   test_build.10kb.fimo.motifDB.mapping.cor04()
   test_build.10kb.fimo.tfclass.mapping.cor04()
   test_reproduceCorysTrem2model()
   test_tfPoolOption()
   test_modelExpressionDataWithEnsgIDs()
   test_noGenesAboveExpressionCorrelationThreshold()
   test_stagedExecution()
   test_parallel_stagedExecution.footprints()

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
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FastFootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)

   checkTrue("FastFootprintDatabaseModelBuilder" %in% is(fpBuilder))

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
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FastFootprintDatabaseModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
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
   fpBuilder.2 <- FastFootprintDatabaseModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
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
   fpBuilder.3 <- FastFootprintDatabaseModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
   x.3 <- build(fpBuilder.3)

   checkEquals(names(x.3), c("model", "regulatoryRegions"))
   tbl.regions.3 <- x.3$regulatoryRegions
   tbl.model.3 <- x.3$model
   expected.tfs <- c("IKZF1", "CEBPA", "IRF8", "TAL1", "NR6A1", "IRF2")

   checkTrue(length(intersect(tbl.model.3$gene, expected.tfs)) > 3)  # allow for some stochasticity
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
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.2,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FastFootprintDatabaseModelBuilder("hg38", "TREM2", build.spec, quiet=TRUE)
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
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FastFootprintDatabaseModelBuilder("hg38", "TREM2", build.spec, quiet=TRUE)
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
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="TFClass",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   fpBuilder <- FastFootprintDatabaseModelBuilder("hg38", "TREM2", build.spec, quiet=TRUE)
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
test_reproduceCorysTrem2model <- function()
{
   printf("--- test_reproduceCorysTrem2model")
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand

   build.spec <- list(title="fp.enhancers",
                      type="footprint.database",
                      regions=tbl.enhancers,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=c("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping=c("TFClass", "MotifDb"),
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FastFootprintDatabaseModelBuilder("hg38", "TREM2", build.spec, quiet=FALSE)
   x <- build(fpBuilder)
   tbl.regions <- x$regulatoryRegions
   tbl.model <- x$model
   checkTrue(nrow(tbl.model) > 20)
   top.tfs <- subset(tbl.model, rfScore >= 4)$gene
   checkTrue(all(top.tfs %in% tbl.regions$geneSymbol))
   checkTrue(all(c("IKZF1", "IRF5", "LYL1", "SPI1") %in% top.tfs))

      # now compare to cory's previous model, which has been designated the reference
   top.tfs.bothModels <- intersect(rownames(subset(tbl.trena, pearsonCoeff >= 0.4))[1:10], tbl.model$gene[1:10])
   checkTrue(length(top.tfs.bothModels) >= 8)  # usually exactly 9, with only BHLHE41 missing

      # the missing tf is BHLHE41, which is a poor-ish fit (~80%) to MA0636.1
      # and a better fit for the more relaxed MA0692.1
      # to see this:
      #    x <- query(MotifDb, "jaspar2018", c("MA0636.1", "MA0692.1"))
      #    motifStack(lapply(names(x), function(mName) new("pfm", x[[mName]], name=mName)))
      # we had previously mapped the permissive MA0692.1 to BHLHE41, but
      # this (mysteriously, but appropriately) is no longer the case:
      #
      # motifToGene(MotifDb, "MA0692.1", "MotifDb")
      #      motif geneSymbol pubmedID organism  source
      # 1 MA0692.1       TFEB 24194598 Hsapiens MotifDb
      # motifToGene(MotifDb, "MA0692.1", "TFClass")
      #       motif geneSymbol pubmedID  source organism
      # 1  MA0692.1       TFEB 23180794 TFClass Hsapiens
      # 2  MA0692.1       MXI1 23180794 TFClass Hsapiens
      # 3  MA0692.1        MNT 23180794 TFClass Hsapiens
      # 4  MA0692.1        MLX 23180794 TFClass Hsapiens
      # 5  MA0692.1        MAX 23180794 TFClass Hsapiens
      # 6  MA0692.1       USF1 23180794 TFClass Hsapiens
      # 7  MA0692.1       USF2 23180794 TFClass Hsapiens
      # 8  MA0692.1       TFE3 23180794 TFClass Hsapiens
      # 9  MA0692.1       TFEC 23180794 TFClass Hsapiens
      # 10 MA0692.1       MITF 23180794 TFClass Hsapiens
      # 11 MA0692.1       MYCN 23180794 TFClass Hsapiens
      # and the single consistent mappoing for MA0636.1
      #   motifToGene(MotifDb, "MA0636.1", "MotifDb")
      #       motif geneSymbol pubmedID organism  source
      #  1 MA0636.1    BHLHE41 24194598 Hsapiens MotifDb
      #  motifToGene(MotifDb, "MA0636.1", "TFClass")
      #       motif geneSymbol pubmedID  source organism
      #  1 MA0636.1    BHLHE41 23180794 TFClass Hsapiens

   checkEquals(setdiff(rownames(tbl.trena)[1:10], tbl.model$gene[1:10]), "BHLHE41")

      # now take a deeper look.   agreement is still pretty good

   checkEquals(length(intersect(tbl.model$gene, rownames(subset(tbl.trena, pearsonCoeff >= 0.4)))), 19)

} # test_reproduceCorysTrem2model
#------------------------------------------------------------------------------------------------------------------------
viz.corysTrem2Model <- function()
{
   library(igvR)
   library(motifStack)
   igv <- igvR()
   setGenome(igv, "hg38")
   if(!exists("tbl.enhancers"))
      print(load(system.file(package="trenaSGM", "extdata", "enhancers.TREM2.RData")))

   bigRegion <- with(tbl.enhancers, sprintf("%s:%d-%d", unique(chrom), min(start)-2000, max(end)+2000))
   showGenomicRegion(igv, bigRegion)
   track <- DataFrameAnnotationTrack("GeneHancer 4.6", tbl.enhancers, color="black")
   displayTrack(igv, track)

        #---- v4.7 enhancers
    load("../extdata/enhancers.v47.TREM2.RData")
    trackv47 <- DataFrameAnnotationTrack("GeneHancer 4.7", tbl.enhancers, color="purple")
    displayTrack(igv, trackv47)

     #--- the LYL1 mystery: 45 footprints which trenaSGM have not yet found
   tbl.lyl1 <- subset(tbl.footprints, tf == "LYL1")
   track.lyl1 <- DataFrameAnnotationTrack("LYL1", tbl.lyl1[, c("chrom", "fp_start", "fp_end", "motifName")], color="darkred")
   displayTrack(igv, track.lyl1)
   long.motif.names <- unique(tbl.lyl1$motifName)   # just 4
   motifStack(lapply(long.motif.names, function(mName) new("pfm", MotifDb[[mName]], name=mName)))
   short.motif.names <- mcols(MotifDb[long.motif.names])$providerId   # [1] "MA0091.1" "MA0140.2" "MA0092.1" "MA0048.2"


   tbl.ikzf1 <- subset(tbl.footprints, tf=="IKZF1")
   dups <- which(duplicated(tbl.ikzf1$loc))
   if(length(dups) > 0)
      tbl.ikzf1 <- tbl.ikzf1[-dups,]
   track <- DataFrameAnnotationTrack("cory IKZF1", tbl.ikzf1[, c("chrom", "start", "endpos")], color="green")
   displayTrack(igv, track)
   dim(subset(tbl.footprints, tf=="IKZF1"))  # [1] 40 20

} # viz.corysTrem2Model
#------------------------------------------------------------------------------------------------------------------------
queryDB <- function()
{
   library(RPostgreSQL)
   dbConnection <- dbConnect(PostgreSQL(), user="trena", password="trena", host="khaleesi", dbname="brain_hint_20")
   tbl.regions <- data.frame(chrom="chr6", start=41200271, end=41200306, stringsAsFactors=FALSE)
   tbl.hits <- trenaSGM:::.multiQueryFootprints(dbConnection, tbl.regions)

} # queryDB
#------------------------------------------------------------------------------------------------------------------------
find.tf.bindingSites <- function(tf)
{
   mdb.tf <- query(MotifDb, c("hsapiens", "jaspar2018", "BHLHE41"))
   stopifnot(length(mdb.tf) > 0)
   motifs <- mcols(mdb.tf)$providerId

   print(load(system.file(package="trenaSGM", "extdata", "enhancers.TREM2.RData")))  # tbl.enhancers
   bigRegion <- with(tbl.enhancers, sprintf("%s:%d-%d", unique(chrom), min(start)-2000, max(end)+2000))
   mm <- MotifMatcher("hg38", pfms=as.list(mdb.tf))

   tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.enhancers, pwmMatchMinimumAsPercentage=80)
   track <- DataFrameAnnotationTrack("MA0636.1", tbl.motifs[, c("chrom", "motifStart", "motifEnd")], color="green")
   displayTrack(igv, track)

   tbl.bhlhe41 <- subset(tbl.footprints, tf=="BHLHE41")
   track <- DataFrameAnnotationTrack("fimo-BHLHE41", tbl.bhlhe41[, c("chrom", "fp_start", "fp_end")], color="blue")
   displayTrack(igv, track)


} # find.tf.bindingSites
#------------------------------------------------------------------------------------------------------------------------
test_tfPoolOption <- function()
{
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
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup, all known tfs.
     #------------------------------------------------------------

   builder <- FastFootprintDatabaseModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
   x <- build(builder)

   build.spec$tfPool <- c("IKZF1", "TAL1")
   builder <- FastFootprintDatabaseModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
   x <- build(builder)
   checkEquals(sort(x$model$gene), sort(build.spec$tfPool))
   checkEquals(sort(unique(x$regulatoryRegions$geneSymbol)), sort(build.spec$tfPool))

} # test_tfPoolOption
#------------------------------------------------------------------------------------------------------------------------
# the footprint databases (as of jun 20 2018) name footprint/motif-associated TFs by their HGNC gene symbols
# to support trena analysis on expression data expressed with ensg (ensembl gene) IDs, we need to translate
# those footprint tfs from geneSymbol to ensembl ids.  develop and test that here
#
# use this mapping
#    SYMBOL         ENSEMBL
# 1   PRDM2 ENSG00000116731
# 2   EOMES ENSG00000163508
# 3    ELF2 ENSG00000109381
# 4   TREM2 ENSG00000095970
# 5   TBPL1 ENSG00000028839
# 6    IRF5 ENSG00000128604
# 7    SPI1 ENSG00000066336
# 8   PRDM4 ENSG00000110851
# 9    EBF4 ENSG00000088881
# 10   LYL1 ENSG00000104903
# 11  IKZF1 ENSG00000185811
# 12   TFEC ENSG00000105967
# 13  DMRT2 ENSG00000173253
test_modelExpressionDataWithEnsgIDs <- function()
{
   printf("--- test_ensgIDs")

      # create a smallish expression matrix with ENSEMBL gene ids as rownames
      # starting with mtx.tcx for which we have a good gene-symbol-centric model

   tfs.highRanked <- c("IKZF1", "IRF5",  "LYL1",  "SPI1",  "TFEC")
   set.seed(17)
   tfs.randomUnranked <- allKnownTFs()[sample(1:1000, 100)]
   all.genes <- intersect(rownames(mtx), c("TREM2", tfs.highRanked, tfs.randomUnranked))
   suppressMessages(tbl.map <- select(org.Hs.eg.db, keys=all.genes, keytype="SYMBOL", columns="ENSEMBL"))
   dups <- which(duplicated(tbl.map$SYMBOL))
   if(length(dups) > 0)
      tbl.map <- tbl.map[-dups,]

   checkTrue(all(tbl.map$SYMBOL %in% rownames(mtx)))
   mtx.sub <- mtx[tbl.map$SYMBOL,]
   indices <- match(tbl.map$SYMBOL, rownames(mtx.sub))
   rownames(mtx.sub) <- tbl.map$ENSEMBL[indices]

   genome <- "hg38"
   targetGene <- "ENSG00000095970"  # "TREM2"
   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186

      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.ensemblTest",
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx.sub,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfMapping="MotifDB",
                      tfPool=allKnownTFs(identifierType="ensemblGeneID"),
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup, all known tfs.
     #------------------------------------------------------------

   builder <- FastFootprintDatabaseModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
   x <- build(builder)
   checkTrue("ENSG00000185811" %in% x$model$gene)

} # test_ensgIDs
#------------------------------------------------------------------------------------------------------------------------
test_noGenesAboveExpressionCorrelationThreshold <- function()
{
   printf("--- test_noGenesAboveExpressionCorrelationThreshold")

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
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPool=allKnownTFs(),
                      tfPrefilterCorrelation=0.99,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup, all known tfs.
     #------------------------------------------------------------

   builder <- FastFootprintDatabaseModelBuilder(genome, targetGene, build.spec, quiet=TRUE)
   x <- build(builder)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   checkEquals(nrow(x$model), 0)
   checkEquals(nrow(x$regulatoryRegions), 0)

} # test_noGenesAboveExpressionCorrelationThreshold()
#------------------------------------------------------------------------------------------------------------------------
test_stagedExecution <- function()
{
   printf("--- test_stagedExecution")

   genome <- "hg38"
   targetGene <- "TREM2"
   targetGene.ensembl <- "ENSG00000095970"
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
                      geneSymbol=targetGene.ensembl,
                      tss=tss,
                      matrix=mtx.ensembl,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      tfPool=allKnownTFs(identifierType="ensemblGeneID"),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.2,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   # stageDir <- "stage"
   stageDir <- tempdir()

   fpBuilder <- FastFootprintDatabaseModelBuilder(genome, targetGene.ensembl, build.spec, quiet=FALSE,
                                              stagedExecutionDirectory=stageDir)
   fp.filename <- staged.fast.build(fpBuilder, stage="find.footprints")
   checkTrue(file.exists(fp.filename))

   fp.tfMapped.filename <- staged.fast.build(fpBuilder, stage="associateTFs")
   checkTrue(file.exists(fp.tfMapped.filename))

   models.filename <- staged.fast.build(fpBuilder, stage="build.models")

   load(models.filename)
   checkEquals(names(tbls), c("model", "regulatoryRegions"))
   tbl.regions <- tbls$regulatoryRegions
   tbl.model <- tbls$model
   tbl.model <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE),]
   top.tfs <- head(tbl.model$gene)
       # allow for stochasticity, check just 3 of those 6 top.tfs
   checkTrue("ENSG00000185811" %in% head(tbl.model$gene))
   tbl.regions <- tbls$regulatoryRegions
   tfs.in.regions <- unique(tbl.regions$geneSymbol)
      # we keep only regions with tfs in model.  check that
   checkTrue(all(tfs.in.regions %in% tbl.model$gene))

} # test_stagedExecution
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
