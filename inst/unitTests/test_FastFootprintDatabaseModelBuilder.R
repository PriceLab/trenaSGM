library(RUnit)
library(trenaSGM)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
Sys.setlocale("LC_ALL", "C")

if(!exists("mtx")){
    load(system.file(package="trenaSGM", "extdata", "mayo.tcx.new.RData"))
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
   suppressMessages(tbl.map <-  select(org.Hs.eg.db, keys=ensembl.ids, keytype="ENSEMBL", columns=c("SYMBOL", "ENSEMBL")))
   tbl.map <- tbl.map[-which(duplicated(tbl.map$ENSEMBL)),]   # two dups
   rownames(tbl.trena) <- tbl.map$SYMBOL
   }

#----------------------------------------------------------------------------------------------------
targetGene <- "TREM2"
targetGene.ensembl <- "ENSG00000095970"
genome <- "hg38"
targetGene <- "TREM2"
chromosome <- "chr6"
upstream <- 2000
downstream <- 200
tss <- 41163186

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_stagedExecution()
   test_withManyEnhancers()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

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

   fpBuilder <- FastFootprintDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)

   checkTrue("FastFootprintDatabaseModelBuilder" %in% is(fpBuilder))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_stagedExecution <- function()
{
   printf("--- test_stagedExecution")

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
                      tfPrefilterCorrelation=0.1,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

     #------------------------------------------------------------
     # use the above build.spec: a small region, high correlation
     # required, MotifDb for motif/tf lookup
     #------------------------------------------------------------

   # stageDir <- "stage"
   stageDir <- "stage" # tempdir()

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
test_withManyEnhancers <- function()
{
   printf("--- test_withManyEnhancers")

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      regions=tbl.enhancers,    # 16kb
                      geneSymbol=targetGene.ensembl,
                      tss=tss,
                      matrix=mtx.ensembl,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
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

} # test_withManyEnhancers
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
