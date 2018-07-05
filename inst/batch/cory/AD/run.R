# test-10.R: experiment with 10 ensembl genes
#------------------------------------------------------------------------------
library(RUnit)
library(trenaSGM)
library(MotifDb)
library(motifStack)
library(biomaRt)
library(tibble)
library(BiocParallel)
library(futile.logger)
library(RPostgreSQL)
library(org.Hs.eg.db)
#----------------------------------------------------------------------------------------------------
if(!exists("mtx"))
    load("Scaled_Winsorized_MayoRNAseq_TCX.mtx.RData")

if(!exists("enhancer.list"))
    load("enhancer.list.RData")

if(!exists("tbl.geneInfo"))
    load("tbl.geneInfo.RData")

incremental.data.directory <- "./results"
#----------------------------------------------------------------------------------------------------
SOLVERS <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
DBS     <- c("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16")
GENE_EXPRESSION_MINIMUM_CORRELATION <- 0.4
STAGEDIR <- "STAGE.AD"
#----------------------------------------------------------------------------------------------------
ad.genes <-c("TREM2", "CR1", "BIN1", "CD2AP", "EPHA1", "CLU", "MS4A6A", "PICALM", "ABCA7",
             "CD33", "HLA-DRB5", "HLA-DRB1", "PTK2B", "SORL1", "SLC24A4", "RIN3",
             "INPP5D", "MEF2C",  "ZCWPW1", "CELF1", "FERMT2", "CASS4", "APOE", "TOMM40")
          # "NME8", "DSG2",

if(!exists("targetGenes")){   # a named list, geneSymbols as names, ensg as content
   indices <- match(ad.genes, tbl.geneInfo$geneSymbol)
   deleters <- which(is.na(indices))
   if(length(deleters) > 0){
      ad.genes <- ad.genes[-deleters]
      indices <- indices[-deleters]
      }
   targetGenes <- rownames(tbl.geneInfo[indices,])
   names(targetGenes) <- ad.genes
   }
#----------------------------------------------------------------------------------------------------
runStagedSGM.footprints <- function(short.spec)
{
   required.fields <- c("targetGene", "geneSymbol", "regionsMode", "correlationThreshold", "solvers", "dbs")
   missing.fields <- setdiff(required.fields, names(short.spec))
   if(length(missing.fields) > 0){
      msg <- sprintf("runStagedSGM.footprings finds fields missing in short.spec: %s", paste(missing.fields, collapse=", "))
      stop(msg)
      }

   printf("-- runSGM(%s)", short.spec$targetGene)

   genomeName <- "hg38"
   targetGene <- short.spec$targetGene
   geneSymbol <- short.spec$geneSymbol
   tbl.geneLoc <- tbl.geneInfo[targetGene,]
   chromosome <- tbl.geneLoc$chrom
   tss <- tbl.geneLoc$tss

   tbl.regions <- switch(short.spec$regionsMode,
                         "enhancers" = {enhancer.list[[targetGene]][, c("chrom", "start", "end")]},
                         "tiny" = {data.frame(chrom=chromosome, start=tss-1000, end=tss+1000, stringsAsFactors=FALSE)
                         })

   tbl.regions <- switch(short.spec$regionsMode,
                         "enhancers" = {enhancer.list[[targetGene]][, c("chrom", "start", "end")]},
                         "tiny" = {data.frame(chrom=chromosome, start=tss-1000, end=tss+1000, stringsAsFactors=FALSE)
                         })

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      geneSymbol=geneSymbol,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=DBS,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.2,
                      orderModelByColumn="rfScore",
                      solverNames=SOLVERS)


   fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene, build.spec, quiet=FALSE,
                                              stagedExecutionDirectory=STAGEDIR)
   fp.filename <- staged.build(fpBuilder, stage="find.footprints")
   checkTrue(file.exists(fp.filename))

} # runStagedSGM.footprints
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.footprints <- function()
{
   short.specs <- lapply(ad.genes,
                          function(gene)
                            list(targetGene=targetGenes[[gene]],
                                 geneSymbol=gene,
                                 regionsMode="tiny",
                                 correlationThreshold=GENE_EXPRESSION_MINIMUM_CORRELATION,
                                 solvers=SOLVERS,
                                 dbs=DBS))

   names(short.specs) <- as.character(targetGenes)

   runStagedSGM.footprints(short.specs[[1]])
    # x is true/false, named by target gene
   x.fp <- bplapply(short.specs, runStagedSGM.footprints)
   #x.tfMap <- bplapply(short.specs, runStagedSGM.tfMapping)

} # do.runStagedSGM.footprints
#----------------------------------------------------------------------------------------------------
runStagedSGM.associateTFs <- function(short.spec)
{
   required.fields <- c("targetGene", "geneSymbol", "regionsMode", "correlationThreshold", "solvers", "dbs")
   missing.fields <- setdiff(required.fields, names(short.spec))
   if(length(missing.fields) > 0){
      msg <- sprintf("runStagedSGM.footprings finds fields missing in short.spec: %s", paste(missing.fields, collapse=", "))
      stop(msg)
      }

   printf("-- runSGM(%s)", short.spec$targetGene)

   genomeName <- "hg38"
   targetGene <- short.spec$targetGene
   geneSymbol <- short.spec$geneSymbol
   tbl.geneLoc <- tbl.geneInfo[targetGene,]
   chromosome <- tbl.geneLoc$chrom
   tss <- tbl.geneLoc$tss

   tbl.regions <- switch(short.spec$regionsMode,
                         "enhancers" = {enhancer.list[[targetGene]][, c("chrom", "start", "end")]},
                         "tiny" = {data.frame(chrom=chromosome, start=tss-1000, end=tss+1000, stringsAsFactors=FALSE)
                         })

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      geneSymbol=geneSymbol,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.2,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene, build.spec, quiet=FALSE,
                                              stagedExecutionDirectory=STAGEDIR)
   fp.filename <- staged.build(fpBuilder, stage="associateTFs")
   checkTrue(file.exists(fp.filename))

} # runStagedSGM.associateTFs
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.associateTFs <- function()
{
   short.specs <- lapply(ad.genes,
                          function(gene)
                            list(targetGene=targetGenes[[gene]],
                                 geneSymbol=gene,
                                 regionsMode="tiny",
                                 correlationThreshold=GENE_EXPRESSION_MINIMUM_CORRELATION,
                                 solvers=SOLVERS,
                                 dbs=DBS))
   names(short.specs) <- as.character(targetGenes)

   runStagedSGM.associateTFs(short.specs[[1]])
   x.tfMap <- bplapply(short.specs, runStagedSGM.associateTFs)

} # do.runStagedSGM.associateTFs
#----------------------------------------------------------------------------------------------------
runStagedSGM.buildModels <- function(short.spec)
{
   required.fields <- c("targetGene", "geneSymbol", "regionsMode", "correlationThreshold", "solvers", "dbs")
   missing.fields <- setdiff(required.fields, names(short.spec))
   if(length(missing.fields) > 0){
      msg <- sprintf("runStagedSGM.footprings finds fields missing in short.spec: %s", paste(missing.fields, collapse=", "))
      stop(msg)
      }

   printf("-- runSGM(%s)", short.spec$targetGene)

   genomeName <- "hg38"
   targetGene <- short.spec$targetGene
   geneSymbol <- short.spec$geneSymbol
   tbl.geneLoc <- tbl.geneInfo[targetGene,]
   chromosome <- tbl.geneLoc$chrom
   tss <- tbl.geneLoc$tss

   tbl.regions <- switch(short.spec$regionsMode,
                         "enhancers" = {enhancer.list[[targetGene]][, c("chrom", "start", "end")]},
                         "tiny" = {data.frame(chrom=chromosome, start=tss-1000, end=tss+1000, stringsAsFactors=FALSE)
                         })

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      geneSymbol=geneSymbol,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="builtinFimo",
                      tfPool=short.spec$tfPool,
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.2,
                      orderModelByColumn="rfScore",
                      solverNames=short.spec$solvers)


   fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene, build.spec, quiet=FALSE,
                                              stagedExecutionDirectory=STAGEDIR)
   fp.filename <- staged.build(fpBuilder, stage="build.models")
   checkTrue(file.exists(fp.filename))

} # runStagedSGM.buildModels
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.buildModels <- function()
{
   ensembl.tfPool <- allKnownTFs(identifierType="ensemblGeneID")

   short.specs <- lapply(ad.genes,
                          function(gene)
                            list(targetGene=targetGenes[[gene]],
                                 geneSymbol=gene,
                                 regionsMode="tiny",
                                 correlationThreshold=GENE_EXPRESSION_MINIMUM_CORRELATION,
                                 tfPool=ensembl.tfPool,
                                 solvers=SOLVERS,
                                 dbs=DBS))

   names(short.specs) <- as.character(targetGenes)

   runStagedSGM.buildModels(short.specs[[1]])
   x.tfMap <- bplapply(short.specs, runStagedSGM.buildModels)

} # do.runStagedSGM.buildModels
#----------------------------------------------------------------------------------------------------
run <- function()
{

    #    file <- "Scaled_Winsorized_MayoRNAseq_TCX.csv"
    #     # file <- "test.csv"
    #     # file <- "test2.csv"   # cut -d , -f 1-5 test.csv > test2.csv
    #    tbl.counts <- read.csv(file)
    #    rownames(tbl.counts) <- tbl.counts$X
    #    tbl.counts <- tbl.counts[, -1]
    #    mtx <- as.matrix(t(tbl.counts))
    #    stopifnot(class(mtx[1,1]) == "numeric")
    #    dim(mtx)
    #    save(mtx, file="Scaled_Winsorized_MayoRNAseq_TCX.mtx.RData")
    #    print(load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData")))
    #    dim(mtx)
   load("Scaled_Winsorized_MayoRNAseq_TCX.mtx.RData")
   print(load("enhancer.list.RData"))
   # get a list of all the tss's
   all.tss <- getTSS.ensembl(row.names(mtx))
   tbl.tss <- all.tss
   rownames(tbl.tss) <- NULL
   save(tbl.tss, file="tbl.tss-allEnsemblIDsInMatrix.RData")

    # test in direct serial mode
   x <- runSGM(rownames(mtx)[1])

    # test the above function with 10 genes
   names.10 <- head(row.names(mtx), 10)
   num.cores <- 3
   register(MulticoreParam(workers = num.cores, stop.on.error = FALSE, log = TRUE), default = TRUE)
   system.time(test.SGM.3 <- bplapply(names.10, runSGM))

   #   stderr and stdout:
   #
   #   <simpleError in .runTrenaWithRegulatoryRegions(obj@genomeName, s$tfPool, obj@targetGene,
   #            tbl.fp, s$matrix, s$tfPrefilterCorrelation, s$solverNames,     obj@quiet):
   #            NO genes have expression >= 0.400000 correlated with targetGene 'ENSG00000227232'>
   #
   #   Error: BiocParallel errors
   #      element index: 2, 3, 4, 5, 6, 7, ...
   #     first error: database disk image is malformed

} # run
#------------------------------------------------------------------------------------------------------------------------
demo.bplapply.bug <- function()
{
    library(org.Hs.eg.db)
    lookup <- function(geneSymbols){
       tbl.map <- select(org.Hs.eg.db, keys=geneSymbols, keytype="SYMBOL", columns="ENSEMBL")
       }

    tf.entrezIDs <- unique(unlist(get("GO:0003700", envir=org.Hs.egGO2ALLEGS)), use.names=FALSE)
    tf.geneSymbols <- unique(unlist(mget(tf.entrezIDs, envir=org.Hs.egSYMBOL), use.names=FALSE))

    geneSymbols.list <- list(g1=tf.geneSymbols, g2=tf.geneSymbols, g3=tf.geneSymbols)
    x <- bplapply(geneSymbols.list, lookup)
        
} # demo.bplapply.bug
#------------------------------------------------------------------------------------------------------------------------
demo.bplapply.fix <- function()
{
  library(AnnotationDbi)
  library(BiocParallel)
  library(org.Hs.eg.db)
  
  lookup <- function(geneSymbols, dbfile) {
     db <- AnnotationDbi::loadDb(dbfile)
     on.exit(RSQLite::dbDisconnect(dbconn(db)))
     tbl.map <- AnnotationDbi::select(db, keys=geneSymbols, keytype="SYMBOL", columns="ENSEMBL")
     }
  
  # tf.entrezIDs <- unique(unlist(get("GO:0003700", envir=org.Hs.egGO2ALLEGS)), use.names=FALSE)
  # tf.geneSymbols <- unique(unlist(mget(tf.entrezIDs, envir=org.Hs.egSYMBOL), use.names=FALSE))
  tf.geneSymbols <- unique( select(org.Hs.eg.db, "GO:0003700", "SYMBOL", "GOALL")$SYMBOL)
  geneSymbols.list <- list(g1=tf.geneSymbols, g2=tf.geneSymbols, g3=tf.geneSymbols)
  x <- bplapply(geneSymbols.list, lookup, dbfile(org.Hs.eg.db))

} # demo.bplapply.fix
#------------------------------------------------------------------------------------------------------------------------

