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

#----------------------------------------------------------------------------------------------------
SOLVERS <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
DBS     <- c("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16")
GENE_EXPRESSION_MINIMUM_CORRELATION <- 0.2
#---------------------------------------------------------------------------------------------------
OUTPUTDIR <- "OUTPUT.3.v2"
fp.logDir <- file.path(OUTPUTDIR, "log.fp")
tfMapping.logDir <- file.path(OUTPUTDIR, "log.tfMap")
model.logDir <- file.path(OUTPUTDIR, "log.model")
if(!file.exists(OUTPUTDIR)) dir.create(OUTPUTDIR)
if(!file.exists(fp.logDir)) dir.create(fp.logDir)
if(!file.exists(tfMapping.logDir)) dir.create(tfMapping.logDir)
if(!file.exists(model.logDir)) dir.create(model.logDir)
#---------------------------------------------------------------------------------------------------
DBHOST <- "khaleesi.systemsbiology.net"
#----------------------------------------------------------------------------------------------------
basic.build.spec <- list(title="footprint-based-tf-model-builder-for-ADgenes",
                         type="footprint.database",
                         stageDirectory=OUTPUTDIR,
                         genomeName="hg38",
                         matrix=mtx,
                         db.host=DBHOST,
                         databases=DBS,
                         annotationDbFile=dbfile(org.Hs.eg.db),
                         motifDiscovery="builtinFimo",
                         tfPool=allKnownTFs(identifierType="ensemblGeneID"),
                         tfMapping=c("MotifDB", "TFClass"),
                         tfPrefilterCorrelation=GENE_EXPRESSION_MINIMUM_CORRELATION,
                         correlationThreshold=GENE_EXPRESSION_MINIMUM_CORRELATION,
                         orderModelByColumn="rfScore",
                         solverNames=SOLVERS,
                         solvers=SOLVERS,
                         dbs=DBS)

#----------------------------------------------------------------------------------------------------
ad.genes <-c ("TREM2", "CR1", "BIN1", "CD2AP", "EPHA1", "CLU", "MS4A6A", "PICALM", "ABCA7",
              "CD33", "HLA-DRB5", "HLA-DRB1", "PTK2B", "SORL1", "SLC24A4", "RIN3",
              "INPP5D", "MEF2C",  "ZCWPW1", "CELF1", "FERMT2", "CASS4", "APOE", "TOMM40")
          # "NME8", "DSG2",

goi <- ad.genes
goi <- rownames(mtx)[1:3]
#----------------------------------------------------------------------------------------------------
if(!exists("targetGenes")){   # a named list, geneSymbols as names, ensg as content
   indices <- match(goi, rownames(tbl.geneInfo))
   deleters <- which(is.na(indices))
   if(length(deleters) > 0){
      goi <- goi[-deleters]
      indices <- indices[-deleters]
      }
   targetGenes <- rownames(tbl.geneInfo[indices,])
   names(targetGenes) <- goi
   }
#----------------------------------------------------------------------------------------------------
runStagedSGM.footprints <- function(short.spec)
{
   required.fields <- c("targetGene", "geneSymbol", "regionsMode")
   missing.fields <- setdiff(required.fields, names(short.spec))
   if(length(missing.fields) > 0){
      msg <- sprintf("runStagedSGM.footprints finds fields missing in short.spec: %s", paste(missing.fields, collapse=", "))
      stop(msg)
      }

   printf("-- runSGM(%s)", short.spec$targetGene)

   targetGene <- short.spec$targetGene
   geneSymbol <- short.spec$geneSymbol
   tbl.geneLoc <- tbl.geneInfo[targetGene,]
   chromosome <- tbl.geneLoc$chrom
   tss <- tbl.geneLoc$tss

   tbl.regions <- switch(short.spec$regionsMode,
                         "enhancers" = {enhancer.list[[targetGene]][, c("chrom", "start", "end")]},
                         "tiny" = {data.frame(chrom=chromosome, start=tss-1000, end=tss+1000, stringsAsFactors=FALSE)
                         })


   spec <- basic.build.spec
   spec$regions <- tbl.regions
   spec$targetGene <- targetGene
   spec$tss=tss
   spec$geneSymbol=geneSymbol
   spec

   fpBuilder <- FootprintDatabaseModelBuilder(spec$genomeName, targetGene, spec,
                                              quiet=FALSE,
                                              stagedExecutionDirectory=spec$stageDir)
   fp.filename <- staged.build(fpBuilder, stage="find.footprints")
   checkTrue(file.exists(fp.filename))

} # runStagedSGM.footprints
#----------------------------------------------------------------------------------------------------
runStagedSGM.associateTFs <- function(short.spec)
{
   required.fields <- c("targetGene", "geneSymbol") #, "regionsMode", "correlationThreshold", "solvers", "dbs")
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


   
   spec <- basic.build.spec
   spec$targetGene <- targetGene
   spec$tss=tss
   spec$regions <- tbl.regions
   spec$geneSymbol=geneSymbol

   fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene, spec, quiet=FALSE,
                                              stagedExecutionDirectory=OUTPUTDIR)
   fp.filename <- staged.build(fpBuilder, stage="associateTFs")
   checkTrue(file.exists(fp.filename))

} # runStagedSGM.associateTFs
#----------------------------------------------------------------------------------------------------
runStagedSGM.buildModels <- function(short.spec)
{
   required.fields <- c("targetGene", "geneSymbol", "regionsMode") # , "correlationThreshold", "solvers", "dbs")
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
   
   spec <- basic.build.spec
   spec$targetGene <- targetGene
   spec$tss=tss
   spec$regions <- tbl.regions
   spec$geneSymbol=geneSymbol
   #spec$correlationThreshold <-
   #spec$solvers <- 
   #spec$dbs <- 

   fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene, spec, quiet=FALSE,
                                              stagedExecutionDirectory=OUTPUTDIR)
   fp.filename <- staged.build(fpBuilder, stage="build.models")
   checkTrue(file.exists(fp.filename))

} # runStagedSGM.buildModels
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.footprints <- function()
{
   short.specs <- lapply(goi,
                          function(gene)
                            list(targetGene=gene,
                                 geneSymbol=tbl.geneInfo[gene,"geneSymbol"],
                                 regionsMode="tiny"))
   names(short.specs) <- as.character(goi)

     # runStagedSGM.footprints(short.specs[[1]])
   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=fp.logDir, level="INFO")
   results.fp <<- bptry({bplapply(short.specs, runStagedSGM.footprints, BPPARAM=bp.params)})

} # do.runStagedSGM.footprints
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.associateTFs <- function()
{
   short.specs <- lapply(goi,
                          function(gene)
                            list(targetGene=gene,
                                 geneSymbol=tbl.geneInfo[gene,"geneSymbol"],
                                 regionsMode="tiny"))
                                 #correlationThreshold=GENE_EXPRESSION_MINIMUM_CORRELATION,
                                 #solvers=SOLVERS,
                                 #dbs=DBS))
   names(short.specs) <- as.character(goi)

     # runStagedSGM.associateTFs(short.specs[[3]])
   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=tfMapping.logDir, level="INFO")
   results.tfMap <<- bptry({bplapply(short.specs, runStagedSGM.associateTFs, BPPARAM=bp.params)})



} # do.runStagedSGM.associateTFs
#----------------------------------------------------------------------------------------------------
do.runStagedSGM.buildModels <- function()
{
   #ensembl.tfPool <- allKnownTFs(identifierType="ensemblGeneID")

   short.specs <- lapply(goi,
                          function(gene)
                            list(targetGene=gene,
                                 geneSymbol=tbl.geneInfo[gene,"geneSymbol"],
                                 regionsMode="tiny"))
                                 #correlationThreshold=GENE_EXPRESSION_MINIMUM_CORRELATION,
                                 #tfPool=ensembl.tfPool,
                                 #solvers=SOLVERS,
                                 #dbs=DBS))

   names(short.specs) <- as.character(goi)

      # runStagedSGM.buildModels(short.specs[[1]])
   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=model.logDir, level="INFO")
   results.buildModels <<- bptry({bplapply(short.specs, runStagedSGM.buildModels, BPPARAM=bp.params)})

} # do.runStagedSGM.buildModels
#----------------------------------------------------------------------------------------------------
printf("starting run of %d goi", length(goi))
do.runStagedSGM.footprints()
do.runStagedSGM.associateTFs()
do.runStagedSGM.buildModels()

