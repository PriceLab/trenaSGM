# test-10.R: experiment with 10 ensembl genes
#------------------------------------------------------------------------------
library(RUnit)
library(trenaSGM)
library(MotifDb)
library(motifStack)
library(biomaRt)
library(tibble)
library(BiocParallel)
#library(BatchJobs)
library(futile.logger)
library(RPostgreSQL)
#----------------------------------------------------------------------------------------------------
if(!exists("mtx"))
    load("Scaled_Winsorized_MayoRNAseq_TCX.mtx.RData")

if(!exists("enhancer.list"))
    load("enhancer.list.RData")

if(!exists("tbl.geneInfo"))
    load("tbl.geneInfo.RData")

incremental.data.directory <- "./results"
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
runSGM <- function(spec)
{
   printf("-- runSGM(%s)", spec$targetGene)
   targetGene <- spec$targetGene
   tbl.geneLoc <- tbl.geneInfo[targetGene,]
   geneSymbol <- tbl.geneLoc$geneSymbol
   chromosome <- tbl.geneLoc$chrom
   tss <- tbl.geneLoc$tss

   flog.info(sprintf("runSGM on %s, %s", geneSymbol, targetGene))
   flog.info(sprintf("  assigning regions in mode %s", spec$regionsMode))

   tbl.regions <- switch(spec$regionsMode,
                         "enhancers" = {
                             enhancer.list[[targetGene]][, c("chrom", "start", "end")]
                         },
                         "tiny" = {
                            data.frame(chrom=chromosome, start=tss-1000,
                                       end=tss+1000, stringsAsFactors=FALSE)
                         })
   flog.info(sprintf("  assigning regions in mode %s, %d regions", spec$regionsMode, nrow(tbl.regions)))

   genome <- "hg38"

   flog.info("constructing trenaSGM")
   sgm <- trenaSGM(genome, targetGene, quiet=FALSE)
   flog.info("   after ctor")

   build.spec <- list(title=targetGene,
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=spec$db,
                      motifDiscovery="builtinFimo",
                      tfMapping=c("MotifDB","TFClass"),
                      tfPrefilterCorrelation=spec$correlationThreshold,
                      tfPool=allKnownTFs(identifierType="ensemblGeneID"),
                      orderModelByColumn="pearsonCoeff",
                      solverNames=spec$solvers
                      )

   flog.info("   after build.spec assignement")

   strategies <- list(one=build.spec)

   flog.info("calling calculate")

   model <- calculate(sgm, strategies)$one

   flog.info("calculate complete")


   filename <- sprintf("%s/%s-%s.RData", incremental.data.directory, targetGene,
                       tbl.geneInfo[targetGene, "geneSymbol"])

   flog.info(sprintf("saving model for %s (%s): %d tfs", targetGene, tbl.geneInfo[targetGene,]$geneSymbol,
                     nrow(model$model)))

   save(model, file=filename)
   flog.info("save complete")

   return(model)

} # runSGM
#----------------------------------------------------------------------------------------------------
test_runSGM <- function()
{
   ad.genes <-c("TREM2", "CR1", "BIN1", "CD2AP", "EPHA1", "CLU", "MS4A6A", "PICALM", "ABCA7",
                "CD33", "HLA-DRB5", "HLA-DRB1", "PTK2B", "SORL1", "SLC24A4", "RIN3", "DSG2",
                "INPP5D", "MEF2C", "NME8", "ZCWPW1", "CELF1", "FERMT2", "CASS4", "APOE", "TOMM40")



   mini.recipes <- lapply(ad.genes,
                          function(gene)
                              list(targetGene=targetGenes[[gene]],
                                   regionsMode="tiny",
                                   correlationThreshold=0.5,
                                   solvers= c("pearson", "spearman"),
                                   dbs="brain_hint_20"))
   names(mini.recipes) <- as.character(targetGenes)

  # spec.trem2 <- list(targetGene=targetGenes[["TREM2"]],
  #                      regionsMode="tiny",
  #                    correlationThreshold=0.7,
  #                    solvers= c("pearson", "spearman"),
  #                    dbs="brain_hint_20",
  #                    gene.info=tbl.geneInfo
  #                    )
   x.trem2 <- runSGM(mini.recipes[["TREM2"]])

   x.trem2 <- runSGM(spec.trem2)
   x.rin3 <- runSGM(spec.rin3)
   x.mef2c <- runSGM(spec.mef2c)

   dim(x.trem2$model)
   dim(x.rin3$one$model)
   dim(x.mef2c$one$model)
   checkTrue(all(x.mef2c$one$regulatoryRegions$geneSymbol %in% x.mef2c$one$model$gene))
   checkTrue(all(x.rin3$one$regulatoryRegions$geneSymbol %in% x.rin3$one$model$gene))
   checkTrue(all(x.trem2$one$regulatoryRegions$geneSymbol %in% x.trem2$one$model$gene))

   fast.genes <- c(1,5,7,10,11,12)   # relatively few footprints
   slow.genes <- c(2,3,4)
   system.time(x <- lapply(mini.recipes[1:12], runSGM))
   system.time(y <- bptry({
       result <- bplapply(mini.recipes[goi],
                runSGM,
                #BPPARAM=SerialParam()
                BPPARAM=MulticoreParam(workers=2) #, stop.on.error=FALSE)
                )
       })) #, error=identify))

   bpok(y)

   goi <- c(1,5)
   goi <- fast.genes
   goi <- slow.genes
   goi <- seq_len(length(ad.genes))
   param <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir="./")
   sgm.out <- bptry(bplapply(mini.recipes[goi],runSGM, BPPARAM=param))
   bpok(sgm.out)



    # riptide    2   812
    # khaleesi   2   558     10 genes
    # khaleesi   2   707     24 genes
    # khaleesi   2   710     24 genes
    # khaleesi   4    50     ~6 genes, filaed on "error in serialize, ignoring SIGPIPE signal"
    # khaleesi   2   60
    # khaleesi   2   562
    # khaleesi   3   461
    # khaleesi   4   331     5/6 successful

} # test_runSGM
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
   tbl.geneLoc <- tbl.geneInfo[targetGene,]
   geneSymbol <- tbl.geneLoc$geneSymbol
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
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.2,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   stageDir <- "stage" # tempdir()
   fpBuilder <- FootprintDatabaseModelBuilder(genomeName, targetGene, build.spec, quiet=FALSE,
                                              stagedExecutionDirectory=stageDir)
   fp.filename <- staged.build(fpBuilder, stage="find.footprints")
   checkTrue(file.exists(fp.filename))

} # runStagedSGM.footprints
#----------------------------------------------------------------------------------------------------
test_runStagedSGM.footprints <- function()
{
   short.specs <- lapply(ad.genes,
                          function(gene)
                            list(targetGene=targetGenes[[gene]],
                                 geneSymbol=gene,
                                 regionsMode="tiny",
                                 correlationThreshold=0.5,
                                 solvers= c("pearson", "spearman"),
                                 dbs="brain_hint_20"))
   names(short.specs) <- as.character(targetGenes)

   runStagedSGM.footprints(short.specs[[1]])
   x <- bplapply(short.specs, runStagedSGM.footprints)

} # test_runStagedSGM.footprints
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
