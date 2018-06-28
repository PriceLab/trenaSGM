# test-10.R: experiment with 10 ensembl genes
#------------------------------------------------------------------------------
library(RUnit)
library(trenaSGM)
library(MotifDb)
library(motifStack)
library(biomaRt)
library(tibble)
library(BiocParallel)
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
runSGM <- function(spec)
{
   print(1)
   targetGene <- spec$targetGene
   print(2)
   tbl.geneLoc <- tbl.geneInfo[targetGene,]
   print(3)
   chromosome <- tbl.geneLoc$chrom
   print(4)
   tss <- tbl.geneLoc$tss

   print(5)
   tbl.regions <- switch(spec$regionsMode,
                         "enhancers" = {
                             enhancer.list[[targetGene]][, c("chrom", "start", "end")]
                         },
                         "tiny" = {
                            data.frame(chrom=chromosome, start=tss-1000,
                                       end=tss+1000, stringsAsFactors=FALSE)
                         })

   print(6)
   genome <- "hg38"
   print(7)
   sgm <- trenaSGM(genome, targetGene, quiet=FALSE)
   print(8)

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

   print(9)
   strategies <- list(one=build.spec)
   print(10)
   model <- calculate(sgm, strategies)$one
   print(11)

   filename <- sprintf("%s/%s-%s.RData", incremental.data.directory, targetGene,
                       tbl.geneInfo[targetGene, "geneSymbol"])
   printf("saving model for %s (%s): %d tfs", targetGene, tbl.geneInfo[targetGene,]$geneSymbol,
          nrow(model$model))
   save(model, file=filename)

   return(model)

} # runSGM
#----------------------------------------------------------------------------------------------------
test_runSGM <- function()
{
   ad.genes <-c("TREM2", "CR1", "BIN1", "CD2AP", "EPHA1", "CLU", "MS4A6A", "PICALM", "ABCA7",
                "CD33", "HLA-DRB5", "HLA-DRB1", "PTK2B", "SORL1", "SLC24A4", "RIN3", "DSG2",
                "INPP5D", "MEF2C", "NME8", "ZCWPW1", "CELF1", "FERMT2", "CASS4", "APOE", "TOMM40")


   indices <- match(ad.genes, tbl.geneInfo$geneSymbol)
   targetGenes <- rownames(tbl.geneInfo[indices,])
   names(targetGenes) <- ad.genes

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
   goi <- c(1,5)
   goi <- fast.genes
   goi <- slow.genes
   goi <- seq_len(length(ad.genes))
   system.time(y <- bptry({
       result <- bplapply(mini.recipes[goi],
                runSGM,
                #BPPARAM=SerialParam()
                BPPARAM=MulticoreParam(stop.on.error=FALSE)
                )
       })) #, error=identify))

    bpok(y)


} # test_runSGM
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
