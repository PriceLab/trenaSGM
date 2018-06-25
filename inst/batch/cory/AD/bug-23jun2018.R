#------------------------------------------------------------------------------
library(RUnit)
library(trenaSGM)
library(MotifDb)
library(motifStack)
library(dplyr)
library(biomaRt)
library(tibble)
library(BiocParallel)
library(RPostgreSQL)

# function that runs a single trenaSGM
#----------------------------------------------------------------------------------------------------
runSGM <- function(specs)
{
   targetGene <- specs$targetGene

   tbl.tssLoc <- subset(tbl.tss, ensembl_gene_id == targetGene)
   tss <- tbl.tssLoc$transcription_start_site

   tbl.regions <- switch(specs$regionsMode,
                         "enhancers" = {
                            enhancer.list[[targetGene]][, c("chrom", "start", "end")]
                            },
                         "tiny" = {
                            with(tbl.tssLoc, data.frame(chrom=chromosome_name, start=tss-1000,
                                                        end=tss+1000, stringsAsFactors=FALSE))
                         })

  # define the genome to be used
  genome <- "hg38"

  # construct the sgm
  sgm <- trenaSGM(genome, targetGene, quiet=FALSE)

  build.spec <- list(title=targetGene,
                     type="footprint.database",
                     regions=tbl.regions,
                     tss=tss,
                     matrix=mtx,
                     db.host="khaleesi.systemsbiology.net",
                     databases=specs$db,
                     motifDiscovery="builtinFimo",
                     tfMapping=c("MotifDB","TFClass"),
                     tfPrefilterCorrelation=specs$correlationThreshold,
                     tfPool=allKnownTFs(identifierType="ensemblGeneID"),
                     orderModelByColumn="rfScore",
                     solverNames=spec$solvers
                     )

  strategies <- list(one=build.spec)
  model.spear <- calculate(sgm, strategies)
  #model.spear$one$model
  return(model.spear)

} #runSGM
#----------------------------------------------------------------------------------------------------
test_runSGM <- function()
{
   ad.genes <- c("TREM2", "RIN3", "MEF2C")
   targetGenes <- tbl.tss[match(ad.genes, tbl.tss$hgnc_symbol), "ensembl_gene_id"]
   names(targetGenes) <- ad.genes

   spec.trem2 <- list(targetGene=targetGenes[["TREM2"]],
                      regionsMode="tiny",
                      correlationThreshold=0.7,
                      solvers= c("pearson", "spearman"),
                      dbs="brain_hint_20"
                      )
   spec.rin3 <- spec.trem2
   spec.rin3$targetGene <- targetGenes[["RIN3"]]

   spec.mef2c <- spec.trem2
   spec.mef2c$targetGene <- targetGenes[["MEF2C"]]

   x.trem2 <- runSGM(spec.trem2)
   dim(x.trem2$one$model)
   checkTrue(all(x.trem2$one$regulatoryRegions$geneSymbol %in% x.trem2$one$model$gene))

   x.rin3 <- runSGM(spec.rin3)
   dim(x.rin3$one$model)
   checkTrue(all(x.rin3$one$regulatoryRegions$geneSymbol %in% x.rin3$one$model$gene))

   x.mef2c <- runSGM(spec.mef2c)
   dim(x.mef2c$one$model)
   checkTrue(all(x.mef2c$one$regulatoryRegions$geneSymbol %in% x.mef2c$one$model$gene))

   system.time(x <- lapply(list(trem2=spec.trem2, rin3=spec.rin3, mef2c=spec.mef2c), runSGM))
   system.time(y <- bplapply(list(trem2=spec.trem2, rin3=spec.rin3, mef2c=spec.mef2c), runSGM))

} # test_runSGM
#----------------------------------------------------------------------------------------------------
getTSS.ensembl <- function(ensemblIDs){

  browser()
  # Switch the name of the database and filter we use
  db.name <-  "hsapiens_gene_ensembl"

  filter.name <- "ensembl_gene_id"

  my.mart <- biomaRt::useMart(biomart="ensembl", dataset= db.name)

  # ensemblIDs <- head(ensemblIDs)

  tbl.geneInfo <- biomaRt::getBM(attributes=c("chromosome_name",
                                              "transcription_start_site",
                                              "transcript_tsl",
                                              "hgnc_symbol",
                                              filter.name),
                                 filters=filter.name, value=ensemblIDs, mart=my.mart)

  if(nrow(tbl.geneInfo) == 0)
    return(NA)

  # Sort by hgnc_symbol and transcript_tsl, then pull the first entry for each gene
  tbl.geneInfo <- tbl.geneInfo[order(tbl.geneInfo[[filter.name]],
                                     tbl.geneInfo$transcript_tsl),]
  tbl.geneInfo <- tbl.geneInfo[match(unique(tbl.geneInfo[[filter.name]]),
                                     tbl.geneInfo[[filter.name]]),]

  # remove contigs and check to make sure it's just 1 chromosome
  tbl.geneInfo <- subset(tbl.geneInfo, chromosome_name %in% c(1:22, "X", "Y", "MT"))
  tbl.geneInfo$chromosome_name <- sprintf("chr%s", tbl.geneInfo$chromosome_name)

  return (tbl.geneInfo)

} # getTSS.ensembl
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
