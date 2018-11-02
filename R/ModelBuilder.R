.ModelBuilder <- setClass("ModelBuilder",
                       slots = c(
                          genomeName="character",
                          targetGene="character",
                          strategy="list",
                          quiet="logical",
                          state="environment"
                          ),
                       )


#------------------------------------------------------------------------------------------------------------------------
setGeneric('build', signature='obj', function (obj, expression.matrix) standardGeneric ('build'))
#------------------------------------------------------------------------------------------------------------------------
#' Constructor for ModelBuilder
#'
#' @name ModelBuilder
#' @rdname ModelBuilder-class
#'
#' @param genomeName hg38, mm10, ...
#' @param targetGene in same vocabulary as rownames in the expression matrix
#' @param strategy a list of options defining the model to be built
#' @param quiet do or do not print progress information
#'
#' @return An object of the FootprintDatabaseModelBuilder class
#'
#' @return An object of class ModelBuilder
#'
#' @export
ModelBuilder <- function(genomeName, targetGene, strategy, quiet=TRUE)
{
   if(!targetGene %in% rownames(strategy$matrix)){
       msg <- sprintf("your targetGene '%s' is not in the rownames of your expression matrix", targetGene);
       stop(msg)
       }

   obj <- .ModelBuilder(genomeName=genomeName,
                        targetGene=targetGene,
                        strategy=strategy,
                        quiet=quiet,
                        state=new.env(parent=emptyenv()))
   obj

} # ModelBuilder
#----------------------------------------------------------------------------------------------------
# remove NA geneSymbols.  uppercase mouse gene symbols to look like humans
# map to, and replace with, ensembl gene ids.   if multiple - alas, and for now -  just use the first.
# use "NA" rather than NA on output
.replaceGeneSymbolsWithEnsemblGeneIDsInDataFrame <- function(tbl.regRegions, annotationDbFile)
{
   nas <- which(is.na(tbl.regRegions$geneSymbol))
   if(length(nas) > 0)
      tbl.regRegions <- tbl.regRegions[-nas,]

   syms <- tbl.regRegions$geneSymbol
   syms <- unique(toupper(syms))

   printf("----------- select(org.Hs.eg.db, ...) in .replaceGeneSymbolsWithEnsemblGeneIDs")

   annotationDb <- AnnotationDbi::loadDb(annotationDbFile)
   on.exit(RSQLite::dbDisconnect(dbconn(annotationDb)))

   suppressMessages(tbl.map <- AnnotationDbi::select(annotationDb, keys=syms, keytype="SYMBOL", columns="ENSEMBL"))

   dups <- which(duplicated(tbl.map$SYMBOL))
   if(length(dups) > 0)
      tbl.map <- tbl.map[-dups,]
   unmapped <- which(is.na(tbl.map$ENSEMBL))
   if(length(unmapped) > 0)
      tbl.map$ENSEMBL[unmapped] <- tbl.map$SYMBOL[unmapped]
   #stopifnot(all(toupper(tbl.regRegions$geneSymbol) %in%  tbl.map$SYMBOL))

   map.list <- tbl.map$ENSEMBL
   names(map.list) <- tbl.map$SYMBOL

   tbl.regRegions$geneSymbol <- map.list[toupper(tbl.regRegions$geneSymbol)]

   invisible(tbl.regRegions)

} # .replaceGeneSymbolsWithEnsemblGeneIDsInDataFrame
#----------------------------------------------------------------------------------------------------
# remove NA geneSymbols.  uppercase mouse gene symbols to look like humans
# map to, and replace with, ensembl gene ids.   if multiple - alas, and for now -  just use the first.
# use "NA" rather than NA on output
.replaceGeneSymbolsWithEnsemblGeneIDsInList <- function(gene.list, annotationDbFile)
{
   nas <- which(is.na(gene.list))
   if(length(nas) > 0)
      gene.list <- gene.list[-nas,]

   syms <- unique(toupper(gene.list))

   printf("----------- select(org.Hs.eg.db, ...) in .replaceGeneSymbolsWithEnsemblGeneIDsInList")

   annotationDb <- AnnotationDbi::loadDb(annotationDbFile)
   on.exit(RSQLite::dbDisconnect(dbconn(annotationDb)))

   suppressMessages(tbl.map <- AnnotationDbi::select(annotationDb, keys=syms, keytype="SYMBOL", columns="ENSEMBL"))

   dups <- which(duplicated(tbl.map$SYMBOL))
   if(length(dups) > 0)
      tbl.map <- tbl.map[-dups,]
   unmapped <- which(is.na(tbl.map$ENSEMBL))
   if(length(unmapped) > 0)
      tbl.map$ENSEMBL[unmapped] <- tbl.map$SYMBOL[unmapped]

   invisible(unique(tbl.map$ENSEMBL))

} # .replaceGeneSymbolsWithEnsemblGeneIDsInList
#----------------------------------------------------------------------------------------------------
.runTrenaWithRegulatoryRegions <- function(genomeName, allKnownTFs, targetGene, tbl.regulatoryRegions,
                                           expression.matrix, tfPrefilterCorrelation, solverNames,
                                           annotationDbFile, quiet)
{
   trena <- Trena(genomeName, quiet=quiet)

   all.known.tfs.mtx <- intersect(allKnownTFs, rownames(expression.matrix))
   ensembl.tfs <- length(grep("ENSG0", all.known.tfs.mtx)) > 0
   if(ensembl.tfs)
      tbl.regulatoryRegions <- .replaceGeneSymbolsWithEnsemblGeneIDsInDataFrame(tbl.regulatoryRegions, annotationDbFile)

      # for our purposes, a tf can only remain a candidate if it IS a tf, according
      # to GO, and has expression data

   candidate.tfs <- intersect(all.known.tfs.mtx, tbl.regulatoryRegions$geneSymbol)

      # now reduce the expession matrix to just those tfs + targetGene
      # in preparation for calculating correlated expression
      # and further restricting the candidate tfs, and the matrix, by
      # an abs(cor) threshold

   mtx.sub <- expression.matrix[c(candidate.tfs, targetGene),]
   mtx.cor <- cor(t(mtx.sub))
   xyz <- "ModelBuilder before mtx.cor[targetGene,]"
   high.correlation.genes <- names(which(abs(mtx.cor[targetGene,]) >= tfPrefilterCorrelation))
   if(!quiet){
      message(sprintf("%d high correlation tfs: %d +  %d -",
                      length(high.correlation.genes),
                      length(which(mtx.cor[high.correlation.genes, targetGene] > 0)),
                      length(which(mtx.cor[high.correlation.genes, targetGene] < 0))))
      } # !quiet

   if(length(high.correlation.genes) == 1){  # it can only be the targetGene
      msg <- sprintf("NO genes have expression >= %f correlated with targetGene '%s'",
                     tfPrefilterCorrelation, targetGene)
      stop(msg)
      }

      # note that the target gene is in the list, since its correlation is perfect
   mtx.tfs.filtered <- mtx.sub[high.correlation.genes,]

      # we can leave the target gene in this list because the trena solvers
      # always check and exclude it

   tf.candidates.final <- high.correlation.genes
   tbl.regulatoryRegions.filtered <- subset(tbl.regulatoryRegions, geneSymbol %in% tf.candidates.final)

   stopifnot(all(c(targetGene, tbl.regulatoryRegions.filtered$geneSymbol) %in% rownames(mtx.tfs.filtered)))

   tbl.model <- createGeneModelFromRegulatoryRegions(trena, targetGene, solverNames,
                                                     tbl.regulatoryRegions.filtered,
                                                     mtx.tfs.filtered)

   return(list(model=tbl.model, regulatoryRegions=tbl.regulatoryRegions.filtered))

} # .runTrenaWithRegulatoryRegions
#------------------------------------------------------------------------------------------------------------------------
.runTrenaWithTFsOnly <- function(genomeName, tfPool, targetGene, tfList, expression.matrix,
                                 tfPrefilterCorrelation, solverNames, annotationDbFile, quiet)
{
   trena <- Trena(genomeName, quiet=quiet)
   if(!quiet)
      printf("--- entering .runTrenaWithTFsOnly")

   all.known.tfs.mtx <- intersect(tfPool, rownames(expression.matrix))

   mtx.tfs <- expression.matrix[c(all.known.tfs.mtx, targetGene),]
   mtx.cor <- cor(t(mtx.tfs))

   high.correlation.genes <- names(which(abs(mtx.cor[targetGene,]) >= tfPrefilterCorrelation))
   if(!quiet) {
      printf("ModelBuilder::.runTrenaWithTFsOnly: %d tfs in matrix and correlation > %f",
             length(high.correlation.genes), tfPrefilterCorrelation)
      }

   mtx.tfs.filtered <- expression.matrix[high.correlation.genes,]
   tf.candidates.final <- intersect(high.correlation.genes, tfList)

   if(length(tf.candidates.final) <= 1){
      message(sprintf("none of the supplied TFs reach expression correlation threshold(%4.2f) with targetGene %s",
                      tfPrefilterCorrelation, targetGene))
      return(list(model=data.frame(), regulatoryRegions=data.frame()))
      }

   missing.genes <- setdiff(c(targetGene, tf.candidates.final), rownames(mtx.tfs.filtered))
   if(length(missing.genes) > 0){
      msg <- sprintf("these (target + tfs) not in matrix: %s", paste(missing.genes, collapse=", "))
      stop(msg)
      }

   if(!quiet) printf("calling createGeneModelFromTfList")
   tbl.model <- createGeneModelFromTfList(trena, targetGene, solverNames, tf.candidates.final,
                                          mtx.tfs.filtered)

   return(list(model=tbl.model, regulatoryRegions=data.frame()))

} # .runTrenaWithTFsOnly
#------------------------------------------------------------------------------------------------------------------------
# create regulatory model of the gene, following all the specified options
#
# @rdname build
# @aliases build
#
#setMethod('build', 'ModelBuilder',
#
#   function (obj) {
#      printf("entering ModelBuilder::build")
#      x <- callNextMethod(obj)
#      browser()
#      printf("back in ModelBuilder::build")
#      })
#------------------------------------------------------------------------------------------------------------------------
