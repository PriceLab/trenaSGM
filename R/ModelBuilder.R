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
   obj <- .ModelBuilder(genomeName=genomeName,
                        targetGene=targetGene,
                        strategy=strategy,
                        quiet=quiet,
                        state=new.env(parent=emptyenv()))
   obj

} # ModelBuilder
#----------------------------------------------------------------------------------------------------
.runTrenaWithRegulatoryRegions <- function(genomeName, allKnownTFs, targetGene, tbl.regulatoryRegions,
                                           expression.matrix, tfPrefilterCorrelation, solverNames, quiet)
{
   trena <- Trena(genomeName, quiet=quiet)

   all.known.tfs.mtx <- intersect(allKnownTFs, rownames(expression.matrix))

      # for our purposes, a tf can only remain a candidate if it IS a tf, according
      # to GO, and has expression data
   candidate.tfs <- intersect(all.known.tfs.mtx, tbl.regulatoryRegions$geneSymbol)

      # now reduce the expession matrix to just those tfs + targetGene
      # in preparation for calculating correlated expression
      # and further restricting the candidate tfs, and the matrix, by
      # an abs(cor) threshold

   mtx.sub <- expression.matrix[c(all.known.tfs.mtx, targetGene),]
   mtx.cor <- cor(t(mtx.sub))
   high.correlation.genes <- names(which(abs(mtx.cor[targetGene,]) >= tfPrefilterCorrelation))

      # note that the target gene is in the list, since its correlation is perfect
   mtx.tfs.filtered <- mtx.sub[high.correlation.genes,]

      # we can leave the target gene in this list because the trena solvers
      # always check and exclude it

   tf.candidates.final <- high.correlation.genes
   tbl.regulatoryRegions.filtered <- subset(tbl.regulatoryRegions, geneSymbol %in% tf.candidates.final)

   stopifnot(all(c(targetGene, tbl.regulatoryRegions.filtered$geneSymbol) %in% rownames(mtx.tfs.filtered)))

   tbl.model <- createGeneModelFromRegulatoryRegions(trena, targetGene, solverNames,
                                                     tbl.regulatoryRegions.filtered, mtx.tfs.filtered)

   return(list(model=tbl.model, regulatoryRegions=tbl.regulatoryRegions.filtered))

} # .runTrenaWithRegulatoryRegions
#------------------------------------------------------------------------------------------------------------------------
.runTrenaWithTFsOnly <- function(genomeName, tfPool, targetGene, tfList, expression.matrix,
                                 tfPrefilterCorrelation, solverNames, quiet)
{
   trena <- Trena(genomeName, quiet=quiet)

   all.known.tfs.mtx <- intersect(tfPool, rownames(expression.matrix))
   candidate.tfs <- intersect(all.known.tfs.mtx, tfList)

   mtx.tfs <- expression.matrix[c(candidate.tfs, targetGene),]
   mtx.cor <- cor(t(mtx.tfs))

   high.correlation.genes <- names(which(abs(mtx.cor[targetGene,]) >= tfPrefilterCorrelation))

   mtx.tfs.filtered <- expression.matrix[high.correlation.genes,]
   tf.candidates.final <- intersect(high.correlation.genes, candidate.tfs)

   if(length(tf.candidates.final) <= 1){
      message(sprintf("none of the supplied TFs reach expression correlation threshold(%4.2f) with targetGene %s",
                      tfPrefilterCorrelation, targetGene))
      return(list(model=data.frame(), regulatoryRegions=data.frame()))
      }

   stopifnot(all(c(targetGene, tf.candidates.final) %in% rownames(mtx.tfs.filtered)))

   tbl.model <- createGeneModelFromTfList(trena, targetGene, solverNames, tf.candidates.final, mtx.tfs.filtered)

   return(list(model=tbl.model, regulatoryRegions=data.frame()))

} # .runTrenaWithTFsOnly
#------------------------------------------------------------------------------------------------------------------------
