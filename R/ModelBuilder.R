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
.runTrena <- function(genomeName, targetGene, tbl.regulatoryRegions, expression.matrix,
                      tfPrefilterCorrelation, solverNames, quiet)
{
   trena <- Trena("hg38", quiet=quiet)
   all.known.tfs <- unique(mcols(MotifDb)$geneSymbol)
   all.known.tfs.mtx <- intersect(all.known.tfs, rownames(expression.matrix))
   mtx.tfs <- expression.matrix[c(all.known.tfs.mtx, targetGene),]
   mtx.cor <- cor(t(mtx.tfs))
   tf.candidates <- names(which(abs(mtx.cor["TREM2",]) >= tfPrefilterCorrelation))
   tf.candidates.with.motifs <- intersect(tf.candidates, tbl.regulatoryRegions$geneSymbol)
   tbl.regulatoryRegions.filtered <- subset(tbl.regulatoryRegions, geneSymbol %in% tf.candidates.with.motifs)
   mtx.tfs.filtered <- mtx.tfs[c(targetGene, tf.candidates.with.motifs),]
   tbl.model <- createGeneModel(trena, targetGene, solverNames, tbl.regulatoryRegions.filtered, mtx.tfs.filtered)
   return(list(model=tbl.model, regulatoryRegions=tbl.regulatoryRegions.filtered))

} # .runTrena
#------------------------------------------------------------------------------------------------------------------------
