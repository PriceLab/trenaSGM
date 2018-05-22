.ModelBuilder <- setClass("ModelBuilder",
                       slots = c(
                          genomeName="character",
                          targetGene="character",
                          strategy="list",
                          quiet="logical",
                          allKnownTFs="character",
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
   load(system.file(package="trenaSGM", "extdata", "tfCollections",
                    "GO_00037000_DNAbindingTranscriptionFactorActivityHuman.RData"))

   obj <- .ModelBuilder(genomeName=genomeName,
                        targetGene=targetGene,
                        strategy=strategy,
                        quiet=quiet,
                        allKnownTFs=tfs,
                        state=new.env(parent=emptyenv()))
   obj

} # ModelBuilder
#----------------------------------------------------------------------------------------------------
.runTrenaWithRegulatoryRegions <- function(genomeName, allKnownTFs, targetGene, tbl.regulatoryRegions, expression.matrix,
                                           tfPrefilterCorrelation, solverNames, quiet)
{
   trena <- Trena(genomeName, quiet=quiet)

   all.known.tfs.mtx <- intersect(allKnownTFs, rownames(expression.matrix))
   mtx.tfs <- expression.matrix[c(all.known.tfs.mtx, targetGene),]
   mtx.cor <- cor(t(mtx.tfs))
   tf.candidates <- names(which(abs(mtx.cor[targetGene,]) >= tfPrefilterCorrelation))

   tf.candidates.with.motifs <- intersect(tf.candidates, tbl.regulatoryRegions$geneSymbol)
   tbl.regulatoryRegions.filtered <- subset(tbl.regulatoryRegions, geneSymbol %in% tf.candidates.with.motifs)
   mtx.tfs.filtered <- mtx.tfs[c(targetGene, tf.candidates.with.motifs),]
   tbl.model <- createGeneModelFromRegulatoryRegions(trena, targetGene, solverNames,
                                                      tbl.regulatoryRegions.filtered, mtx.tfs.filtered)

   return(list(model=tbl.model, regulatoryRegions=tbl.regulatoryRegions.filtered))

} # .runTrenaWithRegulatoryRegions
#------------------------------------------------------------------------------------------------------------------------
.runTrenaWithTFsOnly <- function(genomeName, allKnownTFs, targetGene, tfList, expression.matrix,
                                 tfPrefilterCorrelation, solverNames, quiet)
{
   trena <- Trena(genomeName, quiet=quiet)
   all.known.tfs.mtx <- intersect(allKnownTFs, rownames(expression.matrix))
   candidate.tfs <- intersect(all.known.tfs.mtx, tfList)
   mtx.tfs <- expression.matrix[c(candidate.tfs, targetGene),]
   mtx.cor <- cor(t(mtx.tfs))
   tf.candidates <- names(which(abs(mtx.cor[targetGene,]) >= tfPrefilterCorrelation))
   mtx.tfs.filtered <- expression.matrix[tf.candidates,]
   tbl.model <- createGeneModelFromTfList(trena, targetGene, solverNames, tf.candidates, mtx.tfs.filtered)

   return(list(model=tbl.model, regulatoryRegions=data.frame()))

} # .runTrenaWithTFsOnly
#------------------------------------------------------------------------------------------------------------------------
