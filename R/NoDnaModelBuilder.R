#' @importFrom methods new is
#' @import BiocGenerics
#' @import trena
#' @import MotifDb
#' @import RPostgreSQL

#' @name NoDnaModelBuilder-class
#' @rdname NoDnaModelBuilder-class
#' @exportClass NoDnaModelBuilder

.NoDnaModelBuilder <- setClass("NoDnaModelBuilder", contains="ModelBuilder")

#------------------------------------------------------------------------------------------------------------------------
#' Create a NoDnaModelBuilder object
#'
#' @description
#' tell us what we need to know
#'
#' @rdname NoDnaModelBuilder-class
#'
#' @param genomeName hg38, mm10, ...
#' @param targetGene in same vocabulary as rownames in the expression matrix
#' @param quiet do or do not print progress information
#'
#' @return An object of the NoDnaModelBuilder class
#'
#'
#' @examples
#' if(interactive()){
#'
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   spec <- list(title="2000up.200down.rmm",
#'                type="regions.motifMatching",
#'                tss=41163186,
#'                regions=data.frame(chrom="chr6", start=41162986, end=41165186, stringsAsFactors=FALSE)
#'                matrix=mtx,
#'                motifDiscovery="matchPWM",
#'                matchThreshold=90,
#'                tfMapping=list("MotifDB"),
#'                tfPrefilterCorrelation=0.2)
#'   rmmBuilder <- NoDnaModelBuilder("hg38", "TREM2", spec)
#'   build(rmmBuilder)
#'   }
#'
#' @export
NoDnaModelBuilder <- function(genomeName, targetGene, strategy, quiet=TRUE)
{
   required.strategy.fields <- c("title", "type", "matrix",
                                 "tfPrefilterCorrelation", "tfPool", "orderModelByColumn",
                                 "solverNames")

   for(field in required.strategy.fields)
      if(!field %in% names(strategy))
         stop(sprintf("missing '%s' field in strategy", field))

   obj <- .NoDnaModelBuilder(ModelBuilder(genomeName=genomeName,
                                          targetGene=targetGene,
                                          strategy=strategy,
                                          quiet=quiet))

   obj

} # NoDnaModelBuilder
#------------------------------------------------------------------------------------------------------------------------
#' summarize the attributes specifying the creation of a trena gene regulatory model
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class NoDnaModelBuilder
#'
#' @export
setMethod('show', 'NoDnaModelBuilder',

    function(object) {
      msg = sprintf("NoDnaModelBuilder object named '%s'", object@strategy$title)
      cat (msg, '\n', sep='')
      })

#------------------------------------------------------------------------------------------------------------------------
#' create regulatory model of the gene, following all the specified options
#'
#' @rdname build
#' @aliases build
#'
#' @param obj An object of class NoDnaModelBuilder
#' @param strategy a list specifying all the options to build one or more models
#'
#' @return A list with a bunch of tables...
#'
#' @examples
#' if(interactive()){
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   build.spec <- list(title="trem2.noDNA.allTFs",
#'                      type="noDNA_tfsSupplied",
#'                      matrix=mtx,
#'                      tfPrefilterCorrelation=0.4,
#'                      orderByColumn="rfScore",
#'                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))
#'
#'   builder <- NoDnaModelBuilder("hg38", "TREM2", build.spec, quiet=TRUE)
#'   build(builder)
#'   }
#'
#' @export
setMethod('build', 'NoDnaModelBuilder',

   function (obj) {

      s <- obj@strategy
      tbls <- .runTrenaWithTFsOnly(obj@genomeName,
                                   s$tfPool,
                                   obj@targetGene,
                                   s$tfs,
                                   s$matrix,
                                   s$tfPrefilterCorrelation,
                                   s$solverNames,
                                   obj@quiet)

      tbl.model <- tbls$model
      coi <- obj@strategy$orderModelByColumn
      if(coi %in% colnames(tbl.model))
         tbl.model <- tbl.model[order(tbl.model[, coi], decreasing=TRUE),]
      tbls$model <- tbl.model
      invisible(tbls)
      })

#------------------------------------------------------------------------------------------------------------------------
