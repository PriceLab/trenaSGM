#' @importFrom methods new is
#' @import BiocGenerics
#' @import trena
#' @import MotifDb
#' @import RPostgreSQL

#' @import motifStack
#'
#' @name trenaSGM-class
#' @rdname trenaSGM-class
#' @exportClass trenaSGM

.trenaSGM <- setClass("trenaSGM",
                       representation = representation(
                          genomeName="character",
                          targetGene="character",
                          mtx="matrix",
                          pfms="List",
                          strategies="list",
                          tfCorrelationThreshold="numeric",
                          quiet="logical",
                          state="environment"
                          ),
                       )


#------------------------------------------------------------------------------------------------------------------------
setGeneric('calculate', signature='obj', function (obj) standardGeneric ('calculate'))
setGeneric('execute.footprint.strategy', signature='obj', function (obj, strategy) standardGeneric ('execute.footprint.strategy'))
#------------------------------------------------------------------------------------------------------------------------
#' Create a trenaSGM object
#'
#' @description
#' tell us what we need to know
#'
#' @rdname trenaSGM-class
#'
#' @param genomeName hg38, mm10, ...
#' @param targetGene in same vocabulary as rownames in the expression matrix
#' @param mtx a numeric expression matrix
#' @param pfms a subset of MotifDb
#' @param strategies strings encoding the various means to be used in building models
#' @param tfCorrelationThreshold absolute value in pearson correlation selects candidate TFs
#' @param quiet do or do not print progress information
#'
#' @return An object of the trenaSGM class
#'
#'
#' @examples
#' if(interactive()){
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   pfms <- query(MotifDb, c("sapiens", "jaspar2018"))
#'   sgm <- trenaSGM("hg38", "TREM2", mtx, pfms,
#'                   strategies=list(footprints="5000.5000.remapped"),
#'                   tfCorrelation=0.2)
#'   }
#'
#' @export
trenaSGM <- function(genomeName, targetGene, mtx, pfms, strategies, tfCorrelationThreshold, quiet=TRUE)
{
   obj <- .trenaSGM(genomeName=genomeName,
                    targetGene=targetGene,
                    mtx=mtx,
                    pfms=pfms,
                    strategies=strategies,
                    tfCorrelationThreshold=tfCorrelationThreshold,
                    quiet=quiet,
                    state=new.env(parent=emptyenv()))
   obj

} # trenaSGM
#------------------------------------------------------------------------------------------------------------------------
#' create regulatory model of the gene, following all the specified options
#'
#' @rdname calculate
#' @aliases calculate
#'
#' @param obj An object of class trenaSGM
#'
#' @return A list with a bunch of tables...
#'
#' @examples
#' if(interactive()){
#'   pfms <- query2(MotifDb, c("sapiens", "jaspar2018"))
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   sgm <- trenaSGM("hg38", "TREM2", mtx, pfms,
#'                   strategies=list(footprints="5000.5000.remapped"),
#'                   tfCorrelation=0.3)
#'   calculate(sgm)
#'   }
#'
#' @export
setMethod('calculate', 'trenaSGM',

   function (obj) {
      footprint.strategies <- grep("footprints", names(strategies))
      for(footprint.strategy in footpint.strategies){
         execute.footprint.strategy(obj, footprint.strategy)
         }
      return(list())
   })

#------------------------------------------------------------------------------------------------------------------------
#' build a model based on the specified footprint strategy
#'
#'
#' @rdname  execute.footprint.strategy
#' @aliases execute.footprint.strategy
#'
#' @param obj An object of class trenaSGM
#' @param footprint.strategy a list of parameters for footprint-driven model building
#'
#'
#' @return A list with two data.frames: a model, and motif regions
#'
#' @examples
#' if(interactive()){
#'   }
#'
#' @export
setMethod('execute.footprint.strategy', 'trenaSGM',

   function (obj, strategy) {
      return(list(model=data.frame(), regions=data.frame()))
      })
#------------------------------------------------------------------------------------------------------------------------
