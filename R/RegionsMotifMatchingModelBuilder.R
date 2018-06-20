#' @importFrom methods new is
#' @import BiocGenerics
#' @import trena
#' @import MotifDb
#' @import RPostgreSQL

#' @name RegionsMotifMatchingModelBuilder-class
#' @rdname RegionsMotifMatchingModelBuilder-class
#' @exportClass RegionsMotifMatchingModelBuilder

.RegionsMotifMatchingModelBuilder <- setClass("RegionsMotifMatchingModelBuilder", contains="ModelBuilder")

#------------------------------------------------------------------------------------------------------------------------
#' Create a RegionsMotifMatchingModelBuilder object
#'
#' @description
#' tell us what we need to know
#'
#' @rdname RegionsMotifMatchingModelBuilder-class
#'
#' @param genomeName hg38, mm10, ...
#' @param targetGene in same vocabulary as rownames in the expression matrix
#' @param quiet do or do not print progress information
#'
#' @return An object of the RegionsMotifMatchingModelBuilder class
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
#'   rmmBuilder <- RegionsMotifMatchingModelBuilder("hg38", "TREM2", spec)
#'   build(rmmBuilder)
#'   }
#'
#' @export
RegionsMotifMatchingModelBuilder <- function(genomeName, targetGene, strategy, quiet=TRUE)
{
   required.strategy.fields <- c("title", "type", "regions", "tss", "matrix", "motifDiscovery",
                                 "tfPool", "tfMapping", "tfPrefilterCorrelation",
                                 "orderModelByColumn", "solverNames")

   for(field in required.strategy.fields)
      if(!field %in% names(strategy))
         stop(sprintf("missing '%s' field in strategy", field))

   obj <- .RegionsMotifMatchingModelBuilder(ModelBuilder(genomeName=genomeName,
                                                         targetGene=targetGene,
                                                         strategy=strategy,
                                                         quiet=quiet))

   obj

} # RegionsMotifMatchingModelBuilder
#------------------------------------------------------------------------------------------------------------------------
#' summarize the attributes specifying the creation of a trena gene regulatory model
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class RegionsMotifMatchingModelBuilder
#'
#' @export
setMethod('show', 'RegionsMotifMatchingModelBuilder',

    function(object) {
      msg = sprintf("RegionsMotifMatchingModelBuilder object named '%s'", object@strategy$title)
      cat (msg, '\n', sep='')
      })

#------------------------------------------------------------------------------------------------------------------------
#' create regulatory model of the gene, following all the specified options
#'
#' @rdname build
#' @aliases build
#'
#' @param obj An object of class RegionsMotifMatchingModelBuilder
#' @param strategy a list specifying all the options to build one or more models
#'
#' @return A list with a bunch of tables...
#'
#' @examples
#' if(interactive()){
#'   fp.specs <- list(title="fp.2000up.200down",
#'                    type="database.footprints",
#'                    chrom="chr6",
#'                    tss=41163186,
#'                    start=(tss-2000),
#'                    downstream=(tss+200),
#'                    matrix
#'                    db.host="khaleesi.systemsbiology.net",
#'                    databases=list("brain_hint_20"),
#'                    motifDiscovery="builtinFimo",
#'                    tfMapping=list("TFClass", "MotifDB"),
#'                    tfPrefilterCorrelation=0.2)
#'   fpBuilder <- RegionsMotifMatchingModelBuilder("hg38", "TREM2", fp.specs, quiet=TRUE)
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   build(fpBuilder, mtx)
#'   }
#'
#' @export
setMethod('build', 'RegionsMotifMatchingModelBuilder',

   function (obj) {
      mm <- MotifMatcher(obj@genomeName, as.list(obj@strategy$pfms))
      tbl.regions <- obj@strategy$regions
      tbl.motifs <- findMatchesByChromosomalRegion(mm,
                                                   tbl.regions,
                                                   pwmMatchMinimumAsPercentage=obj@strategy$matchThreshold)
      if(nrow(tbl.motifs) == 0){
         warning(sprintf("failure modeling %s: no motifs found in region for %d pfms",
                          obj@targetGene, length(obj@strategy$pfms)))
         return(list(model=data.frame, regulatoryRegions=data.frame()))
         }
      mappers <- tolower(obj@strategy$tfMapping)
      stopifnot(all(mappers %in% c("motifdb", "tfclass")))

      tbl.motifs.mapped <- associateTranscriptionFactors(MotifDb, tbl.motifs, source=obj@strategy$tfMapping, expand.rows=TRUE)

      s <- obj@strategy
      tbls <- .runTrenaWithRegulatoryRegions(obj@genomeName,
                                             allKnownTFs(),    # from trenaSGM
                                             obj@targetGene,
                                             tbl.motifs.mapped,
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
