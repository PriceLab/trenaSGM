#' @importFrom methods new is
#' @import BiocGenerics
#' @import trena
#' @import MotifDb
#' @import RPostgreSQL
##  @importFrom FimoClient FimoClientClass
#' @import FimoClient

#' @name RegionsFimoModelBuilder-class
#' @rdname RegionsFimoModelBuilder-class
#' @exportClass RegionsFimoModelBuilder

.RegionsFimoModelBuilder <- setClass("RegionsFimoModelBuilder",
                                     contains="ModelBuilder",
                                     slots=c(fimo="FimoClientClass")
                                     )

#------------------------------------------------------------------------------------------------------------------------
#' Create a RegionsFimoModelBuilder object
#'
#' @description
#' tell us what we need to know
#'
#' @rdname RegionsFimoModelBuilder-class
#'
#' @param genomeName hg38, mm10, ...
#' @param targetGene in same vocabulary as rownames in the expression matrix
#' @param quiet do or do not print progress information
#' @param fimoClient a live connection to a FimoServer
#'
#' @return An object of the RegionsFimoModelBuilder class
#'
#' @examples
#' if(interactive()){
#'
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   spec <- list(title="2000up.200down.rmm",
#'                type="regions.fimo",
#'                tss=41163186,
#'                regions=data.frame(chrom="chr6", start=41162986, end=41165186, stringsAsFactors=FALSE)
#'                matrix=mtx,
#'                motifDiscovery="matchPWM",
#'                matchThreshold=90,
#'                tfMapping=list("MotifDB"),
#'                tfPrefilterCorrelation=0.2)
#'   rmmBuilder <- RegionsFimoModelBuilder("hg38", "TREM2", spec)
#'   build(rmmBuilder)
#'   }
#'
#' @export
RegionsFimoModelBuilder <- function(genomeName, targetGene, strategy, fimoClient, quiet=TRUE)
{
   required.strategy.fields <- c("title", "type", "regions", "tss", "matrix", "motifDiscovery",
                                 "tfMapping", "tfPrefilterCorrelation", "orderModelByColumn", "solverNames")

   for(field in required.strategy.fields)
      if(!field %in% names(strategy))
         stop(sprintf("missing '%s' field in strategy", field))

   obj <- .RegionsFimoModelBuilder(ModelBuilder(genomeName=genomeName,
                                                targetGene=targetGene,
                                                strategy=strategy,
                                                quiet=quiet),
                                   fimo=fimoClient)

   obj

} # RegionsFimoModelBuilder
#------------------------------------------------------------------------------------------------------------------------
#' summarize the attributes specifying the creation of a trena gene regulatory model
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class RegionsFimoModelBuilder
#'
#' @export
setMethod('show', 'RegionsFimoModelBuilder',

    function(object) {
      msg = sprintf("RegionsFimoModelBuilder object named '%s'", object@strategy$title)
      cat (msg, '\n', sep='')
      })

#------------------------------------------------------------------------------------------------------------------------
#' create regulatory model of the gene, following all the specified options
#'
#' @rdname build
#' @aliases build
#'
#' @param obj An object of class RegionsFimoModelBuilder
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
#'   fpBuilder <- RegionsFimoModelBuilder("hg38", "TREM2", fp.specs, quiet=TRUE)
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   build(fpBuilder, mtx)
#'   }
#'
#' @export
setMethod('build', 'RegionsFimoModelBuilder',

   function (obj) {
      mm <- MotifMatcher(obj@genomeName, list())
      tbl.regions <- obj@strategy$regions
      tbl.motifs <- data.frame()
      for(r in seq_len(nrow(tbl.regions))){
        tbl.seq <- getSequence(mm, tbl.regions[r,])
        seq.list <- list(tbl.seq$seq)
        tbl.fimo <- requestMatch(obj@fimo, seq.list)
        tbl.fimo.filtered <- subset(tbl.fimo, p.value <= obj@strategy$matchThreshold)
        if(nrow(tbl.fimo.filtered) > 0){
           browser()
           tbl.fimo.filtered$start <- tbl.fimo.filtered$start + tbl.regions$start[r]
           tbl.fimo.filtered$stop <- tbl.fimo.filtered$stop + tbl.regions$start[r]
           tbl.fimo.filtered$chrom <- tbl.regions$chrom[r]
           tbl.motifs <- rbind(tbl.motifs, tbl.fimo.filtered)
           } # if filtered
        } # for nrow(tbl(regions
      tfs <- mcols(MotifDb[tbl.motifs$motif])$geneSymbol
      tbl.motifs$geneSymbol <- tfs
      tbls <- .runTrenaWithRegulatoryRegions(obj@genomeName,
                                             allKnownTFs(),    # from trenaSGM
                                             obj@targetGene,
                                             tbl.motifs,
                                             obj@strategy$matrix,
                                             obj@strategy$tfPrefilterCorrelation,
                                             obj@strategy$solverNames,
                                             obj@quiet)

      tbl.model <- tbls$model
      coi <- obj@strategy$orderModelByColumn
      if(coi %in% colnames(tbl.model))
         tbl.model <- tbl.model[order(tbl.model[, coi], decreasing=TRUE),]
      tbls$model <- tbl.model
      invisible(tbls)
      })

#------------------------------------------------------------------------------------------------------------------------
