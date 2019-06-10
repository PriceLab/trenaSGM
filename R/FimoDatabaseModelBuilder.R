#' @importFrom methods new is
#' @import BiocGenerics
#' @import trena
#' @import MotifDb
#' @import RSQLite

#' @name FimoDatabaseModelBuilder-class
#' @rdname FimoDatabaseModelBuilder-class
#' @exportClass FimoDatabaseModelBuilder

.FimoDatabaseModelBuilder <- setClass("FimoDatabaseModelBuilder",
                                           contains="ModelBuilder",
                                           slots=c(pvalueThreshold="numeric"))

#------------------------------------------------------------------------------------------------------------------------
#' Create a FimoDatabaseModelBuilder object
#'
#' @description
#' tell us what we need to know
#'
#' @rdname FimoDatabaseModelBuilder-class
#'
#' @param genomeName hg38, mm10, ...
#' @param targetGene in same vocabulary as rownames in the expression matrix
#' @param quiet do or do not print progress information
#'
#' @return An object of the FimoDatabaseModelBuilder class
#'
#'
#' @examples
#' if(interactive()){
#'
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   fp.specs <- list(title="2000up.200down.fp",
#'                    type="database.fimos",
#'                    chrom="chr6",
#'                    tss=41163186,
#'                    upstream=2000,
#'                    downstream=200,
#'                    db.host="khaleesi.systemsbiology.net",
#'                    databases=list("brain_hint_16", "brain_hint_20", "brain_wellington_16", "brain_wellington_20"),
#'                    motifDiscovery="builtinFimo",
#'                    tfMapping=c("TFClass", "MotifDB"),
#'                    tfPrefilterCorrelation=0.2)
#'   fpc <- FimoDatabaseModelBuilder("hg38", "TREM2", fp.specs)
#'   }
#'
#' @export
FimoDatabaseModelBuilder <- function(genomeName, targetGene, strategy, quiet=TRUE)
{
   if(!quiet) message("constructing FimoDatabaseModelBuilder")

   required.strategy.fields <- c("title", "type", "regions", "tss", "geneSymbol","matrix",
                                 "db.host", "db.port", "databases", "motifDiscovery","tfMapping", "tfPool",
                                 "tfPrefilterCorrelation", "maxModelSize", "fimoPValueThreshold",
                                 "orderModelByColumn", "solverNames", "annotationDbFile")

   for(field in required.strategy.fields)
      if(!field %in% names(strategy))
         stop(sprintf("missing '%s' field in strategy", field))

   if(!quiet) message(sprintf("constructing FimoDatabaseModelBuilder"))
   obj <- .FimoDatabaseModelBuilder(ModelBuilder(genomeName=genomeName,
                                                 targetGene=targetGene,
                                                 strategy=strategy,
                                                 quiet=quiet))

   obj

} # FimoDatabaseModelBuilder
#------------------------------------------------------------------------------------------------------------------------
#' summarize the attributes specifying the creation of a trena gene regulatory model
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class FimoDatabaseModelBuilder
#'
#' @export
setMethod('show', 'FimoDatabaseModelBuilder',

    function(object) {
      msg = sprintf("FimoDatabaseModel object named '%s'", object@strategy$title)
      cat (msg, '\n', sep='')
      for(database in object@strategy$databases){
         msg = sprintf("database: %s", database)
         cat (msg, '\n', sep='')
         }
      })

#------------------------------------------------------------------------------------------------------------------------
#' create regulatory model of the gene, following all the specified options
#'
#' @rdname build
#' @aliases build
#'
#' @param obj An object of class FimoDatabaseModelBuilder
#' @param strategy a list specifying all the options to build one or more models
#'
#' @return A list with a bunch of tables...
#'
#' @examples
#' if(interactive()){
#'   fp.specs <- list(title="fp.2000up.200down",
#'                    type="database.fimos",
#'                    chrom="chr6",
#'                    tss=41163186,
#'                    start=(tss-2000),
#'                    downstream=(tss+200),
#'                    matrix
#'                    db.host="khaleesi.systemsbiology.net",
#'                    databases=list("brain_hint_20"),
#'                    motifDiscovery="builtinFimo",
#'                    tfMapping=c("TFClass", "MotifDB"),
#'                    tfPrefilterCorrelation=0.2)
#'   fpBuilder <- FimoDatabaseModelBuilder("hg38", "TREM2", fp.specs, quiet=TRUE)
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   build(fpBuilder, mtx)
#'   }
#'
#' @export
#'
setMethod('build', 'FimoDatabaseModelBuilder',

   function (obj) {
      xyz <- "starting FimoDatabaseModelBuilder::build"
      db <- dbConnect(SQLite(), obj@strategy$databases[[1]])
      tbls.hits <- list()
      tbl.regions <- obj@strategy$regions
      pvalThreshold <- 0.001
      for(r in 1:nrow(tbl.regions)){
         query <- sprintf("select * from fimoBindingSites where chrom='%s' and start >= %d and end <= %d and pValue < %f",
                          tbl.regions$chrom[r], tbl.regions$start[r], tbl.regions$end[r], pvalThreshold)
         tbl.tmp <- dbGetQuery(db, query)
         tbls.hits[[r]] <- tbl.tmp
         }
      tbl.hits <- unique(do.call(rbind, tbls.hits))
      tbl.hits <- subset(tbl.hits, pValue <= obj@strategy$fimoPValueThreshold)
      tfs <- mcols(MotifDb[(tbl.hits$motif)])$geneSymbol   # TODO: put this in the database fill (pshannon, 9 jun 2019)
      tbl.hits$tf <- tfs
      signature <- with(tbl.hits, sprintf("%s:%d-%d:%s", chrom, start, end, tf))
      dups <- which(duplicated(signature))
      if(length(dups) > 0)
         tbl.hits <- tbl.hits[-dups,]
      tfs <- unique(tbl.hits$tf)
      recognizedCandidateTFs <- intersect(rownames(obj@strategy$matrix), tfs)
      s <- obj@strategy
      tbls.out <- .runTrenaWithTFsOnly(obj@genomeName,
                                       recognizedCandidateTFs,
                                       obj@targetGene,
                                       s$matrix,
                                       s$tfPrefilterCorrelation,
                                       s$solverNames,
                                       s$annotationDbFile,
                                       obj@quiet)
      tbl.model <- tbls.out$model
      coi <- obj@strategy$orderModelByColumn
      if(coi %in% colnames(tbl.model))
         tbl.model <- tbl.model[order(abs(tbl.model[, coi]), decreasing=TRUE),]
      if(nrow(tbl.model) > s$maxModelSize)
         tbl.model <- head(tbl.model, n=s$maxModelSize)
      tfs.in.model <- tbl.model$gene
      tbl.hits <- subset(tbl.hits, tf %in% tfs.in.model)
      bindingSites <- as.list(table(tbl.hits$tf))
      tbl.model$bindingSites <- as.integer(bindingSites[tbl.model$gene])
      tbls <- list(model=tbl.model, regulatoryRegions=tbl.hits)
      dbDisconnect(db)
      invisible(tbls)
      }) # build

#------------------------------------------------------------------------------------------------------------------------
