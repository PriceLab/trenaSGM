#' @importFrom methods new is
#' @import BiocGenerics
#' @import trena
#' @import MotifDb
#' @import RPostgreSQL

#' @name FootprintDatabaseModelBuilder-class
#' @rdname FootprintDatabaseModelBuilder-class
#' @exportClass FootprintDatabaseModelBuilder

.FootprintDatabaseModelBuilder <- setClass("FootprintDatabaseModelBuilder", contains="ModelBuilder")

#------------------------------------------------------------------------------------------------------------------------
#' Create a FootprintDatabaseModelBuilder object
#'
#' @description
#' tell us what we need to know
#'
#' @rdname FootprintDatabaseModelBuilder-class
#'
#' @param genomeName hg38, mm10, ...
#' @param targetGene in same vocabulary as rownames in the expression matrix
#' @param quiet do or do not print progress information
#'
#' @return An object of the FootprintDatabaseModelBuilder class
#'
#'
#' @examples
#' if(interactive()){
#'
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   fp.specs <- list(title="2000up.200down.fp",
#'                    type="database.footprints",
#'                    chrom="chr6",
#'                    tss=41163186,
#'                    upstream=2000,
#'                    downstream=200,
#'                    db.host="khaleesi.systemsbiology.net",
#'                    databases=list("brain_hint_16", "brain_hint_20", "brain_wellington_16", "brain_wellington_20"),
#'                    motifDiscovery="builtinFimo",
#'                    tfMapping=c("TFClass", "MotifDB"),
#'                    tfPrefilterCorrelation=0.2)
#'   fpc <- FootprintDatabaseModelBuilder("hg38", "TREM2", fp.specs)
#'   }
#'
#' @export
FootprintDatabaseModelBuilder <- function(genomeName, targetGene, strategy, quiet=TRUE)
{
   required.strategy.fields <- c("title", "type", "regions", "tss", "matrix", "db.host", "databases", "motifDiscovery",
                                 "tfMapping", "tfPool", "tfPrefilterCorrelation", "orderModelByColumn", "solverNames")

   for(field in required.strategy.fields)
      if(!field %in% names(strategy))
         stop(sprintf("missing '%s' field in strategy", field))

   obj <- .FootprintDatabaseModelBuilder(ModelBuilder(genomeName=genomeName,
                                                      targetGene=targetGene,
                                                      strategy=strategy,
                                                      quiet=quiet))

   obj

} # FootprintDatabaseModelBuilder
#------------------------------------------------------------------------------------------------------------------------
#' summarize the attributes specifying the creation of a trena gene regulatory model
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class FootprintDatabaseModelBuilder
#'
#' @export
setMethod('show', 'FootprintDatabaseModelBuilder',

    function(object) {
      msg = sprintf("FootprintDatabaseModel object named '%s'", object@strategy$title)
      cat (msg, '\n', sep='')
      })

#------------------------------------------------------------------------------------------------------------------------
#' create regulatory model of the gene, following all the specified options
#'
#' @rdname build
#' @aliases build
#'
#' @param obj An object of class FootprintDatabaseModelBuilder
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
#'                    tfMapping=c("TFClass", "MotifDB"),
#'                    tfPrefilterCorrelation=0.2)
#'   fpBuilder <- FootprintDatabaseModelBuilder("hg38", "TREM2", fp.specs, quiet=TRUE)
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   build(fpBuilder, mtx)
#'   }
#'
#' @export
setMethod('build', 'FootprintDatabaseModelBuilder',

   function (obj) {
      tbl.fp <- .assembleFootprints(obj@strategy, obj@quiet)
      if(obj@strategy$motifDiscovery == "builtinFimo"){
         tbl.fp$motifName <- tbl.fp$name
         mapper <- tolower(obj@strategy$tfMapping)
         stopifnot(all(mapper %in% c("motifdb", "tfclass")))
         tbl.fp <- associateTranscriptionFactors(MotifDb, tbl.fp, source=obj@strategy$tfMapping, expand.rows=TRUE)

         s <- obj@strategy
         xyz <- "FootprintDatabaseModelBuilder, build"
         tbls <- .runTrenaWithRegulatoryRegions(obj@genomeName,
                                                s$tfPool,
                                                obj@targetGene,
                                                tbl.fp,
                                                s$matrix,
                                                s$tfPrefilterCorrelation,
                                                s$solverNames,
                                                obj@quiet)

        } # motifDisocvery, builtinFimo
      tbl.model <- tbls[[1]]
      coi <- s$orderModelByColumn
      if(coi %in% colnames(tbl.model)){
         tbl.model <- tbl.model[order(tbl.model[, coi], decreasing=TRUE),]
         tbls[[1]] <- tbl.model
         }
      return(tbls)
      })

#------------------------------------------------------------------------------------------------------------------------
.assembleFootprints <- function(strategy, quiet=TRUE)
{
   s <- strategy # for lexical brevity

   dbMain <- dbConnect(PostgreSQL(), user="trena", password="trena", host=s$db.host, dbname="hg38")
   all.available <- all(s$databases %in% dbGetQuery(dbMain, "select datname from pg_database")$datname, v=TRUE, ignore.case=TRUE)
   dbDisconnect(dbMain)
   stopifnot(all.available)

   dbConnections <- list()
   fps <- list()

   for(dbName in s$databases){
      if(!quiet) printf("--- opening connection %s", dbName)
      dbConnection <- dbConnect(PostgreSQL(), user="trena", password="trena", host=s$db.host, dbname=dbName)
      if(!quiet) printf("--- querying %s for footprints across %d regions totaling %d bases",
                        dbName, nrow(s$regions), with(s$regions, sum(end-start)))
      tbl.hits <- .multiQueryFootprints(dbConnection, s$regions)
      tbl.hits$chrom <- unlist(lapply(strsplit(tbl.hits$loc, ":"), "[",  1))
      tbl.hits.clean <- tbl.hits # [, c("chrom", "fp_start", "fp_end", "name", "score2", "method")]
      fps[[dbName]] <- tbl.hits.clean
      tbl.hits.clean$database = dbName
      dbDisconnect(dbConnection)
      }

   tbl.fp <- do.call(rbind, fps)
   if(!quiet) printf(" combined tbl.fp: %d %d", nrow(tbl.fp), ncol(tbl.fp))
   tbl.fp$shortMotif <- NA
   missing <- which(!tbl.fp$name %in% names(MotifDb))
   matched <- which(tbl.fp$name %in% names(MotifDb))
   x <- match(tbl.fp$name[matched], names(MotifDb))
   tbl.fp$shortMotif[matched] <- mcols(MotifDb[x])[, "providerName"]
      # our odd convention: MotifDb:associationTranscriptFactors uses BOTH columns, one
      # for the MotifDb mapping, one "shortMotif" for the TFClass mapping
      # TODO (14 may 2018): fix this
   invisible(tbl.fp)

} # .assembleFootprint
#------------------------------------------------------------------------------------------------------------------------
.queryFootprints <- function(db, chrom, start, stop)
{
    query.p0 <- "select loc, chrom, start, endpos from regions"
    query.p1 <- sprintf("where chrom='%s' and start > %d and endpos < %d", chrom, start, stop)

    query.regions <- paste(query.p0, query.p1)
    tbl.regions <- dbGetQuery(db, query.regions)
    if(nrow(tbl.regions) == 0)
       return(data.frame())
    loc.set <- sprintf("('%s')", paste(tbl.regions$loc, collapse="','"))
    query.hits <- sprintf("select * from hits where loc in %s", loc.set)
    invisible(dbGetQuery(db, query.hits))

} # .queryFootprints
#------------------------------------------------------------------------------------------------------------------------
.multiQueryFootprints <- function(db, tbl.regions)
{
   tbl.fpRegions <- data.frame()

   for(r in seq_len(nrow(tbl.regions))){
      chrom <- tbl.regions$chrom[r]
      start <- tbl.regions$start[r]
        end <- tbl.regions$end[r]
      query.p0 <- "select loc, chrom, start, endpos from regions"
      query.p1 <- sprintf("where chrom='%s' and start > %d and endpos < %d", chrom, start, end)
      query.regions <- paste(query.p0, query.p1)
      tbl.fpRegions.new <- dbGetQuery(db, query.regions)
      tbl.fpRegions <- rbind(tbl.fpRegions, tbl.fpRegions.new)
      }

   loc.set <- unique(sprintf("('%s')", paste(tbl.fpRegions$loc, collapse="','")))
   query.hits <- sprintf("select * from hits where loc in %s", loc.set)

   invisible(dbGetQuery(db, query.hits))

} # .multiQueryFootprints
#------------------------------------------------------------------------------------------------------------------------
