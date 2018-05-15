#' @importFrom methods new is
#' @import BiocGenerics
#' @import trena
#' @import MotifDb
#' @import RPostgreSQL

#' @import motifStack
#'
#' @name FootprintDatabaseModelBuilder-class
#' @rdname FootprintDatabaseModelBuilder-class
#' @exportClass FootprintDatabaseModelBuilder

.FootprintDatabaseModelBuilder <- setClass("FootprintDatabaseModelBuilder",
                       representation = representation(
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
#'                    tfMapping=list("TFClass", "MotifDB"),
#'                    tfPrefilterCorrelation=0.2)
#'   fpc <- FootprintDatabaseModelBuilder("hg38", "TREM2", fp.specs)
#'   }
#'
#' @export
FootprintDatabaseModelBuilder <- function(genomeName, targetGene, strategy, quiet=TRUE)
{
   obj <- .FootprintDatabaseModelBuilder(genomeName=genomeName,
                                         targetGene=targetGene,
                                         strategy=strategy,
                                         quiet=quiet,
                                         state=new.env(parent=emptyenv()))
   obj

} # FootprintDatabaseModelBuilder
#------------------------------------------------------------------------------------------------------------------------
#' create regulatory model of the gene, following all the specified options
#'
#' @rdname build
#' @aliases build
#'
#' @param obj An object of class FootprintDatabaseModelBuilder
#' @param strategy a list specifying all the options to build one or more models
#' @param expression.matrix a numeric matrix; row names are genes, column names are sample IDs
#'
#' @return A list with a bunch of tables...
#'
#' @examples
#' if(interactive()){
#'   fp.specs <- list(title="fp.2000up.200down",
#'                    type="database.footprints",
#'                    chrom="chr6",
#'                    tss=41163186,
#'                    upstream=2000,
#'                    downstream=200,
#'                    db.host="khaleesi.systemsbiology.net",
#'                    databases=list("brain_hint_20"),
#'                    motifDiscovery="builtinFimo",
#'                    tfMapping=list("TFClass", "MotifDB"),
#'                    tfPrefilterCorrelation=0.2)
#'   fpBuilder <- FootprintDatabaseModelBuilder("hg38", "TREM2", fp.specs, quiet=TRUE)
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   build(fpBuilder, mtx)
#'   }
#'
#' @export
setMethod('build', 'FootprintDatabaseModelBuilder',

   function (obj, expression.matrix) {
      tbl.fp <- .assembleFootprints(obj@strategy, obj@quiet)
      if(obj@strategy$motifDiscovery == "builtinFimo"){
         tbl.fp$motifName <- tbl.fp$name
         tbl.fp <- associateTranscriptionFactors(MotifDb, tbl.fp, source=obj@strategy$tfMapping, expand.rows=TRUE)
         tbls <- .runTrena(obj@genomeName,
                           obj@targetGene,
                            tbl.fp,
                            expression.matrix,
                            obj@strategy$tfPrefilterCorrelation,
                            obj@strategy$solverNames)
        } # motifDisocvery, builtinFimo
      return(tbls)
      })

#------------------------------------------------------------------------------------------------------------------------
.runTrena <- function(genomeName, targetGene, tbl.regulatoryRegions, expression.matrix,
                      tfPrefilterCorrelation, solverNames)
{
   trena <- Trena("hg38", quiet=FALSE)
   all.known.tfs <- unique(mcols(MotifDb)$geneSymbol)
   all.known.tfs.mtx <- intersect(all.known.tfs, rownames(expression.matrix))
   mtx.tfs <- expression.matrix[c(all.known.tfs.mtx, targetGene),]
   mtx.cor <- cor(t(mtx.tfs))
   tf.candidates <- names(which(mtx.cor["TREM2",] >= tfPrefilterCorrelation))
   tf.candidates.with.motifs <- intersect(tf.candidates, tbl.regulatoryRegions$geneSymbol)
   tbl.regulatoryRegions.filtered <- subset(tbl.regulatoryRegions, geneSymbol %in% tf.candidates.with.motifs)
   mtx.tfs.filtered <- mtx.tfs[c(targetGene, tf.candidates.with.motifs),]
   tbl.model <- createGeneModel(trena, targetGene, solverNames, tbl.regulatoryRegions.filtered, mtx.tfs.filtered)
   return(list(model=tbl.model, regulatoryRegions=tbl.regulatoryRegions.filtered))

} # .runTrena
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

   for(dbName in unlist(s$databases)){
      if(!quiet) printf("--- attempting to open %s", dbName)
      dbConnection <- dbConnect(PostgreSQL(), user="trena", password="trena", host=s$db.host, dbname=dbName)
      printf("--- querying %s for footprints in region of %d bases", dbName, 1 + s$end - s$start)
      tbl.hits <- .queryFootprints(dbConnection, s$chrom, s$start, s$end)
      tbl.hits$chrom <- unlist(lapply(strsplit(tbl.hits$loc, ":"), "[",  1))
      tbl.hits.clean <- tbl.hits # [, c("chrom", "fp_start", "fp_end", "name", "score2", "method")]
      fps[[dbName]] <- tbl.hits.clean
      tbl.hits.clean.database = dbName
      dbDisconnect(dbConnection)
      }

   tbl.fp <- do.call(rbind, fps)
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
