#' @importFrom methods new is
#' @import BiocGenerics
#' @import trena
#' @import MotifDb
#' @import RPostgreSQL

#' @name FootprintDatabaseModelBuilder-class
#' @rdname FootprintDatabaseModelBuilder-class
#' @exportClass FootprintDatabaseModelBuilder

.FootprintDatabaseModelBuilder <- setClass("FootprintDatabaseModelBuilder",
                                           contains="ModelBuilder",
                                           slots=c(stagedExecutionDirectory="character",
                                                   motifSpeciesRestriction="character"))

#------------------------------------------------------------------------------------------------------------------------
setGeneric('staged.build', signature='obj', function (obj, stage) standardGeneric ('staged.build'))
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
FootprintDatabaseModelBuilder <- function(genomeName, targetGene, strategy, stagedExecutionDirectory=NA_character_, quiet=TRUE)
{
   if(!quiet) message("constructing FootprintDatabaseModelBuilder")

   required.strategy.fields <- c("title", "type", "regions", "tss", "geneSymbol","matrix",
                                 "db.host", "db.port", "databases", "motifDiscovery","tfMapping", "tfPool",
                                 "tfPrefilterCorrelation", "orderModelByColumn", "solverNames",
                                 "annotationDbFile")

   for(field in required.strategy.fields)
      if(!field %in% names(strategy))
         stop(sprintf("missing '%s' field in strategy", field))

   motifSpeciesRestriction <- NA_character_
   if("motifSpeciesRestriction" %in% names(strategy))
      motifSpeciesRestriction <- strategy$motifSpeciesRestriction

   printf("FPDBbuilder ctor, quiet: %s", quiet)
   obj <- .FootprintDatabaseModelBuilder(ModelBuilder(genomeName=genomeName,
                                                      targetGene=targetGene,
                                                      strategy=strategy,
                                                      quiet=quiet),
                                         stagedExecutionDirectory=stagedExecutionDirectory,
                                         motifSpeciesRestriction=motifSpeciesRestriction)

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
#'
setMethod('build', 'FootprintDatabaseModelBuilder',

   function (obj) {
      tbls <- tryCatch({
        if(!obj@quiet) message(sprintf("FootprintDatabaseModleBuilder::build"))
        tbl.fp <- .queryFootprintsFromDatabase(obj@strategy, obj@quiet)
        if(nrow(tbl.fp) == 0){
           stop(base::simpleError("no footprints found"))
           }
        if(!is.na(obj@motifSpeciesRestriction)){  # filter out all but the motifs named
           keepers <- grep(obj@motifSpeciesRestriction, tbl.fp$name, ignore.case=TRUE)
           if(length(keepers) > 0)
              tbl.fp <- tbl.fp[keepers,]
           } # motifSpeciesRestriction
        if(nrow(tbl.fp) == 0){
           stop(base::simpleError("no footprints survived the motif species filter"))
           }
        if(obj@strategy$motifDiscovery == "builtinFimo"){
           if(!obj@quiet) message(sprintf("motifDiscovery: bulitinFimo"))
           tbl.fp$motifName <- tbl.fp$name
           mapper <- tolower(obj@strategy$tfMapping)
           stopifnot(all(mapper %in% c("motifdb", "tfclass", "motifdb+tfclass")))
           #if(mapper == "motifdb+tfclass")
           #   mapper <- c("motifdb", "tfclass")
           if(!obj@quiet){
              message(sprintf("associateTranscriptFactors (with motifs), mappers: %s", paste(mapper, collapse=", ")))
              message(sprintf("  tbl.fp for mapping: %d rows", nrow(tbl.fp)))
              }
           tbl.fp <- associateTranscriptionFactors(MotifDb, tbl.fp, source=obj@strategy$tfMapping, expand.rows=TRUE)
           if(!obj@quiet){
              message(sprintf("after associateTranscriptFactors"))
              message(sprintf("  tbl.fp after mapping: %d rows", nrow(tbl.fp)))
              }
             # an ad hoc processing step: if compound ensembl|geneSymbol identifers are used, which
             # we can find out by checking the target gene and the matrix, then we want to make the candidate
             # tfs into an identifier in the same style
           #ensembl.geneSymbol.names.in.use <- grepl("|", obj@targetGene, fixed=TRUE)
           #if(ensembl.geneSymbol.names.in.use){
           #   ensembl.geneSymbols <- make.ensembl.geneSymbol.identifiers(tbl.fp$geneSymbol)
           #   tbl.fp$geneSymbol <- ensembl.geneSymbols
           #   }
           s <- obj@strategy
           xyz <- "FootprintDatabaseModelBuilder, build"
           tbls <- .runTrenaWithRegulatoryRegions(obj@genomeName,
                                                  s$tfPool,
                                                  obj@targetGene,
                                                  tbl.fp,
                                                  s$matrix,
                                                  s$tfPrefilterCorrelation,
                                                  s$solverNames,
                                                  s$annotationDbFile,
                                                  obj@quiet)

           } # motifDisocvery, builtinFimo
        tbl.model <- tbls[[1]]
        tbl.regRegions <- tbls[[2]]
           # recalculate the bindingSites count, which are likely to be inflated
           # by identical motif/tf/loc matches from multiple samples in the footprints
        tbl.regRegions.unique <- unique(tbl.regRegions[, c("loc", "geneSymbol")])
        tbl.tf.counts <- as.data.frame(table(tbl.regRegions.unique$geneSymbol))
        colnames(tbl.tf.counts) <- c("gene", "bindingSites")
        gene.order <- match(tbl.model$gene, tbl.tf.counts$gene)
        if(!obj@quiet)
           message(sprintf("bindingSite reduction, from total %d to %d across %d tfs in model",
                  sum(tbl.model$bindingSites), sum(tbl.tf.counts$bindingSites), nrow(tbl.model)))
        tbl.model$bindingSites <- tbl.tf.counts$bindingSites[gene.order]
        coi <- s$orderModelByColumn
        if(coi %in% colnames(tbl.model)){
           tbl.model <- tbl.model[order(tbl.model[, coi], decreasing=TRUE),]
           tbls[[1]] <- tbl.model
           }
        tbls
        }, error=function(e){
           message(e)
           return(list(model=data.frame(), regulatoryRegions=data.frame()))
           })

      return(tbls)
      }) # build

#------------------------------------------------------------------------------------------------------------------------
.queryFootprintsFromDatabase <- function(strategy, quiet)
{
   s <- strategy # for lexical brevity

   if(!quiet) {
      message(sprintf("=========================================================================="))
      message(sprintf("FootprintDatabaseModelBuilder::.queryFootprintsFromDatabase, quiet? %s", quiet))
      message(sprintf("==========================================================================="))
      message(sprintf("  s$db.host: %s", s$db.host))
      }

      # dbConnect requires us to specify a database to connect to
      # having done that, we can then query postgres for all of the databases it holds
      # we take these two steps first, then close that initial connection
      # it seems (n=2) that every postgres installation has a database called "postgres":
      # that's what we use for the initial "what databases do you have?" query

   dbMain <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="postgres", host=s$db.host, port=s$db.port)
   availableDatabases <- dbGetQuery(dbMain, "select datname from pg_database")$datname
   requestedDatabases <- unlist(s$databases)

   if(!quiet){
      message(sprintf("==== available databases:"))
      print(availableDatabases)
      message(sprintf("==== requested databases:"))
      print(requestedDatabases)
      }

   all.available <- all(requestedDatabases %in% availableDatabases)
   dbDisconnect(dbMain)
   stopifnot(all.available)

   dbConnections <- list()
   fps <- list()

   for(dbName in s$databases){
      if(!quiet) message(sprintf("--- opening connection %s", dbName))
      dbConnection <- dbConnect(PostgreSQL(), user="trena", password="trena", host=s$db.host, dbname=dbName, port=s$db.port)

      if(!quiet){
         options(digits.secs = 6)
         message(Sys.time())
         message(sprintf("--- querying %s for footprints across %d regions totaling %d bases",
                         dbName, nrow(s$regions), with(s$regions, sum(end-start))))
         }
      tbl.hits <- .multiQueryFootprints(dbConnection, s$regions, quiet)
      if(!quiet){
         message(Sys.time())
         message(sprintf("--- back from multiQueryFootprints, nrow: %d", nrow(tbl.hits)))
         }
      if(nrow(tbl.hits) == 0){
         message(sprintf("--- no footprints found in regions in db '%s'", dbName))
      } else {
         tbl.hits$chrom <- unlist(lapply(strsplit(tbl.hits$loc, ":"), "[",  1))
         tbl.hits.clean <- tbl.hits
         tbl.hits.clean$database = dbName
         fps[[dbName]] <- tbl.hits.clean
         }
      dbDisconnect(dbConnection)
      }

   if(length(fps) == 0){
      message("no footprints found, returning empty data.frame")
      return(data.frame())
      }

   if(!quiet) message(sprintf("about to do.call rbind, length(fps): %d", length(fps)))
   tbl.fp <- do.call(rbind, fps)
   if(!quiet) message(sprintf(" combined tbl.fp: %d %d", nrow(tbl.fp), ncol(tbl.fp)))
   tbl.fp$shortMotif <- NA
   # missing <- which(!tbl.fp$name %in% names(MotifDb))
   matched <- which(tbl.fp$name %in% names(MotifDb))
   x <- match(tbl.fp$name[matched], names(MotifDb))
   tbl.fp$shortMotif[matched] <- mcols(MotifDb[x])[, "providerName"]
      # our odd convention: MotifDb:associationTranscriptFactors uses BOTH columns, one
      # for the MotifDb mapping, one "shortMotif" for the TFClass mapping
      # TODO (14 may 2018): fix this
   invisible(tbl.fp)

} # .queryFootprintsFromDatabase
#------------------------------------------------------------------------------------------------------------------------
#' create regulatory model of the gene, following all the specified options, one stage at a time, saving intermedate data
#'
#' @rdname staged.build
#' @aliases staged.build
#'
#' @param obj An object of class FootprintDatabaseModelBuilder
#' @param strategy a list specifying all the options to build one or more models
#' @param stage a character string, one of "find.footprints", "associateTFs", "build.models"
#'
#' @return A list with two data.frames, "model" and "regulatoryRegions"
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
#'   staged.build(fpBuilder, mtx, stage="find.footprints")
#'   }
#'
#' @export
#'
setMethod('staged.build', 'FootprintDatabaseModelBuilder',

   function (obj, stage=c("find.footprints", "associateTFs", "build.models")) {

      stopifnot(dir.exists(obj@stagedExecutionDirectory))
      stopifnot(stage %in% c("find.footprints", "associateTFs", "build.models"))
      subdirName <- sprintf("%s-%s", obj@targetGene, obj@strategy$geneSymbol)
      subdir.path <- file.path(obj@stagedExecutionDirectory, subdirName)

      if(!file.exists(subdir.path))
         dir.create(subdir.path)
      footprint.filename <-           file.path(subdir.path, "tbl.fp.RData")
      footprints.tfMapped.filename <- file.path(subdir.path, "tbl.fpMapped.RData")
      models.filename <-              file.path(subdir.path, "models.RData")

      tbls <- tryCatch({

        if(stage == "find.footprints"){
           tbl.fp <- .assembleFootprints(obj@strategy, obj@quiet)
           save(tbl.fp, file=footprint.filename)
           if(!obj@quiet)
              printf("saving %d footprints to %s", nrow(tbl.fp), footprint.filename)
           return(footprint.filename)
           } # find.footprints

        if(obj@strategy$motifDiscovery == "builtinFimo" & stage=="associateTFs"){
           stopifnot(file.exists(footprint.filename))
           load(footprint.filename)
           tbl.fp$motifName <- tbl.fp$name
           mapper <- tolower(obj@strategy$tfMapping)
           stopifnot(all(mapper %in% c("motifdb", "tfclass")))
           tbl.fp <- associateTranscriptionFactors(MotifDb, tbl.fp, source=obj@strategy$tfMapping, expand.rows=TRUE)
           save(tbl.fp, file=footprints.tfMapped.filename)
           if(!obj@quiet)
              printf("saving %d tf-mapped footprints to %s", nrow(tbl.fp), footprints.tfMapped.filename)
           return(footprints.tfMapped.filename)
           } # associateTFs

         if(stage == "build.models"){
           load(footprints.tfMapped.filename)
           s <- obj@strategy
           tbls <- .runTrenaWithRegulatoryRegions(obj@genomeName,
                                                  s$tfPool,
                                                  obj@targetGene,
                                                  tbl.fp,
                                                  s$matrix,
                                                  s$tfPrefilterCorrelation,
                                                  s$solverNames,
                                                  s$annotationDbFile,
                                                  obj@quiet)

           tbl.model <- tbls[[1]]
           coi <- s$orderModelByColumn
           if(coi %in% colnames(tbl.model)){
                # note: all columns are ordered by their absolute value.
                #       this could, in principle, cause problems, but current
                #       trena columms (pearson, spearman, beta lasso, beta ridge, random forest)
                 #       are all safely treated this way
              tbl.model <- tbl.model[order(abs(tbl.model[, coi]), decreasing=TRUE),]
              tbls[[1]] <- tbl.model
              }
           save(tbls, file=models.filename)
           if(!obj@quiet)
              printf("saving %d model tfs, %d regulatoryRegions, in %s",
                     nrow(tbls$model), nrow(tbls$regulatoryRegions), models.filename)
           return(models.filename)
           } # build.models
         }, error=function(e){
           print(e)
           if(!obj@quiet)
              printf("saving failed %d model tfs, %d regulatoryRegions, in %s", 0, 0, models.filename)
           tbls <- list(model=data.frame(), regulatoryRegions=data.frame())
           save(tbls, file=models.filename)
           return(models.filename)
            }
         ) # tryCatch

      return(tbls)
      }) # staged.build

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
.multiQueryFootprints <- function(db, tbl.regions, quiet=TRUE)
{
   if(!quiet){
      message(sprintf("--- entering FootprintDatabaseModelBuilder, .multiQueryFootprints, regions: %d", nrow(tbl.regions)))
      message(Sys.time())
      }

   tbls <- list()

   for(r in seq_len(nrow(tbl.regions))){
      chrom <- tbl.regions$chrom[r]
      start <- tbl.regions$start[r]
        end <- tbl.regions$end[r]
      query.p0 <- "select loc, chrom, start, endpos from regions"
      query.p1 <- sprintf("where chrom='%s' and start > %d and endpos < %d", chrom, start, end)
      query.regions <- paste(query.p0, query.p1)
      if(!quiet){
         message(sprintf("--- .multiQueryFootprints, about to query db: %s", query.regions))
         message(Sys.time())
         }

      tbl.fpRegions.new <- dbGetQuery(db, query.regions)
      if(!quiet){
         message(sprintf("--- .multiQueryFootprints, after query.regions, rows returned: %d", nrow(tbl.fpRegions.new)))
         message(Sys.time())
         }

      tbls[[r]] <- tbl.fpRegions.new
      }

   tbl.fpRegions <- do.call(rbind, tbls)
   loc.set <- unique(sprintf("('%s')", paste(tbl.fpRegions$loc, collapse="','")))
   query.hits <- sprintf("select * from hits where loc in %s", loc.set)

   invisible(dbGetQuery(db, query.hits))

} # .multiQueryFootprints
#------------------------------------------------------------------------------------------------------------------------
.assembleFootprints <- function(strategy, quiet)
{
   s <- strategy # for lexical brevity

   if(!quiet) {
      message(sprintf("=========================================================================="))
      message(sprintf("FastFootprintDatabaseModelBuilder::.assembleFootprints, quiet? %s", quiet))
      message(sprintf("==========================================================================="))
      message(sprintf("  s$db.host: %s", s$db.host))
      }

      # dbConnect requires us to specify a database to connect to
      # having done that, we can then query postgres for all of the databases it holds
      # we take these two steps first, then close that initial connection
      # it seems (n=2) that every postgres installation has a database called "postgres":
      # that's what we use for the initial "what databases do you have?" query

   dbMain <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="postgres", host=s$db.host, port=s$db.port)
   availableDatabases <- dbGetQuery(dbMain, "select datname from pg_database")$datname
   requestedDatabases <- unlist(s$databases)

   if(!quiet){
      message(sprintf("==== available databases:"))
      print(availableDatabases)
      message(sprintf("==== requested databases:"))
      print(requestedDatabases)
      }

   all.available <- all(requestedDatabases %in% availableDatabases)
   dbDisconnect(dbMain)
   stopifnot(all.available)

   dbConnections <- list()
   fps <- list()

   for(dbName in s$databases){
      if(!quiet) message(sprintf("--- opening connection %s", dbName))
      dbConnection <- dbConnect(PostgreSQL(), user="trena", password="trena", host=s$db.host, dbname=dbName, port=s$db.port)
      if(!quiet) message(sprintf("--- querying %s for footprints across %d regions totaling %d bases",
                        dbName, nrow(s$regions), with(s$regions, sum(end-start))))
      tbl.hits <- .multiQueryFootprints(dbConnection, s$regions)
      if(nrow(tbl.hits) == 0){
         message(sprintf("--- no footprints found in regions in db '%s'", dbName))
      } else {
         tbl.hits$chrom <- unlist(lapply(strsplit(tbl.hits$loc, ":"), "[",  1))
         tbl.hits.clean <- tbl.hits
         tbl.hits.clean$database = dbName
         fps[[dbName]] <- tbl.hits.clean
         }
      dbDisconnect(dbConnection)
      }

   if(length(fps) == 0){
      message("no footprints found, returning empty data.frame")
      return(data.frame())
      }

   tbl.fp <- do.call(rbind, fps)
   if(!quiet) message(printf(" combined tbl.fp: %d %d", nrow(tbl.fp), ncol(tbl.fp)))
   tbl.fp$shortMotif <- NA
   missing <- which(!tbl.fp$name %in% names(MotifDb))
   matched <- which(tbl.fp$name %in% names(MotifDb))
   x <- match(tbl.fp$name[matched], names(MotifDb))
   tbl.fp$shortMotif[matched] <- mcols(MotifDb[x])[, "providerName"]
      # our odd convention: MotifDb:associationTranscriptFactors uses BOTH columns, one
      # for the MotifDb mapping, one "shortMotif" for the TFClass mapping
      # TODO (14 may 2018): fix this
   invisible(tbl.fp)

} # .assembleFootprints
#------------------------------------------------------------------------------------------------------------------------
