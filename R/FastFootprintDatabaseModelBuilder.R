#' @importFrom methods new is
#' @import BiocGenerics
#' @import trena
#' @import MotifDb
#' @import RPostgreSQL

#' @name FastFootprintDatabaseModelBuilder-class
#' @rdname FastFootprintDatabaseModelBuilder-class
#' @exportClass FastFootprintDatabaseModelBuilder

.FastFootprintDatabaseModelBuilder <- setClass("FastFootprintDatabaseModelBuilder",
                                           contains="ModelBuilder",
                                           slots=c(stagedExecutionDirectory="character"))

#------------------------------------------------------------------------------------------------------------------------
setGeneric('staged.fast.build', signature='obj', function (obj, stage) standardGeneric ('staged.fast.build'))
#------------------------------------------------------------------------------------------------------------------------
#' Create a FastFootprintDatabaseModelBuilder object
#'
#' @description
#' tell us what we need to know
#'
#' @rdname FastFootprintDatabaseModelBuilder-class
#'
#' @param genomeName hg38, mm10, ...
#' @param targetGene in same vocabulary as rownames in the expression matrix
#' @param quiet do or do not print progress information
#'
#' @return An object of the FastFootprintDatabaseModelBuilder class
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
#'   fpc <- FastFootprintDatabaseModelBuilder("hg38", "TREM2", fp.specs)
#'   }
#'
#' @export
FastFootprintDatabaseModelBuilder <- function(genomeName, targetGene, strategy, stagedExecutionDirectory=NA_character_, quiet=TRUE)
{
    required.strategy.fields <- c("title", "type", "regions", "tss", "geneSymbol","matrix",
                                  "db.host", "databases", "motifDiscovery","tfMapping", "tfPool",
                                  "tfPrefilterCorrelation", "orderModelByColumn", "solverNames",
                                  "annotationDbFile")

   for(field in required.strategy.fields)
      if(!field %in% names(strategy))
         stop(sprintf("missing '%s' field in strategy", field))

   obj <- .FastFootprintDatabaseModelBuilder(ModelBuilder(genomeName=genomeName,
                                                      targetGene=targetGene,
                                                      strategy=strategy,
                                                      quiet=quiet),
                                         stagedExecutionDirectory=stagedExecutionDirectory)

   obj

} # FastFootprintDatabaseModelBuilder
#------------------------------------------------------------------------------------------------------------------------
#' summarize the attributes specifying the creation of a trena gene regulatory model
#'
#' @rdname show
#' @aliases show
#'
#' @param obj An object of class FastFootprintDatabaseModelBuilder
#'
#' @export
setMethod('show', 'FastFootprintDatabaseModelBuilder',

    function(object) {
      msg = sprintf("FastFootprintDatabaseModel object named '%s'", object@strategy$title)
      cat (msg, '\n', sep='')
      })

#------------------------------------------------------------------------------------------------------------------------
#' create regulatory model of the gene, following all the specified options
#'
#' @rdname build
#' @aliases build
#'
#' @param obj An object of class FastFootprintDatabaseModelBuilder
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
#'   fpBuilder <- FastFootprintDatabaseModelBuilder("hg38", "TREM2", fp.specs, quiet=TRUE)
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   build(fpBuilder, mtx)
#'   }
#'
#' @export
setMethod('build', 'FastFootprintDatabaseModelBuilder',

   function (obj) {
      tbls <- tryCatch({
        tbl.fp <- .assembleFootprints(obj@strategy, obj@quiet)
        if(obj@strategy$motifDiscovery == "builtinFimo"){
           tbl.fp$motifName <- tbl.fp$name
           mapper <- tolower(obj@strategy$tfMapping)
           stopifnot(all(mapper %in% c("motifdb", "tfclass")))
           tbl.fp <- associateTranscriptionFactors(MotifDb, tbl.fp, source=obj@strategy$tfMapping, expand.rows=TRUE)

           s <- obj@strategy
           xyz <- "FastFootprintDatabaseModelBuilder, build"
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
        coi <- s$orderModelByColumn
        if(coi %in% colnames(tbl.model)){
           tbl.model <- tbl.model[order(tbl.model[, coi], decreasing=TRUE),]
           tbls[[1]] <- tbl.model
           }
        tbls
        }, error=function(e){
           print(e)
           return(list(model=data.frame(), regulatoryRegions=data.frame()))
           })

      return(tbls)
      }) # build

#------------------------------------------------------------------------------------------------------------------------
#' create regulatory model of the gene, following all the specified options, one stage at a time, saving intermedate data
#'
#' @rdname staged.fast.build
#' @aliases staged.fast.build
#'
#' @param obj An object of class FastFootprintDatabaseModelBuilder
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
#'   fpBuilder <- FastFootprintDatabaseModelBuilder("hg38", "TREM2", fp.specs, quiet=TRUE)
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   staged.fast.build(fpBuilder, mtx, stage="find.footprints")
#'   }
#'
#' @export
#'
setMethod('staged.fast.build', 'FastFootprintDatabaseModelBuilder',

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
           tbl.regions <- obj@strategy$regions
           save(tbl.fp, tbl.regions, file=footprint.filename)
           if(!obj@quiet)
              printf("saving %d footprints to %s", nrow(tbl.fp), footprint.filename)
           return(footprint.filename)
           } # find.footprints

        if(obj@strategy$motifDiscovery == "builtinFimo" & stage=="associateTFs"){
           stopifnot(file.exists(footprint.filename))
           load(footprint.filename)
           motifs.to.map <- unique(tbl.fp$name)
           printf("--- mapping %d motifs to genes", length(motifs.to.map))
           print(system.time(tbl.motif2tf <- motifToGene(MotifDb, unique(tbl.fp$name), c("MotifDb", "TFClass"))))
           tf.candidates <- sort(unique(tbl.motif2tf$geneSymbol))
           if(grep("^ENSG", obj@targetGene)){  # rather crude hack for now.  TODO: refactor ensg/geneSymol navigation
              suppressMessages(
                 tf.candidates.ensg <- select(org.Hs.eg.db, keys=tf.candidates, keytype="SYMBOL",
                                              columns=c("SYMBOL", "ENSEMBL"))$ENSEMBL
                 )
              if(!obj@quiet){
                 message(sprintf("target gene is ENSG, converted %d tf geneSymbols from MotifDb to %d tf ensg",
                                 length(tf.candidates), length(tf.candidates.ensg)))
                 }
               tf.candidates <- tf.candidates.ensg
               } # if targetGene has an ensembl (ENSG) identifier
           printf("associateTFs stage found %d tf.candidates", length(tf.candidates))
           save(tbl.fp, tbl.motif2tf, tf.candidates, file=footprints.tfMapped.filename)
           if(!obj@quiet)
              printf("saving %d tf-mapped footprints to %s", nrow(tbl.fp), footprints.tfMapped.filename)
           return(footprints.tfMapped.filename)
           } # associateTFs

         if(stage == "build.models"){
           load(footprints.tfMapped.filename)
           s <- obj@strategy
           tbls <- .runTrenaWithTFsOnly(obj@genomeName,
                                        s$tfPool,
                                        obj@targetGene,
                                        tf.candidates,
                                        s$matrix,
                                        s$tfPrefilterCorrelation,
                                        s$solverNames,
                                        s$annotationDbFile,
                                        obj@quiet)

           # tbls <- .runTrenaWithRegulatoryRegions(obj@genomeName,
           #                                       s$tfPool,
           #                                       obj@targetGene,
           #                                       tf.candidates,
           #                                       s$matrix,
           #                                       s$tfPrefilterCorrelation,
           #                                       s$solverNames,
           #                                       s$annotationDbFile,
           #                                       obj@quiet)

           tbl.model <- tbls[[1]]
           coi <- s$orderModelByColumn
           if(coi %in% colnames(tbl.model)){
              tbl.model <- tbl.model[order(tbl.model[, coi], decreasing=TRUE),]
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
      }) # staged.fast.build

#------------------------------------------------------------------------------------------------------------------------
.assembleFootprints <- function(strategy, quiet)
{
   s <- strategy # for lexical brevity

   printf("=============================================")
   printf(".assembleFootprints, quiet? %s", quiet)
   printf("=============================================")

   dbMain <- dbConnect(PostgreSQL(), user="trena", password="trena", host=s$db.host, dbname="hg38")
   all.available <- all(s$databases %in% dbGetQuery(dbMain, "select datname from pg_database")$datname, v=TRUE, ignore.case=TRUE)
   dbDisconnect(dbMain)
   stopifnot(all.available)

   dbConnections <- list()
   fps <- list()

   for(dbName in s$databases){
      if(!quiet) message(sprintf("--- opening connection %s", dbName))
      dbConnection <- dbConnect(PostgreSQL(), user="trena", password="trena", host=s$db.host, dbname=dbName)
      if(!quiet) message(sprintf("--- querying %s for footprints across %d regions totaling %d bases",
                        dbName, nrow(s$regions), with(s$regions, sum(end-start))))
      tbl.hits <- .multiQueryFootprints(dbConnection, s$regions)
      tbl.hits$chrom <- unlist(lapply(strsplit(tbl.hits$loc, ":"), "[",  1))
      tbl.hits.clean <- tbl.hits # [, c("chrom", "fp_start", "fp_end", "name", "score2", "method")]
      fps[[dbName]] <- tbl.hits.clean
      tbl.hits.clean$database = dbName
      dbDisconnect(dbConnection)
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
