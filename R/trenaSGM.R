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
                          quiet="logical",
                          state="environment"
                          ),
                       )


#------------------------------------------------------------------------------------------------------------------------
setGeneric('calculate', signature='obj', function (obj, strategies) standardGeneric ('calculate'))
setGeneric('summarizeModels', signature='obj', function (obj, orderBy, maxTFpredictors) standardGeneric ('summarizeModels'))
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
trenaSGM <- function(genomeName, targetGene, quiet=TRUE)
{
   obj <- .trenaSGM(genomeName=genomeName,
                    targetGene=targetGene,
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
#' @param strategies a list of lists, specifying all the options to build one or more models
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

   function (obj, strategies) {
      obj@state$strategies <- strategies
      models <- list()
      for(name in names(strategies)){
         strategy <- strategies[[name]]
         if(strategy$type == "footprint.database"){
            builder <- FootprintDatabaseModelBuilder(obj@genomeName, obj@targetGene, strategy, quiet=obj@quiet)
            models[[name]] <- build(builder)
            } #
         } # for strategy
      obj@state$models <- models
      return(models)
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
      browser()
      return(list(model=data.frame(), regions=data.frame()))
      })
#------------------------------------------------------------------------------------------------------------------------
#' summarize the gene regulatory models preciously calculated by this trensSGM instance
#'
#' @rdname summarizeModels
#' @aliases summarizeModels
#'
#' @param obj An object of class trenaSGM
#' @param orderBy a characters string, the name of the column in the standard trena model data.frame, typically pcaMax, rfScore, pearsonCoeff
#' @param maxTFpredictors an integer, the number of tfs to extract from each model (when present)
#'
#' @return A data.frame summarizing tfs by model type, by count (the number of models in which it appears) and rank in each model
#'
#' @examples
#' if(interactive()){
#'   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#'   genome <- "hg38"
#'   targetGene <- "TREM2"
#'   sgm <- trenaSGM(genome, targetGene)
#'   chromosome <- "chr6"
#'   upstream <- 2000
#'   downstream <- 200
#'   tss <- 41163186
#'
#'      # strand-aware start and end: trem2 is on the minus strand
#'   start <- tss - downstream
#'   end   <- tss + upstream
#'
#'   build.spec <- list(title="fp.2000up.200down.04",
#'                      type="footprint.database",
#'                      chrom=chromosome,
#'                      start=start,
#'                      end=end,
#'                      tss=tss,
#'                      matrix=mtx,
#'                      db.host="khaleesi.systemsbiology.net",
#'                      databases=list("brain_hint_20"),
#'                      motifDiscovery="builtinFimo",
#'                      tfMapping="MotifDB",
#'                      tfPrefilterCorrelation=0.4,
#'                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))
#'
#'
#'   build.spec.2 <- build.spec
#'   build.spec.2$title <- "fp.2000up.200down.02"
#'   build.spec.2$tfPrefilterCorrelation=0.2
#'
#'   strategies <- list(one=build.spec, two=build.spec.2)
#'   models <- calculate(sgm)
#'   tbl.summary <- summarizeModels(sgm, orderBy="rfScore", maxTFpredictors=5)
#'   }
#'
#' @export
setMethod('summarizeModels', 'trenaSGM',

   function(obj, orderBy, maxTFpredictors) {

      models <- obj@state$models
         #  collect the top tfs from each model
      top.tfs <- c()
      for(model in models){
        tbl.model <- model$model
        stopifnot(orderBy %in% colnames(tbl.model))
        tbl.model <- tbl.model[order(tbl.model[, orderBy], decreasing=TRUE),]
        tfs.this.model <- tbl.model$gene
        if(nrow(tbl.model) > maxTFpredictors)
           tfs.this.model <- tfs.this.model[1:maxTFpredictors]
        top.tfs <- c(top.tfs, tfs.this.model)
        } # for model
      top.tfs <- sort(unique(top.tfs))

         #------------------------------------------------------------
         # build up an empty table
         #------------------------------------------------------------

      strategies <- obj@state$strategies
      model.names <- unlist(lapply(strategies, function(strategy) strategy$title), use.names=FALSE)
      count <- length(top.tfs)
      empty.column <- rep(0, count)
      tbl <- data.frame(bogus.column=empty.column)

      for(i in seq_len(length(model.names))){
         tbl <- cbind(tbl, empty.column, stringsAsFactors=FALSE)
         } # for i
      bogus.column.position <- grep("bogus.column", colnames(tbl))
      tbl <- tbl[, -bogus.column.position]
      model.columns <- grep("empty.column", colnames(tbl))
      colnames(tbl)[model.columns] <- model.names
      rownames(tbl) <- top.tfs

         #------------------------------------------------------------
         # now populate the table
         #------------------------------------------------------------

      failed.match <- NA
      stopifnot(length(model.names) == length(models))
      for(i in seq_len(length(model.names))){
         model.name <- model.names[i]
         tbl.model <- models[[i]]$model
         tbl[top.tfs, model.name] <- match(top.tfs, tbl.model$gene, nomatch=failed.match)
         } # for model.name

      rank.sum <- apply(tbl, 1, function(row) sum(row, na.rm=TRUE))
      observed  <- apply(tbl, 1, function(row) length(which(!is.na(row))))
      tbl$rank.sum <- rank.sum
      tbl$observed <- observed
      tbl <- tbl[order(tbl$rank.sum, decreasing=FALSE),]
      tbl
      })

#------------------------------------------------------------------------------------------------------------------------

  # max <- 10
  # top.tfs <- sort(unique(c(tbl.model.fpReduced.mdb$gene[1:max],
  #                          tbl.model.fpReduced.tfc$gene[1:max],
  #                          tbl.model.enhancerAll.mdb$gene[1:max],
  #                          tbl.model.enhancerAll.tfc$gene[1:max],
  #                          tbl.model.dhsAll.tfc$gene[1:max],
  #                          tbl.model.dhsAll.tfc$gene[1:max])))
  # count <- length(top.tfs)
  # tbl <- data.frame(fp2200rank.mdb=rep(0, count),
  #                   fp2200rank.tfc=rep(0, count),
  #                   enhancers.mdb=rep(0, count),
  #                   enhancers.tfc=rep(0, count),
  #                   dhs.mdb=rep(0, count),
  #                   dhs.tfc=rep(0, count),
  #                   stringsAsFactors=FALSE)
  # rownames(tbl) <- top.tfs
  #
  # failed.match <- NA
  # tbl[top.tfs, "fp2200rank.mdb"] <- match(top.tfs, tbl.model.fpReduced.mdb$gene, nomatch=failed.match)
  # tbl[top.tfs, "fp2200rank.tfc"] <- match(top.tfs, tbl.model.fpReduced.tfc$gene, nomatch=failed.match)
  # tbl[top.tfs, "enhAncers.mdb"]  <- match(top.tfs, tbl.model.enhancerAll.mdb$gene, nomatch=failed.match)
  # tbl[top.tfs, "enhancers.tfc"]  <- match(top.tfs, tbl.model.enhancerAll.tfc$gene, nomatch=failed.match)
  # tbl[top.tfs, "dhs.mdb"]  <- match(top.tfs, tbl.model.dhsAll.mdb$gene, nomatch=failed.match)
  # tbl[top.tfs, "dhs.tfc"]  <- match(top.tfs, tbl.model.dhsAll.tfc$gene, nomatch=failed.match)
  #
  # tbl$rank.sum <- apply(tbl, 1, function(row) sum(row, na.rm=TRUE))
  # tbl$observed  <- apply(tbl, 1, function(row) length(which(!is.na(row))))
  # tbl <- subset(tbl, observed >= 3)
  # tbl <- tbl[order(tbl$rank.sum, decreasing=FALSE),]
  #
