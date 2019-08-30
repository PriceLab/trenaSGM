library(RUnit)
library(trenaSGM)
library(org.Hs.eg.db)
library(TrenaProjectBrainCell)
#------------------------------------------------------------------------------------------------------------------------
Sys.setlocale("LC_ALL", "C")
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.new.RData"))
if(!exists("tbl.enhancers"))
  load(system.file(package="trenaSGM", "extdata", "enhancers.TREM2.RData"))
if(!exists("tp"))
   tp <- TrenaProjectBrainCell()
#------------------------------------------------------------------------------------------------------------------------
library(futile.logger)
flog.appender(appender.file("timing-postgresFootprintQueries.log"))

logTimingInfo <- function(msg, timingInfo){
   string <- sprintf("%50s: %6.2f %6.2f %6.2f %6.2f %6.2f", msg,
                     timingInfo[[1]], timingInfo[[2]], timingInfo[[3]], timingInfo[[4]], timingInfo[[5]])
   flog.info(string, "timingLog")
   }
#------------------------------------------------------------------------------------------------------------------------
test_brain <- function()
{
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand

   recipe <- list(title="fp.enhancers",
                      type="footprint.database",
                      geneSymbol="TREM2",
                      regions=tbl.enhancers,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      db.port=5432,
                      databases=c("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                      motifDiscovery="builtinFimo",
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      tfPool=allKnownTFs(),
                      tfMapping=c("TFClass", "MotifDb"),
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   region.count <- nrow(tbl.enhancers)
   region.totalBases <- sum(apply(tbl.enhancers, 1, function(row) 1 + as.numeric(row[["end"]]) - as.numeric(row[["start"]])))
   msg <- sprintf("queryFootprints, %3d regions, %8d bases", region.count, region.totalBases)
   logTimingInfo(msg, system.time(tbl.fp <- trenaSGM:::.queryFootprintsFromDatabase(recipe, FALSE)))

   shoulder <- 10000
   tbl.huge <- data.frame(chrom="chr6", start=tss-shoulder, end=tss+shoulder, stringsAsFactors=FALSE)
   recipe$regions <- tbl.huge
   region.count <- nrow(recipe$regions)
   region.totalBases <- sum(apply(recipe$regions, 1, function(row) 1 + as.numeric(row[["end"]]) - as.numeric(row[["start"]])))
   msg <- sprintf("queryFootprints, %3d regions, %8d bases", region.count, region.totalBases)
   logTimingInfo(msg, system.time(tbl.fp <- trenaSGM:::.queryFootprintsFromDatabase(recipe, FALSE)))

   shoulder <- 100000
   tbl.huge <- data.frame(chrom="chr6", start=tss-shoulder, end=tss+shoulder, stringsAsFactors=FALSE)
   recipe$regions <- tbl.huge
   region.count <- nrow(recipe$regions)
   region.totalBases <- sum(apply(recipe$regions, 1, function(row) 1 + as.numeric(row[["end"]]) - as.numeric(row[["start"]])))
   msg <- sprintf("queryFootprints, %3d regions, %8d bases", region.count, region.totalBases)
   logTimingInfo(msg, system.time(tbl.fp <- trenaSGM:::.queryFootprintsFromDatabase(recipe, FALSE)))


     #----------------------------
     # 6 regions, each 500bp long
     #----------------------------

   starts <- tss-seq(0, 5000, by=1000)
   ends   <- starts + 500
   tbl.regions <- data.frame(chrom=rep("chr6", 6),
                             start=starts,
                             end=ends,
                             stringsAsFactors=FALSE)

   recipe$regions <- tbl.regions
   region.count <- nrow(recipe$regions)
   region.totalBases <- sum(apply(recipe$regions, 1, function(row) 1 + as.numeric(row[["end"]]) - as.numeric(row[["start"]])))
   msg <- sprintf("queryFootprints, %3d regions, %8d bases", region.count, region.totalBases)
   logTimingInfo(msg, system.time(tbl.fp <- trenaSGM:::.queryFootprintsFromDatabase(recipe, FALSE)))

     #----------------------------
     # 20 regions, each 500bp long
     #----------------------------

   starts <- tss-seq(from=0, by=1000, length.out=20)
   ends   <- starts + 500
   tbl.regions <- data.frame(chrom=rep("chr6", 20),
                             start=starts,
                             end=ends,
                             stringsAsFactors=FALSE)

   recipe$regions <- tbl.regions
   region.count <- nrow(recipe$regions)
   region.totalBases <- sum(apply(recipe$regions, 1, function(row) 1 + as.numeric(row[["end"]]) - as.numeric(row[["start"]])))
   msg <- sprintf("queryFootprints, %3d regions, %8d bases", region.count, region.totalBases)
   logTimingInfo(msg, system.time(tbl.fp <- trenaSGM:::.queryFootprintsFromDatabase(recipe, FALSE)))

     #--------------------------------------------------------------------------------------------------
     # 5 regions, each 5000bp long, ignore overlap: they should be equally demanding db lookups anyway
     #--------------------------------------------------------------------------------------------------

   starts <- tss-seq(from=0, by=5000, length.out=5)
   ends   <- starts + 4950
   tbl.regions <- data.frame(chrom=rep("chr6", 5),
                             start=starts,
                             end=ends,
                             stringsAsFactors=FALSE)

   recipe$databases <- c("brain_hint_20")
   recipe$regions <- tbl.regions
   region.count <- nrow(recipe$regions)
   region.totalBases <- sum(apply(recipe$regions, 1, function(row) 1 + as.numeric(row[["end"]]) - as.numeric(row[["start"]])))
   msg <- sprintf("queryFootprints(%15s), %3d regions, %8d bases", recipe$database, region.count, region.totalBases)
   logTimingInfo(msg, system.time(tbl.fp <- trenaSGM:::.queryFootprintsFromDatabase(recipe, FALSE)))

     #----------------------------------------------------------------------------------
     # now use placenta_fp, should take about same time as brain_hint_20
     #----------------------------------------------------------------------------------


   recipe$databases <- c("placenta_fp")
   recipe$regions <- tbl.regions
   region.count <- nrow(recipe$regions)
   region.totalBases <- sum(apply(recipe$regions, 1, function(row) 1 + as.numeric(row[["end"]]) - as.numeric(row[["start"]])))
   msg <- sprintf("queryFootprints(%15s), %3d regions, %8d bases", recipe$database, region.count, region.totalBases)
   logTimingInfo(msg, system.time(tbl.fp <- trenaSGM:::.queryFootprintsFromDatabase(recipe, FALSE)))

} # test_brain
#------------------------------------------------------------------------------------------------------------------------
explore_placenta_fp_speed <- function()
{
   trem2.tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand

   tbl.regions <- tbl.enhancers[1,]

   recipe <- list(title="fp.enhancers",
                      type="footprint.database",
                      geneSymbol="TREM2",
                      regions=tbl.regions,
                      tss=trem2.tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      db.port=5432,
                      databases=c("placenta_fp"),
                      motifDiscovery="builtinFimo",
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      tfPool=allKnownTFs(),
                      tfMapping=c("TFClass", "MotifDb"),
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   region.count <- nrow(tbl.regions)
   region.totalBases <- sum(apply(tbl.regions, 1, function(row) 1 + as.numeric(row[["end"]]) - as.numeric(row[["start"]])))
   recipe$databases <- "placenta_fp"
   msg <- sprintf("queryFootprints(%15s), %3d regions, %8d bases", recipe$database, region.count, region.totalBases)
   print(msg)
   #gc()
   #gcinfo(TRUE)
   #Rprof(filename="Rprof.placenta", interval=0.01)
   logTimingInfo(msg, system.time(tbl.fp <- trenaSGM:::.queryFootprintsFromDatabase(recipe, FALSE)))
   dim(tbl.fp)
   #Rprof(NULL)
   #summaryRprof("Rprof.placenta")

   recipe$databases <- "brain_hint_20"
   msg <- sprintf("queryFootprints(%15s), %3d regions, %8d bases", recipe$database, region.count, region.totalBases)
   print(msg)
   #gc()
   #gcinfo(TRUE)
   #Rprof(filename="Rprof.brain", interval=0.01)
   logTimingInfo(msg, system.time(tbl.fp <- trenaSGM:::.queryFootprintsFromDatabase(recipe, FALSE)))
   dim(tbl.fp)
   #Rprof(NULL)
   #summaryRprof("Rprof.brain")

} # explore_placenta_fp_speed
#------------------------------------------------------------------------------------------------------------------------
