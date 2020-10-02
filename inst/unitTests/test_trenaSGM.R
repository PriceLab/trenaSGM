library(RUnit)
library(trenaSGM)
library(MotifDb)
library(org.Hs.eg.db)
library(org.At.tair.db)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx")){
   filename <- system.file(package="trenaSGM", "extdata", "mayo.tcx.new.RData")
   if(!file.exists(filename))
      filename <- file.path("../extdata", "mayo.tcx.new.RData")
   stopifnot(file.exists(filename))
   load(filename)
   }


genome <- "hg38"
targetGene <- "TREM2"

if(!exists("sgm"))
   sgm <- trenaSGM(genome, targetGene)

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_trem2_fpdb()
   test_summarizeModels()
   test_allKnownTFs()
   test_trimModel()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   sgm <- trenaSGM("hg38", "TREM2")
   checkEquals(is(sgm), "trenaSGM")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_trem2_fpdb <- function()
{
   printf("--- test_trem2_fpdb")

   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186

      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      regions=tbl.regions,
                      geneSymbol="TREM2",
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      db.port=5432,
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      annotationDbFile=org.Hs.eg.db,
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   build.spec.2 <- build.spec
   build.spec.2$title <- "fp.2000up.200down.02"
   build.spec.2$tfPrefilterCorrelation=0.2
   build.spec.2$orderModelByColumn="pcaMax"
   build.spec.2$tfMapping <- c("MotifDb", "TFClass")

   strategies <- list(one=build.spec, two=build.spec.2)

   models <- calculate(sgm, strategies)
   checkTrue(is.list(models))
   checkEquals(names(models), c("one", "two"))
   checkEquals(names(models$one), c("model", "regulatoryRegions"))
   checkEquals(names(models$two), c("model", "regulatoryRegions"))

      # the two different tfPrefilterCorrelation values produce models with very different row counts
   checkTrue(nrow(models$one$model) < 8)
   checkTrue(nrow(models$two$model) > 20)
   checkTrue(nrow(models$one$regulatoryRegions) < 100)
   checkTrue(nrow(models$two$regulatoryRegions) > 250)

      # check the sort order on the requested orderByColumns

   checkEquals(models$one$model$rfScore, sort(models$one$model$rfScore, decreasing=TRUE))
   checkEquals(models$two$model$pcaMax, sort(models$two$model$pcaMax, decreasing=TRUE))

      # now summarize the two models.
   tbl.summary <- summarizeModels(sgm, orderBy="rfScore", maxTFpredictors=6)
   checkEquals(ncol(tbl.summary), 4)
   checkEquals(colnames(tbl.summary), c("one", "two", "rank.sum", "observed"))
   checkTrue(nrow(tbl.summary) > 5)   # expect 8, but stochasticity may arise

} # test_trem2_fpdb
#------------------------------------------------------------------------------------------------------------------------
test_summarizeModels <- function()
{
   printf("--- test_summarizeModels")

   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186

      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream

   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   spec.1 <- list(title="fp.2000up.200down.04",
                      type="footprint.database",
                      geneSymbol="TREM2",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      db.port=5432,
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      annotationDbFile=org.Hs.eg.db,
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))


   spec.2 <- spec.1
   spec.2$title <- "fp.2000up.200down.02"
   spec.2$tfPrefilterCorrelation=0.2

   spec.3 <- list(title="trem2.rmm.2000up.200down",
                  type="regions.motifMatching",
                  tss=tss,
                  regions=tbl.regions,
                  matrix=mtx,
                  pfms=query(MotifDb, "sapiens", "jaspar2018"),
                  matchThreshold=90,
                  motifDiscovery="matchPWM",
                  tfPool=allKnownTFs(),
                  tfMapping="MotifDB",
                  annotationDbFile=org.Hs.eg.db,
                  quiet=TRUE,
                  tfPrefilterCorrelation=0.4,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   candidate.tfs <- c("IRF5", "IKZF1", "LYL1", "SPI1", "CEBPA", "TFEC",
                      "BHLHE41", "IRF8", "TAL1","ELK3", "POU2F2", "MAFB",
                      "ZBTB18", "bogus")
   spec.4 <- list(title="trem2.noDNA.allTFs",
                  type="noDNA.tfsSupplied",
                  matrix=mtx,
                  geneSymbol="TREM2",
                  candidateTFs=candidate.tfs,
                  tfPool=allKnownTFs(),
                  tfPrefilterCorrelation=0.4,
                  annotationDbFile=org.Hs.eg.db,
                  quiet=TRUE,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   strategies <- list(one=spec.1, two=spec.2, three=spec.3, four=spec.4)

   x <- calculate(sgm, strategies)

   tbl.summary <- summarizeModels(sgm, orderBy="rfScore", maxTFpredictors=8)
     # for a tf to be in the summary, it must appear in at least one model
     # thus there should be a non-NA entry in the first four columns of each row
     # (the fifth and sixth columns are counts, and will always be non-NA: so just
     #  check the first four columns)
   checkTrue(all(as.integer(apply(tbl.summary, 1, function(row) length(which(!is.na(row[1:4]))))) > 0))
      # two tfs not observed in the 0.4 abs(correlation) prefilter:
   checkEquals(ncol(tbl.summary), 6)
   checkEquals(colnames(tbl.summary), c("one", "two", "three", "four", "rank.sum", "observed"))
   checkTrue(nrow(tbl.summary) > 10)   # 17 on (23 may 2018)

} # test_summarizeModels
#------------------------------------------------------------------------------------------------------------------------
test_allKnownTFs <- function()
{
   printf("--- test_allKnownTFs")
   tfs <- allKnownTFs(source="GO:DNAbindingTranscriptionFactorActivity")
   checkTrue(length(tfs) > 1500)
   checkTrue(all(c("AATF", "ADNP", "ADNP2") %in% tfs))

   tfs.ensembl <- allKnownTFs(identifierType="ensemblGeneID")
   checkTrue(length(tfs.ensembl) > length(tfs))  # some identifier expansion
   checkEquals(length(which(is.na(tfs.ensembl))), 0)

   tfs.entrez <- allKnownTFs(identifierType="entrezGeneID")
   checkEquals(length(which(is.na(tfs.entrez))), 0)
   checkEqualsNumeric(length(tfs), length(tfs.entrez), tolerance=5)

   tfs.combined <- allKnownTFs(identifierType="ensembl|geneSymbol")
   checkEquals(head(tfs.combined, n=3), c("ENSG00000275700|AATF", "ENSG00000276072|AATF", "ENSG00000101126|ADNP"))
   checkTrue("ENSG00000175387|SMAD2" %in% tfs.combined)

} # test_allKnownTFs
#------------------------------------------------------------------------------------------------------------------------
test_tair10_frd3 <- function()
{
   load("~/github/trenaShinyApps/mentewab/frd3/shinyApp/data/mtx.zinc.22810x42.RData")
   orf <- "AT3G08040"
   genome <- "tair10"
   sgm <- trenaSGM(genome, orf)

   candidate.tfs <- c("AT5G66940", "AT5G17800", "AT1G51700", "AT1G69570",
                      "AT2G28810", "AT2G01930", "AT3G55370", "AT5G02460",
                      "AT3G53200", "AT3G50410", "AT1G29160", "AT1G55110",
                      "AT1G69560", "AT1G14580", "AT2G02070", "AT5G60850",
                      "AT3G01530", "AT2G02080", "AT3G45610")


   genome <- "tair10"
   targetGene <- orf
   chromosome <- "3"
   tss <- 2569396
      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   build.spec <- list(title="frd3.",
                      type="regions.motifMatching",
                      tss=tss,
                      geneSymbol="FRD3",
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=pfms,
                      matchThreshold=90,
                      motifDiscovery="matchPWM",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.0,
                      annotationDbFile=org.At.tair.db,
                      quiet=TRUE,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)


} # test_tair10_frd3
#------------------------------------------------------------------------------------------------------------------------
test_trimModel <- function()
{
   printf("--- test_trimModel")

   load(system.file(package="trenaSGM", "extdata", "sampleDataForTrimModel.RData"))
   tbl.model.full <- x$model
   tbl.reg.full <- x$regulatoryRegions

      #------------------------------------------------------------
      # first, a very strong model, 3 votes needed
      #------------------------------------------------------------

   tbls.trimmed.3 <- trimModel(tbl.model.full, tbl.reg.full, votesNeeded=3)
   lapply(tbls.trimmed.3, nrow)
   tbl.model <- tbls.trimmed.3$model
   tbl.reg <- tbls.trimmed.3$regulatoryRegions

   checkTrue(nrow(tbl.model) < 7)
   checkTrue(nrow(tbl.reg) > 30)
   checkEquals(sort(tbl.model$gene), sort(unique(tbl.reg$geneSymbol)))

     # validate bindingSites count in tbl.model against UNIQUE sites in tbl.reg
     # first eliminate duplicated sites - which most often (and perhaps exclusively)
     # occur when overlapping footprints from different footprinting runs have
     # different though overlapping coordinates, but identically located motif hits
     # the "loc" column has these motif locations

   duplicated.sites <- which(duplicated(tbl.reg$loc))
   tbl.reg.reduced <- tbl.reg[-duplicated.sites,]
   tbl.siteCounts <- as.data.frame(table(tbl.reg.reduced$geneSymbol))
   colnames(tbl.siteCounts) <- c("gene", "bindingSites")

   list.reg <- tbl.siteCounts$bindingSites
   names(list.reg) <- tbl.siteCounts$gene
   list.reg <- list.reg[sort(names(list.reg))]

   list.model <- tbl.model$bindingSites
   names(list.model) <- tbl.model$gene
   list.model <- list.model[sort(names(list.model))]

   checkEquals(list.reg, list.model)

      #------------------------------------------------------------
      # now, a less stringent model, only 2 votes needed
      #------------------------------------------------------------

   tbls.trimmed.2 <- trimModel(tbl.model.full, tbl.reg.full, votesNeeded=2)
   lapply(tbls.trimmed.2, nrow)
   tbl.model <- tbls.trimmed.2$model
   tbl.reg <- tbls.trimmed.2$regulatoryRegions

   checkTrue(nrow(tbl.model) > 8)
   checkTrue(nrow(tbl.reg) > 60)
   checkEquals(sort(tbl.model$gene), sort(unique(tbl.reg$geneSymbol)))

      #----------------------------------------------------------------
      # now the strict model again, but force keeping EGR4 and ZBTB18
      #----------------------------------------------------------------

   tbls.trimmed.3.plus <- trimModel(tbl.model.full, tbl.reg.full, votesNeeded=3, tf.keepers=c("EGR4", "ZBTB18"))
   lapply(tbls.trimmed.3.plus, nrow)
   tbl.model <- tbls.trimmed.3.plus$model
   tbl.reg <- tbls.trimmed.3.plus$regulatoryRegions

   checkEquals(nrow(tbl.model), 8)
   checkEquals(nrow(tbl.reg), 35)
   checkEquals(sort(tbl.model$gene), sort(unique(tbl.reg$geneSymbol)))

} # test_trimModel
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()

