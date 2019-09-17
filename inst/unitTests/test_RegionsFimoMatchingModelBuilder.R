library(RUnit)
library(trenaSGM)
library(motifStack)
library(FimoClient)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.new.RData"))


# fimo setup for these tests
# cd ~/github/fimoService/server/
# make -f makefile.unitTests
#    PORT=5560
#    FIMO="/Users/paul/meme/bin/fimo"
#    MOTIFS="../pfms/trem2.curated.15.pfms.meme"
#    start:
#      	$(PYTHON) -i $(PYTHON_SCRIPT) $(PORT) $(FIMO) $(MOTIFS)
#    which resolves to  /Users/paul/anaconda3/bin/python -i /Users/paul/github/fimoService/server/runFimoServer.py 5560 "/Users/paul/meme/bin/fimo" "../pfms/trem2.curated.15.pfms.meme"

if(!exists("fimo")){
   FIMO_HOST <- "localhost"
   FIMO_PORT <- 5560
   fimo <-FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   }

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_build.trem2.400bp.motifDB.model()
   test_build.trem2.twoRegions.motifDB.model()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186

      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)

   build.spec <- list(title="trem2.rmm.2000up.200down",
                      type="regions.fimo",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=list(),
                      motifDiscovery="fimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      quiet=TRUE,
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      fimoThreshold=1e-4,
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- RegionsFimoModelBuilder(genome, targetGene,  build.spec, fimo, quiet=TRUE)
   checkTrue(all(c("RegionsFimoModelBuilder", "ModelBuilder") %in%  is(builder)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.400bp.motifDB.model <- function()
{
   printf("--- test_build.trem2.400bp.motifDb.model")

   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+200, stringsAsFactors=FALSE)

   build.spec <- list(title="trem2.rmm.2000up.200down",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, "sapiens", "jaspar2018"),
                      matchThreshold=1e-4,
                      motifDiscovery="fimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      quiet=TRUE,
                      fimoThreshold=1e-4,
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "ridge", "spearman", "pearson", "randomForest", "xgboost"))

   # builder <- RegionsMotifMatchingModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   builder <- RegionsFimoModelBuilder(genome, targetGene,  build.spec, fimo, quiet=TRUE)

   x <- build(builder)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regRegions <- x$regulatoryRegions
   tbl.model <- x$model
   some.expected.tfs <- c("IKZF1", "CEBPA", "TAL1")
   checkTrue(all(some.expected.tfs %in% tbl.model$gene))
   checkTrue(all(some.expected.tfs %in% tbl.regRegions$geneSymbol))
   checkTrue(all(tbl.regRegions$chrom == chromosome))
   checkTrue(all(tbl.regRegions$fmotifStart >= min(tbl.regions$start)))
   checkTrue(all(tbl.regRegions$motifEnd <= max(tbl.regions$end)))

} # test_build.trem2.400bp.motifDB.model
#------------------------------------------------------------------------------------------------------------------------
test_build.trem2.twoRegions.motifDB.model <- function()
{
   printf("--- test_build.trem2.twoRegions.motifDb.model")

   genome <- "hg38"
   targetGene <- "TREM2"
   chromosome <- "chr6"
   tss <- 41163186
      # strand-aware start and end: trem2 is on the minus strand
   tbl.regions <- data.frame(chrom=rep(chromosome, 2),
                             start=c(tss-5000, tss),
                             end=c(tss, tss+5000),
                             stringsAsFactors=FALSE)

   build.spec <- list(title="trem2.fimo.5kup.5kdown",
                      type="regions.motifMatching",
                      tss=tss,
                      regions=tbl.regions,
                      matrix=mtx,
                      pfms=query(MotifDb, "sapiens", "jaspar2018"),
                      matchThreshold=1e-4,
                      motifDiscovery="fimo",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      fimoThreshold=1e-4,
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- RegionsFimoModelBuilder(genome, targetGene,  build.spec, fimo, quiet=TRUE)
   x <- build(builder)
   checkEquals(names(x), c("model", "regulatoryRegions"))
   tbl.regRegions <- x$regulatoryRegions
   tbl.model <- x$model
   some.expected.tfs <- c("IRF5", "IKZF1", "SPI1", "TFEC", "CEBPA", "IRF8")
   checkTrue(all(some.expected.tfs %in% tbl.regRegions$geneSymbol))
   checkTrue(all(tbl.regRegions$chrom == chromosome))
   checkTrue(all(tbl.regRegions$fmotifStart >= min(tbl.regions$start)))
   checkTrue(all(tbl.regRegions$motifEnd <= max(tbl.regions$end)))

   # now make sure that there are motif hits in BOTH of the two, non-overlapping
   # regions in tbl.regions
   region.1.count <- nrow(subset(tbl.regRegions, start >=  tss-5000 & stop <= tss))
   region.2.count <- nrow(subset(tbl.regRegions, start >=  tss & stop <= tss + 5000))
   checkTrue(region.1.count > 0)
   checkTrue(region.2.count > 0)
   checkEquals(region.1.count + region.2.count, nrow(tbl.regRegions))

} # test_build.trem2.twoRegions.motifDB.model
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
