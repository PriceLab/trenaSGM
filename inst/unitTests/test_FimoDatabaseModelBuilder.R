library(RUnit)
library(trenaSGM)
library(org.Hs.eg.db)
library(TrenaProjectErythropoiesis)
#------------------------------------------------------------------------------------------------------------------------
Sys.setlocale("LC_ALL", "C")

if(!exists("tp"))
   tp <- TrenaProjectErythropoiesis()

if(!exists("mtx")){
   mtx.name <- "brandLabDifferentiationTimeCourse-16173x28"
   mtx <- getExpressionMatrix(tp, mtx.name)
   }

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_simpleBuild()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   genome <- "hg38"
   targetGene <- "GATA2"
   chromosome <- "chr3"
   tss <- 128493185
   upstream <- 5000
   downstream <- 5000
      # strand-aware start and end: GATA2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

      # a very temporary fimo database, using sqlite
   db.name <- system.file(package="TrenaProjectErythropoiesis", "extdata", "fimoDBs",
                          "gata2.gh.fimoBindingSites.sqlite")

   checkTrue(file.exists(db.name))

   build.spec <- list(title="fimo.5000up.5000down",
                      type="fimo.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene,
                      tss=tss,
                      matrix=mtx,
                      db.host="localhost",
                      db.port=NA_integer_,
                      databases=list(db.name),
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="fimoDatabase",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      maxModelSize=10,
                      fimoPValueThreshold=0.0002,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- FimoDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)

   checkTrue("FimoDatabaseModelBuilder" %in% is(builder))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_simpleBuild <- function()
{
   printf("--- test_simpleBuild")

   genome <- "hg38"
   targetGene <- "GATA2"
   chromosome <- "chr3"
   tss <- 128493185
   upstream <- 5000
   downstream <- 5000
      # strand-aware start and end: GATA2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   db.name <- system.file(package="TrenaProjectErythropoiesis", "extdata", "fimoDBs", "gata2.gh.fimoBindingSites.sqlite")
   checkTrue(file.exists(db.name))

   build.spec <- list(title="fimo.5000up.5000down",
                      type="fimo.database",
                      regions=tbl.regions,
                      geneSymbol=targetGene,
                      tss=tss,
                      matrix=mtx,
                      db.host="localhost",
                      db.port=NA_integer_,
                      databases=list(db.name),
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      motifDiscovery="fimoDatabase",
                      tfPool=allKnownTFs(),
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      maxModelSize=10,
                      fimoPValueThreshold=0.0002,
                      orderModelByColumn="spearmanCoeff",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   builder <- FimoDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   checkTrue("FimoDatabaseModelBuilder" %in% is(builder))
   x <- build(builder)
   checkTrue(all(c("model", "regulatoryRegions") %in% names(x)))
   tbl.model <- x$model
   checkEquals(nrow(tbl.model), 10)

       # some stochasticity with rfScore.  but with spearmanCoeff ranking, two
       # top contenders seem always to be present

   checkTrue(all(c("CBFB", "ETV4") %in% tbl.model$gene))
   tbl.reg <- x$regulatoryRegions
   checkEquals(sum(tbl.model$bindingSites), nrow(tbl.reg))

      # ensure no duplicated binding sites, which can arise from similar motifs
      # from different motif sources
   bs.sig <- with(tbl.reg, sprintf("%s:%d-%d:%s:%s", chrom, start, end, strand, tf))
   checkEquals(length(which(duplicated(bs.sig))), 0)


      # now a far more stringent fimo pvalue
   build.spec$fimoPValueThreshold <- 1e-5
   build.spec$maxModelSize <- 5
   builder <- FimoDatabaseModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x2 <- build(builder)
   tbl.model <- x2$model
   tbl.reg   <- x2$regulatoryRegions
   checkEquals(nrow(tbl.model), 5)
   checkTrue(all(tbl.reg$pValue <= build.spec$fimoPValueThreshold))

} # test_simpleBuild
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
