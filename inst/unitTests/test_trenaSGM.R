library(RUnit)
library(trenaSGM)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx"))
   load(system.file(package="trenaSGM", "extdata", "mayo.tcx.RData"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_trem2()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")

   sgm <- trenaSGM("hg38", "TREM2")

   checkEquals(is(sgm), "trenaSGM")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_trem2 <- function()
{
   printf("--- test_trem2")

   sgm <- trenaSGM("hg38", "TREM2")

   fp.specs <- list(title="fp.5000up.5000down",
                    type="database.footprints",
                    matrix=mtx,
                    chrom="chr6",
                    tss=41163186,
                    upstream=5000,
                    downstream=5000,
                    db.host="khaleesi.systemsbiology.net",
                    databases=list("brain_hint_16", "brain_hint_20", "brain_wellington_16", "brain_hint_20"),
                    motifDiscovery="builtinFimo",
                    tfMapping=list("TFClass", "MotifDB"),
                    tfPrefilterCorrelation=0.2)

   strategies <- list(footprints=fp.specs)
   x <- calculate(sgm, strategies)
   checkTrue(is.list(x))

} # test_trem2
#------------------------------------------------------------------------------------------------------------------------

