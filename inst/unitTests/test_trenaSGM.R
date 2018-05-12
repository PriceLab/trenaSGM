library(RUnit)
library(trenaSGM)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("pfms"))
   pfms <- query2(MotifDb, c("sapiens", "jaspar2018"))

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

   sgm <- trenaSGM("hg38", "TREM2", mtx, pfms,
                   strategies=list(footprints="5000.5000.remapped"),
                   tfCorrelation=0.2)

   checkEquals(is(sgm), "trenaSGM")

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_trem2 <- function()
{
   printf("--- test_trem2")
   sgm <- trenaSGM("hg38", "TREM2", mtx, pfms,
                   strategies=list(footprints="5000.5000.remapped"),
                   tfCorrelation=0.3)
   x <- calculate(sgm)
   checkTrue(is.list(x))

} # test_trem2
#------------------------------------------------------------------------------------------------------------------------

