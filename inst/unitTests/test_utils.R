library(trenaSGM)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_roundNumericColumns()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_roundNumericColumns <- function()
{

   load(system.file(package="trenaSGM", "extdata", "tbl.model.testFormatting.RData"))
   tbl.fixed <- roundNumericColumns(tbl.model.test, 3,
                                                   exponentialColumnNames = "lassoPValue")
   checkEquals(dim(tbl.fixed), dim(tbl.model.test))
   checkEquals(as.list(tbl.fixed[1,]),
               list(gene="MXI1",
                    betaLasso=0.395,
                    lassoPValue=1.64e-27,
                    pearsonCoeff=0.675,
                    rfScore=10.298,
                    betaRidge=0.106,
                    spearmanCoeff=0.593,
                    concordance=0.518,
                    pcaMax=2.766,
                    bindingSites=17))

} # test_roundNumericColumns
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
