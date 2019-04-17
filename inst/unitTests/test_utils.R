library(trenaSGM)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_roundNumericColumns()
   test_make.ensembl.geneSymbol.identifiers()

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
test_make.ensembl.geneSymbol.identifiers <- function()
{
   printf("--- test_make.ensembl.geneSymbol.identifiers")

   geneSymbols <-  toupper(c("CEBPA", "CEBPA", NA, "znf143", "MYC::MAX", "RARA(var.2)",  "Tcf12", NA))
   ensm.syms <- trenaSGM:::make.ensembl.geneSymbol.identifiers(geneSymbols)

   expected <- c("CEBPA|ENSG00000245848",
                 "CEBPA|ENSG00000245848",
                 "NA|NA",
                 "ZNF143|ENSG00000166478",
                 "MYC::MAX|NA",
                 "RARA(VAR.2)|NA",
                 "TCF12|ENSG00000140262",
                 "NA|NA")

   checkEquals(ensm.syms, expected)

} # test_make.ensembl.geneSymbol.identifiers
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
