# allKnowTFs
#------------------------------------------------------------------------------------------------------------------------
library(trenaSGM)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{

} # runTests
#------------------------------------------------------------------------------------------------------------------------
# see https://www.ebi.ac.uk/QuickGO/term/GO:0003700
test_DNA.binding.transcription.factor.activity <- function()
{
   printf("--- test_DNA.binding.transcription.factor.activity")

   tfs <- allKnownTFs(source="GO:DNAbindingTranscriptionFactorActivity", identifierType="geneSymbol")
   checkTrue(length(tfs) > 1500 & length(tfs) < 1600)  # 1551 in (oct 2018)

} # test_DNA.binding.transcription.factor.activity
#------------------------------------------------------------------------------------------------------------------------


