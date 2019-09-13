quick:  docs install

all:  docs vig build install check biocCheck

docs:
	R -e "devtools::document()"
vig:
	R -e "devtools::build_vignettes()"

build:
	(cd ..; R CMD build trenaSGM)

install:
	(cd ..; R CMD INSTALL --no-test-load trenaSGM)

check:
	(cd ..; R CMD check `ls -t trenaSGM_* | head -1`)

biocCheck:
	(cd ..; R CMD BiocCheck `ls -t trenaSGM_* | head -1`)

test:
	R -f inst/unitTests/test_allKnownTFs.R
	R -f inst/unitTests/test_utils.R
	R -f inst/unitTests/test_trenaSGM.R
	R -f inst/unitTests/test_ModelBuilder.R
	R -f inst/unitTests/test_NoDnaModelBuilder.R
	R -f inst/unitTests/test_FootprintDatabaseModelBuilder.R
	R -f inst/unitTests/test_RegionsFimoMatchingModelBuilder.R
	# R -f inst/unitTests/test_FimoDatabaseModelBuilder.R   # deferred: depends on FimoService running
