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
	for x in inst/unitTests/test_*.R; do echo $x; R -f $$x; done

