library(AnnotationDbi)
library(BiocParallel)
library(org.Hs.eg.db)

lookup <- function(geneSymbols, dbfile) {
   db <- AnnotationDbi::loadDb(dbfile)
   on.exit(RSQLite::dbDisconnect(dbconn(db)))
   tbl.map <- AnnotationDbi::select(db, keys=geneSymbols, keytype="SYMBOL", columns="ENSEMBL")
   }

# tf.entrezIDs <- unique(unlist(get("GO:0003700", envir=org.Hs.egGO2ALLEGS)), use.names=FALSE)
# tf.geneSymbols <- unique(unlist(mget(tf.entrezIDs, envir=org.Hs.egSYMBOL), use.names=FALSE))
tf.geneSymbols <- unique( select(org.Hs.eg.db, "GO:0003700", "SYMBOL", "GOALL")$SYMBOL)
geneSymbols.list <- list(g1=tf.geneSymbols, g2=tf.geneSymbols, g3=tf.geneSymbols)
x <- bplapply(geneSymbols.list, lookup, dbfile(org.Hs.eg.db))


