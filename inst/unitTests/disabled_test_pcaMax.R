library(RUnit)
library(trenaSGM)
library(MotifDb)
library(motifStack)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx")){
   filename <- system.file(package="trenaSGM", "extdata", "mayo.tcx.new.RData")
   if(!file.exists(filename))
      filename <- file.path("../extdata", "mayo.tcx.new.RData")
   stopifnot(file.exists(filename))
   load(filename)
   }


genome <- "hg38"
targetGene <- "TREM2"
tss <- 41163186
chromosome <- "chr6"

if(!exists("sgm"))
   sgm <- trenaSGM(genome, targetGene, quiet=FALSE)

if(!exists("tbl.enhancers")){
   source("~/github/trenaShinyApps/utils/getEnhancers.R")
   tbl.enhancers <- getEnhancerRegions(targetGene)
   }

#------------------------------------------------------------------------------------------------------------------------
bug_coryPcaMaxTREM2 <- function()
{
   printf("--- bug_coryPcaMaxTREM2")

   build.spec <- list(title="trem2.test",
                   type="footprint.database",
                   regions=tbl.enhancers,
                   tss=tss,
                   matrix=mtx,
                   db.host="khaleesi.systemsbiology.net",
                   databases=list("brain_hint_20", "brain_hint_16", "brain_wellington_20", "brain_wellington_16"),
                   motifDiscovery="builtinFimo",
                   tfMapping=c("MotifDB","TFClass"),
                   tfPrefilterCorrelation=0.4,
                   orderModelByColumn="pcaMax",
                   solverNames=c("lasso", "lassopv", "pearson", "sqrtlasso","randomForest", "ridge", "spearman"))

   strategies <- list(one=build.spec)
   test.model <- calculate(sgm, strategies)

   tbl.model <- test.model$one$model

    # problem demonstrated:
    # head(tbl.model, n=10)
    #     gene    betaLasso  lassoPValue pearsonCoeff betaSqrtLasso   rfScore   betaRidge spearmanCoeff concordance    pcaMax bindingSites
    # 19 NR6A1 -0.107859024 4.822947e-16   -0.4488069   -0.13632529  5.506683 -0.08333591    -0.4104492   0.4336171 2.2403429           98
    # 9  FOXP1 -0.057579174 1.383883e-11   -0.4549901   -0.07616693  5.352252 -0.05978206    -0.4441304   0.4258663 2.2321554          100
    # 8  FOXJ3  0.000000000 5.918305e-01   -0.4404496    0.00000000  3.691212 -0.02686905    -0.4200710   0.4507829 2.1419949          109
    # 18 MESP1 -0.039577314 1.533386e-07   -0.4128342   -0.06723244  5.283938 -0.04534640    -0.3832579   0.4222752 2.0942765           98
    # 13  IRF5  0.185352915 2.984581e-52    0.7662607    0.18354351 42.341548  0.10530790     0.7632320   0.4422842 1.2028771          119
    # 11 IKZF1  0.210014912 2.070879e-51    0.7638503    0.22211133 39.761465  0.10287181     0.7643571   0.4442307 1.1969145           62
    # 16  LYL1  0.138439201 4.683068e-36    0.7395869    0.13962218 28.544232  0.08795347     0.7307639   0.3878261 0.9088304           86
    # 22  SPI1  0.074030939 5.011984e-13    0.7075473    0.09973794 24.435249  0.09031790     0.7281323   0.3553523 0.7515626          202
    # 27  TFEC  0.006148505 1.448796e-18    0.7027833    0.00000000 18.992706  0.07090456     0.7045391   0.3216068 0.6368863          111
    # 2  CEBPA  0.103909212 1.524604e-17    0.6160476    0.12839944 15.379453  0.09230454     0.6283211   0.3278179 0.6306446           90

      # an informal test, showing the disagreement between rfScore and pcaMax
   genes.by.rf <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE), "gene"][1:5]
   genes.by.pcaMax <- tbl.model[order(tbl.model$pcaMax, decreasing=TRUE), "gene"][1:5]
   checkTrue(length(intersect(genes.by.rf, genes.by.pcaMax)) == 1)

     # now repeat with a smaller model, for more nimble exploration for the underlying problem

   tbl.regions <- data.frame(chrom=chromosome, start=tss-200, end=tss+2000, stringsAsFactors=FALSE)
   build.spec <- list(title="trem2.test",
                   type="footprint.database",
                   regions=tbl.regions,
                   tss=tss,
                   matrix=mtx,
                   db.host="khaleesi.systemsbiology.net",
                   databases=list("brain_hint_20"),
                   motifDiscovery="builtinFimo",
                   tfMapping=c("MotifDB"),
                   tfPrefilterCorrelation=0.4,
                   orderModelByColumn="pcaMax",
                   solverNames=c("lasso", "lassopv", "pearson", "sqrtlasso","randomForest", "ridge", "spearman"))

   build.spec.noBetaLasso <- build.spec
   build.spec.noBetaLasso$solverNames <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
   strategies <- list(bug=build.spec, fixed=build.spec.noBetaLasso)
   x <- calculate(sgm, strategies)

   tbl.model <- test.model$quick$model
   genes.by.rf <- tbl.model[order(tbl.model$rfScore, decreasing=TRUE), "gene"][1:6]
   genes.by.pcaMax <- tbl.model[order(tbl.model$pcaMax, decreasing=TRUE), "gene"][1:6]
   relative.ordering <- match(genes.by.rf, genes.by.pcaMax)  # 2 4 5 6 1 3


} # bug_coryPcaMaxTREM2
#------------------------------------------------------------------------------------------------------------------------
nobug <- function()
{
   chromosome <- "chr6"
   upstream <- 2000
   downstream <- 200
   tss <- 41163186

      # strand-aware start and end: trem2 is on the minus strand
   start <- tss - downstream
   end   <- tss + upstream
   tbl.regions <- data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

   build.spec <- list(title="fp.2000up.200down",
                      type="footprint.database",
                      regions=tbl.regions,
                      tss=tss,
                      matrix=mtx,
                      db.host="khaleesi.systemsbiology.net",
                      databases=list("brain_hint_20"),
                      motifDiscovery="builtinFimo",
                      tfMapping="MotifDB",
                      tfPrefilterCorrelation=0.4,
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   build.spec.2 <- build.spec
   build.spec.2$title <- "fp.2000up.200down.02"
   build.spec.2$tfPrefilterCorrelation=0.2
   build.spec.2$orderModelByColumn="pcaMax"
   build.spec.2$tfMapping <- c("MotifDb", "TFClass")

   strategies <- list(one=build.spec, two=build.spec.2)

   models <- calculate(sgm, strategies)



} # nobug
#------------------------------------------------------------------------------------------------------------------------
