require(data.table)
require(dplyr)
require(Rmutmod)
source("Rmutpheno.r")
gtfdir <- "~/projects/translateSelection/MC3/producedData/targetGtf/"
vardir <- "~/projects/translateSelection/MC3/producedData/variantAnnot/"
ws <- 1000L
modelpath <- "~/projects/translateSelection/MC3/src/test.rds"
.chr <- "chr1"
.vartype <- "syn"
phenodir <- "~/projects/translateSelection/MC3/producedData/summaryBPP/"
phenoCol <- "d"
minmut <- 3L
n <- 10000L
gspath <- "~/projects/translateSelection/MC3/producedData/geneSets.tsv"

rmutmod <- readRDS(modelpath)
cohort <- cohortGet(rmutmod)
mafdir <- mafdirGet(rmutmod)
cohort <- "HNSC"

flankdt <- maf2rnaFlanks(mafdir, .chr, cohort, .vartype, gtfdir, ws, minmut)
pmutdt <- flankdt2pmutdt(.chr, flankdt, rmutmod, .vartype, vardir, phenodir, phenoCol)

pmutdt <- split(pmutdt, by = "transcript_id")
pmutdt <- lapply(pmutdt, split, by = "mutid")
simdt <- lapply(pmutdt, redistMuts, n)
simdt <- data.table(
    pheno = unlist(simdt),
    transcript_id = rep(names(simdt), sapply(simdt, length)),
    sim = rep(1:n, length(simdt))
)

pvaldt <- phenoPvalue(simdt, mafdir, cohort, .vartype, phenodir, phenocol)
simpdt <- simPvalue(simdt)

setdt <- fread(gspath)
setdt <- setdt[
    ,
    .SD,
    .SDcols = c(
        "Transcript_ID",
        "cgc",
        "cgc-OG",
        "cgc-TSG",
        "wang",
        names(setdt)[grep(cohort, names(setdt))]
    )
]
# pretend the following is a batch of simulations for all genes
simpdt <- simpdt[sim %in% 1:1000]
simqdt <- simpdt[, list("score" = p.adjust(score, method = "fdr"), "transcript_id" = transcript_id), by = c("sim", "scoreType")]
simqdt[scoreType == "rpval", "scoreType" := "rqval"]
simqdt[scoreType == "lpval", "scoreType" := "lqval"]
setcolorder(simqdt, names(simpdt))
simsdt <- rbind(simpdt, simqdt)
simadt <- simAUC(simpdt, setdt)

# pretend simadt is already all aucs for all simulations and pvaldt is already the results for all genes
qvaldt <- pvaldt[
    ,
    list(
        "score" = p.adjust(score, method = "fdr"),
        "transcript_id" = transcript_id,
        "sumPheno" = sumPheno,
        "nmut" = nmut,
        "ntumor" = ntumor
    ),
    by = "scoreType"
]
qvaldt[scoreType == "rpval", "scoreType" := "rqval"]
qvaldt[scoreType == "lpval", "scoreType" := "lqval"]
setcolorder(qvaldt, names(pvaldt))
scoredt <- rbind(pvaldt, qvaldt)
aucdt <- aucPvalue(scoredt, simadt, setdt)

# after this in the pipeline there is only quality of life things left to do
# like put gene symbols and gene set label in a casted score table for easier read in excel or so
