#' @import data.table
#' @importFrom dplyr %>%

#' @export
maf2rnaFlanks <- function(mafdt, rnaGtf, ws) {

    # get ranges of windows to the flank of mutations
    mwRanges <- GenomicRanges::GRanges(
        mafdt$Transcript_ID,
        IRanges::IRanges(mafdt$Start_Position - ws, mafdt$Start_Position + ws),
        mutid = mafdt$id # identifier to track which resulting ranges correspond to each mutation
    )

    # get ranges of RNAs
    rnaRanges <- GenomicRanges::GRanges(
        rnaGtf$transcript_id,
        IRanges::IRanges(rnaGtf$Start_Position, rnaGtf$End_Position)
    )

    # find overlaps between each mutation's window and each exon of the same RNA
    mwrnaPov <- IRanges::findOverlapPairs(mwRanges, rnaRanges)

    # find the intersection of each mutation's window and each intersecting exon
    flankRanges <- GenomicRanges::pintersect(mwrnaPov)

    # organize results in data table
    flankdt <- data.table::data.table(
        mutid = S4Vectors::first(mwrnaPov)$mutid,
        transcript_id = as.character(GenomicRanges::seqnames(flankRanges)),
        start = GenomicRanges::start(flankRanges),
        end = GenomicRanges::end(flankRanges)
    )
    
    return(flankdt)

}


#' @export
flankdt2pmutdt <- function(rmutmod, flankdt, .chr) {

    pmutdt <- Rmutmod::mutdesign(rmutmod, flankdt, .chr)
    flankdt[, "rangeid" := 1:nrow(flankdt)]
    pmutdt[flankdt, ':=' ("mutid" = i.mutid, "transcript_id" = i.transcript_id), on = "rangeid"]

    return(pmutdt)

}

#' @export
annotateVartype <- function(.chr, pmutdt, vardb) {

    vardt <- vardb %>%
        dplyr::filter(seqnames == .chr) %>%
        dplyr::select(dplyr::all_of(c("position", "transcript_id", "pyrimidine", "pyrimidineMut", "type", "codon", "mutCodon"))) %>%
        dplyr::collect()
    
    pmutdt[
        vardt,
        ':=' ("vartype" = i.type, "codon" = i.codon, "mutCodon" = i.mutCodon),
        on = setNames(
            c("position", "transcript_id", "pyrimidine", "pyrimidineMut"),
            c("start", "transcript_id", "ref", "mut")
        )
    ]

    return()

}

#' @export
annotatePheno <- function(.chr, pmutdt, phenodb, phenocol, errorOnMiss = TRUE) {

    phenodt <- phenodb %>%
        dplyr::filter(seqnames == .chr) %>%
        dplyr::select(dplyr::all_of(c("position.abs", "transcript_id", "wt", "snp", phenocol))) %>%
        dplyr::collect()
    names(phenodt)[ncol(phenodt)] <- "pheno"

    # ensure phenotype mutations are pyrimidine oriented
    phenodt <- data.table::as.data.table(phenodt) # this is annoying but without this, strings are inmutable for some reason
    phenodt[
        Rmutmod::isPuri(wt, 0L),
        ':=' (
            "wt" = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(wt)), use.names = FALSE),
            "snp" = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(snp)), use.names = FALSE)
        )
    ]
    
    pmutdt[
        phenodt,
        "pheno" := i.pheno,
        on = setNames(
            c("position.abs", "transcript_id", "wt", "snp"),
            c("start", "transcript_id", "ref", "mut")
        )
    ]

    # handle missing data
    if (any(is.na(pmutdt$pheno))) {

        if (errorOnMiss) {

            stop("At least 1 missing phenotype for possible mutations in queried regions")

        } else {

            warning("At least 1 missing phenotype for possible mutations in queried regions")

        }

    }  

    return()

}

#' @export
sampleMut <- function(x, pmutdt, .cols) {

   UseMethod("sampleMut", x)

}

#' @export
sampleMut.numeric <- function(mu, pmutdt, .cols) {

    pmutdt[, "p" := mu]
    pmutdt[, "pnorm" := p / sum(p), by = "mutid"]

    n <- attr(mu, "n")
    simdt <- pmutdt[
        , .SD[.Internal(sample(.N, n, TRUE, pnorm))],
        by = "mutid",
        .SDcols = setdiff(.cols, "mutid")
    ]
    simdt[, "sim" := rep(1:n, length(unique(pmutdt$mutid)))]

    pmutdt[, ':=' ("p" = NULL, "pnorm" = NULL)]
    return(simdt)

}

#' @export
sampleMut.matrix <- function(mu, pmutdt, .cols) {

    LPS <- lapply(
        split(1:nrow(mu), pmutdt$mutid),
        function(idxs) mu[idxs, , drop = FALSE]
    )
    LPS <- lapply(LPS, function(M) t(M) / colSums(M))

    # redistribution simulation by window
    S <- lapply(
        LPS,
        function(M) apply(M, 1, function(p) .Internal(sample(length(p), 1, TRUE, p)))
    )

    # return simulated indices of pmutdt
    s <- mapply(function(M, idxs) colnames(M)[idxs], LPS, S, SIMPLIFY = FALSE)
    simdt <- pmutdt[as.integer(unlist(s)), .SD, .SDcols = .cols]
    simdt[, "sim" := rep(1:length(s[[1]]), length(s))]

    return(simdt)

}

#' @export
redistMut <- function(x, ...) {

   UseMethod("redistMut", x)

}

#' @export
redistMut.MutMatrixSim <- function(mutMatrixSim, pmutdt, .cols) {

    pmutdt[mutMatrixSim$fixef, "p" := i.density, on = c(mutMatrixSim$features, "mut")]
    mu <- pmutdt$p
    attr(mu, "n") <- mutMatrixSim$n
    simdt <- sampleMut(mu, pmutdt, .cols)
    
    return(simdt)

}

#' @export
redistMut.MonoMAFglmmTMBsim <- function(monoMAFglmmTMBsim, pmutdt, .cols) {

    # get linear predictor of mutation density for each window
    mu <- Rmutmod::linearPredictor(monoMAFglmmTMBsim, pmutdt)
    simdt <- sampleMut(mu, pmutdt, .cols)
    
    return(simdt)

}

#' @export
redistMut.MultiMAFglmmTMBsim <- function(multiMAFglmmTMBsim, pmutdt, .cols) {

    # separate mutations by mononucleotide substitution type
    pmutList <- split(pmutdt, by = c("ref", "mut"))
    names(pmutList) <- sub("[.]", "_", names(pmutList))

    # get linear predictor of mutation density per mutation type per simulation
    LPS <- mapply(
        Rmutmod::linearPredictor,
        multiMAFglmmTMBsim$sims[names(pmutList)],
        pmutList,
        SIMPLIFY = FALSE
    )

    # for now we will only support vectors
    mu <- setNames(
        unlist(LPS, recursive = FALSE, use.names = FALSE),
        unlist(lapply(pmutList, function(.dt) .dt$id), recursive = FALSE, use.names = FALSE)
    )
    n <- attr(LPS[[1]], "n")
    rm(pmutList, LPS)
    mu <- mu[as.character(pmutdt$id)]
    attr(mu, "n") <- n
    simdt <- sampleMut(mu, pmutdt, .cols)

    return(simdt)

}

phenoLoad <- function(phenodb, .cols, .chr = "all", utx = "all") {

    queryChr <- phenodb %>% dplyr::filter(seqnames %in% .chr)
    if (.chr[1] == "all") queryChr <- phenodb

    queryTx <- queryChr %>% dplyr::filter(transcript_id %in% utx)
    if (utx[1] == "all") queryTx <- queryChr

    phenodt <- queryTx %>%
        dplyr::select(dplyr::all_of(.cols)) %>%
        dplyr::collect()
    names(phenodt)[ncol(phenodt)] <- "pheno" # last column always assumed the phenotype column

    return(phenodt)

}

#' @export
obsPheno <- function(utx, mafdb, cohort, .vartype, phenodb, phenocol) {

    # load maf observed data of interest
    obsdt <- Rmutmod::mafLoad(
        mafdb,
        c("Start_Position", "Transcript_ID", "Tumor_Sample_Barcode", "codingRef", "codingMut"),
        cohort = cohort,
        .vartype = .vartype,
        utx = utx
    )
    
    # load phenotype data of interest
    phenodt <- phenoLoad(phenodb, c("position.abs", "transcript_id", "wt", "snp", phenocol), utx = utx)
    
    # get phenotype of observed mutations
    obsdt[
        phenodt,
        "pheno" := i.pheno,
        on = setNames(
            c("position.abs", "transcript_id", "wt", "snp"),
            c("Start_Position", "Transcript_ID", "codingRef", "codingMut")
        )
    ]

    return(obsdt)

}

#' @export
txObsPheno <- function(obsdt) {

    # sum observed phenotypes by tx
    txobsdt <- obsdt[
        ,
        list(
            "pheno" = sum(pheno),
            "nmut" = .N,
            "ntumor" = sum(!duplicated(Tumor_Sample_Barcode))
        ),
        by = "Transcript_ID"
    ]

    return(txobsdt)

}

#' @export
phenoPvalue <- function(simdt, mafdir, cohort, .vartype, phenodir, phenocol) {

    # get transcripts of interest (obviously, only the simulated ones) and number of simulations
    utx <- unique(simdt$transcript_id)
    nsim <- max(simdt$sim)

    # # load maf observed data of interest
    # mafdb <- arrow::open_dataset(mafdir)
    # mafdt <- mafdb %>% 
    #     dplyr::filter(
    #         Cohort == cohort,
    #         Chromosome == .chr,
    #         Transcript_ID %in% utx,
    #         Variant_Classification == .vartype
    #     ) %>%
    #     dplyr::select(
    #         dplyr::all_of(c("Start_Position", "Transcript_ID", "Tumor_Sample_Barcode", "codingRef", "codingMut"))
    #     ) %>%
    #     dplyr::collect()

    # # load phenotype data of interest
    # phenodb <- arrow::open_dataset(phenodir)
    # phenodt <- phenodb %>%
    #     dplyr::filter(seqnames == .chr, transcript_id %in% utx) %>%
    #     dplyr::select(dplyr::all_of(c("position.abs", "transcript_id", "wt", "snp", phenocol))) %>%
    #     dplyr::collect()
    # names(phenodt)[ncol(phenodt)] <- "pheno"

    # # get phenotype of observed mutations
    # mafdt[
    #     phenodt,
    #     "pheno" := i.pheno,
    #     on = setNames(
    #         c("position.abs", "transcript_id", "wt", "snp"),
    #         c("Start_Position", "Transcript_ID", "codingRef", "codingMut")
    #     )
    # ]

    # # sum observed phenotypes by tx
    # obsdt <- mafdt[
    #     ,
    #     list(
    #         "sumPheno" = sum(pheno),
    #         "nmut" = .N,
    #         "ntumor" = sum(!duplicated(Tumor_Sample_Barcode))
    #     ),
    #     by = "Transcript_ID"
    # ]

    # build pvalue table by comparing obesrved phenotypes per tx against their simulations
    simdt[obsdt, "obsPheno" := i.sumPheno, on = setNames("Transcript_ID", "transcript_id")]
    pvaldt <- simdt[
        ,
        list(
            "lpval" = (sum(pheno <= obsPheno) + 1L) / (nsim + 1L),
            "rpval" = (sum(pheno >= obsPheno) + 1L) / (nsim + 1L)
        ),
        by = "transcript_id"
    ]

    # treat each kind of pvalue as a different score scheme, this package expects this melted format in other functions
    pvaldt <- data.table::melt(
        pvaldt,
        id.vars = "transcript_id",
        measure.vars = c("lpval", "rpval"),
        variable.name = "scoreType",
        value.name = "score"
    )

    # add observed data info
    pvaldt[
        obsdt,
        ':=' (
            "sumPheno" = i.sumPheno,
            "nmut" = nmut,
            "ntumor" = ntumor
        ),
        on = setNames("Transcript_ID", "transcript_id")
    ]

    simdt[, "obsPheno" := NULL]
    return(pvaldt)

}

#' @export
simPvalue <- function(simdt) {

    nsim <- max(simdt$sim)
    utx <- unique(simdt$transcript_id)
    ntx <- length(utx)

    # each column in D is the distribution vector of simulated phenotypes for each tx
    D <- matrix(simdt$pheno, nrow = nsim, ncol = ntx)
    simLeft <- list()
    simRight <- list()
    for (ii in 1:nsim) {

        if (ii %% 1000 == 0) {

            gc()
            cat(ii, "/", nsim, "...\n", sep = "")

        }
        
        # treat each of the values of sim ii across tx as the "observed" to get a pvalue under the null per tx
        O <- matrix(
            rep(D[ii, ], each = nsim),
            nrow = nsim,
            ncol = ntx
        )
        # we don't sum +1 to numerator & denominator here because there is a -1 as we should remove from D the value we are testing in each iteration
        # denominator is just nsim but we can do that later for more efficiency
        simRight[[ii]] <- colSums(D <= O)
        simLeft[[ii]] <- colSums(D >= O)

        rm(O)

    }

    simpdt <- rbind(
        data.table::data.table(
            score = unlist(simRight, recursive = FALSE, use.names = FALSE),
            transcript_id = rep(utx, length(simRight)),
            sim = rep(1:nsim, each = ntx),
            scoreType = "rpval"
        ),
        data.table::data.table(
            score = unlist(simLeft, recursive = FALSE, use.names = FALSE),
            transcript_id = rep(utx, length(simLeft)),
            sim = rep(1:nsim, each = ntx),
            scoreType = "lpval"
        )
    )
    simpdt[, "score" := score / nsim]
   
   return(simpdt)

}

prrocWrapper <- function(score, setlabel, curveType) {

    auc <- NA_real_

    if (sum(setlabel) > 0) {

        auc <- switch(
            curveType,
            "roc" = PRROC::roc.curve(
                scores.class0 = score,
                weights.class0 = setlabel,
                sorted = TRUE
            )$auc,
            "pr" = PRROC::pr.curve(
                scores.class0 = score,
                weights.class0 = setlabel,
                sorted = TRUE,
                dg.compute = FALSE
            )$auc.integral
        )

    }

    return(auc)

}

score2auc <- function(sdt, setdt, aggcols) {

    # calculate the observed aucs per gene set
    sdt[, "setlabel" := NA_integer_]
    adt <- list()
    for (ii in 2:ncol(setdt)) { # first col is tx id

        sdt[, "setlabel" := NULL]

        sdt[
            setNames(
                setdt[, .SD, .SDcols = c(1, ii)],
                c("transcript_id", "setlabel")
            ),
            "setlabel" := i.setlabel,
            on = "transcript_id"
        ]

        # calculate AUCs
        adt[[names(setdt)[ii]]] <- sdt[
            ,
            list(
                "roc" = prrocWrapper(score, setlabel, "roc"),
                "prc" = prrocWrapper(score, setlabel, "pr"),
                "TP" = sum(setlabel),
                "TN" = sum(!setlabel)
            ),
            by = aggcols
        ]

    }

    adt <- cbind(
        rbindlist(adt),
        geneSet = rep(names(adt), sapply(adt, nrow))
    )
    adt <- data.table::melt(
        adt,
        measure.vars = c("roc", "prc"),
        variable.name = "curveType",
        value.name = "auc"
    )
    
    sdt[, "setlabel" := NULL]
    return(adt)

}

#' @export
simAUC <- function(simsdt, setdt) {

    simadt <- score2auc(simsdt, setdt, c("sim", "scoreType"))
    simadt[, ':=' ("TP" = NULL, "TN" = NULL)]

    return(simadt)
    
}

#' @export
aucPvalue <- function(scoredt, simadt, setdt) {

    # get number of simulations
    nsim <- max(simadt$sim)

    # get AUCs of observed scores
    obsdt <- score2auc(scoredt, setdt, "scoreType")

    # build pvalue table by comparing obesrved AUCs against simulations
    simadt[
        obsdt,
        "obsAUC" := i.auc,
        on = c("scoreType", "curveType", "geneSet")
    ]
    aucdt <- simadt[
        ,
        list("rpval" = (sum(auc >= obsAUC) + 1L) / (nsim + 1L)),
        by = c("scoreType", "curveType", "geneSet")
    ]

    # add auc effect size and geneSet counts
    aucdt[
        obsdt,
        ':=' (
            "auc" = i.auc,
            "TP" = i.TP,
            "TN" = i.TN
        ),
        on = c("scoreType", "curveType", "geneSet")
    ]

    simadt[, "obsAUC" := NULL]
    return(aucdt)

}
