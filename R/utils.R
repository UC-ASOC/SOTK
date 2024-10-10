#' Import cNMF results then export to an NMF object
#'
#' @param prefix \code{character} The run_name parameter from cNMF run
#' @param metrics \code{character} Either score or TPM, default is score
#' @param minRank \code{numeric} The minimum rank from cNMF run
#' @param maxRank \code{numeric} the maximum rank from cNMF run
#' @param by \code{numeric} increment of the sequence
#' @param denThre \code{character} The density threshold in string with a underscore
#' @param cnmfDir \code{character} Directory to the cNMF results
#'
#' @return An NMF object
#'
#' @examples
#' nmfObj <- importCNMF(prefix = "PREFIX", metrics = "score", minRank = 2, maxRank = 10, by = 1, denThre = "2_0", cnmfDir = "./cNMF")
#'
#' @export
#' @importFrom NMF nmfModel
#' @importFrom methods new
importCNMF <- function(prefix = NULL, metrics = "score",
                       minRank = NULL, maxRank = NULL, by = 1, denThre = NULL,
                       cnmfDir = "NULL", verbose = TRUE) {
        if (is.null(prefix) || is.null(minRank) || is.null(maxRank) || is.null(denThre) || is.null(cnmfDir)) {
                stop("ERROR::Please provide all required arguments: prefix, minRank, maxRank, denThre, cnmfDir")
        }

        if (!any(metrics %in% c("score", "TPM"))) {
                stop("ERROR::Please specify either score or TPM")
        }

        if (!is.numeric(minRank) || !is.numeric(maxRank) || !is.numeric(by) || by < 1) {
                stop("ERROR::Please provide integer for both minRank/maxRank/by and/or by should be greater than 1")
        }

        resObj <- list()
        hitK <- c()
        for (k in seq(minRank, maxRank, by)) {
                if (file.exists(file.path(cnmfDir, paste0(prefix, ".gene_spectra_tpm.k_", k, ".dt_", denThre, ".txt"))) &&
                        file.exists(file.path(cnmfDir, paste0(prefix, ".gene_spectra_score.k_", k, ".dt_", denThre, ".txt"))) &&
                        file.exists(file.path(cnmfDir, paste0(prefix, ".usages.k_", k, ".dt_", denThre, ".consensus.txt")))) {
                        if (verbose) message(paste0("\nLoading cNMF result files at k = ", k))

                        if (metrics == "TPM") {
                                wUsage <- t(read.table(file.path(cnmfDir, paste0(prefix, ".gene_spectra_tpm.k_", k, ".dt_", denThre, ".txt")),
                                        header = T, row.names = 1, sep = "\t", stringsAsFactor = F
                                ))
                        } else if (metrics == "score") {
                                wUsage <- t(read.table(file.path(cnmfDir, paste0(prefix, ".gene_spectra_score.k_", k, ".dt_", denThre, ".txt")),
                                        header = T, row.names = 1, sep = "\t", stringsAsFactor = F
                                ))
                        }

                        colnames(wUsage) <- NULL
                        if (verbose) message(paste("\tUsage (W): ", nrow(wUsage), " genes X ", ncol(wUsage), sep = ""))

                        hScore <- t(read.table(file.path(cnmfDir, paste0(prefix, ".usages.k_", k, ".dt_", denThre, ".consensus.txt")),
                                header = T, row.names = 1, sep = "\t", stringsAsFactor = F
                        ))
                        rownames(hScore) <- NULL
                        if (verbose) message(paste("\tScore (H): ", nrow(hScore), " X ", ncol(hScore), " samples", sep = ""))

                        fitRes <- methods::new("NMFfitX1")
                        attr(fitRes, "fit") <- NMF::nmfModel(model = "NMFstd", W = wUsage, H = hScore)

                        if (length(resObj) < 1) {
                                resObj <- list(fitRes)
                        } else {
                                resObj <- append(resObj, fitRes)
                        }
                        hitK <- c(hitK, k)
                } else {
                        warning(paste0("WARNING::cNMF result files (i.e., *gene_spectra_*, *usage*) do not exist at k = ", k))
                }
        }

        if (length(hitK) > 0) {
                names(resObj) <- as.character(hitK)

                cnmf <- list(fit = resObj)
                class(cnmf) <- "NMF.rank"
        } else {
                stop("ERROR::No result found.")
        }

        return(cnmf)
}

#' Merge multiple NMF objects into one
#'
#' @param nmfObjL \code{list} A list of NMF objects
#'
#' @return An NMF object
#'
#' @examples
#' nmfObj <- mergeNMFObjs(nmfObjL = objL)
#'
#' @export
#' @importFrom NMF compare connectivity
mergeNMFObjs <- function(nmfObjL) {
        measures <- compare(nmfObjL)
        measures <- measures[, c(5:ncol(measures))]
        rownames(measures) <- c(1:nrow(measures))

        consensus <- lapply(nmfObjL, function(res) {
                cons <- connectivity(res)
                attr(cons, "model") <- NULL
                return(cons)
        })

        fit <- nmfObjL

        obj <- list(
                measures = measures,
                consensus = consensus,
                fit = fit
        )
        class(obj) <- "NMF.rank"

        return(obj)
}

#' Get the most variable genes (MVGs)
#'
#' @param profile code{matrix} The expression profile
#' @param coefVar \code{numeric} A threshold value for the coefficient variation, default is 0.1
#' @param no \code{numeric} Number of genes to return, default is 2000
#'
#' @return A vector of genes
#'
#' @examples
#' nmfObj <- getMVGs(profile = expr, coefVar = 0.1, no = 2000)
#'
#' @export
#' @importFrom stringr str_split
#' @importFrom EnvStats cv
getMVGs <- function(profile = NULL, coefVar = 0.1, no = 2000) {
        if (is.null(profile)) {
                stop("ERROR::Please provide a expression profile where rows are genes and columns are samples")
        }

        if (coefVar < 0 || no < 0) {
                stop("ERROR::Please provide a positive value for coefVar")
        }

        coefficientVariation <- apply(profile, 1, cv)
        if (length(which(coefficientVariation < coefVar)) > 0) {        
                profile <- profile[-which(coefficientVariation < coefVar), ]
        }
        if (nrow(profile) < 1) {
                steop("ERROR::Please adjust coefVar")
        }
        
        standardDeviation <- apply(profile, 1, sd)
        profile <- profile[order(standardDeviation, decreasing=T),]
        
        return(rownames(profile)[c(1:no)])
}
