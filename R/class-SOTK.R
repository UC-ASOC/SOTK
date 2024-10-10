setOldClass("igraph")

#' @include class-SOSet.R
#' @import igraph
#' @rdname SOTK
#' @export
#'
setClass(
        "SOTK", 
        slots = c(
                SOSet = "SOSet", # NMF objects/results and correlation matrix in SOSet class
                sample2metagene = "list", # assign samples to metagene based on their max usage values
                corNetwork = "igraph", # Correlation network in igraph class
                unweightedLay = "matrix", # Correlation network layout (unweighted), i.e., X and Y coordinates
                weightedLay = "matrix", # Weighted layout based on community search algorithm, i.e., X and Y coordinates
                commNetwork = "igraph", # Correlation network in igraph class
                commLay = "matrix", # Correlation network layout (unweighted), i.e., X and Y coordinates
                commCols = "vector", # Colors that represent each community
                parameters = "list" # Stash parameters
        )
)
# TO-DOs: setValidity("SOTK", function(object) {}

#' @title SOTK object and constructor
#' 
#' @description An S4 class to contain Mx input and output
#' 
#' @details A
#'
#' @param SOSet A \code{SOSet} of the NMF results in SOSet class.
#' @param coefThre A \code{numeric} of the coefficient threshold value to include node in the correlation network. Default is 0.5.
#' @param seed A \code{numeric} of the seed to generate a correlation network. Default is 1118.
#' @param niter A \code{numeric} of the number of iteration when the igraph package generates a correlation network. Default is 1000.
#' @param drop \code{logical} Whether drop small/isolated island(s) from the correlation network. Default is FALSE.
#' @param searchMet \code{character} of the community search method, either greedy (fast greedy), leiden, betweeness, randomwalk (random walk), eigen, or louvain. Default is greedy.
#' @param commWeight \code{numeric} of the weight value on community info. If users want to rearrange nodes that identified as the same community closer, set a higher value. Default is 100.
#' @param cohortWeight \code{numeric} of the weight value on data/cohort (input). If users want to rearrange nodes closer for each data/cohort, set a higher value. Default is 5.
#' 
#' @return SOTK object
#'
#' @docType class
#' @rdname SOTK
#' @import igraph stringr RColorBrewer NMF
#' @export 
#'
SOTK <- function(SOSet, coefThre = 0.5, seed = 123456, niter = 1000, drop = FALSE, searchMet = "greedy", commWeight = 100, cohortWeight = 10) {
        if (is.null(SOSet) || class(SOSet) != "SOSet") {
                stop("ERROR::Provide an SOSet object.")
        }

        if (coefThre < -1 || coefThre > 1) {
                stop("ERROR::Provide coefThre between -1 and 1.")
        }        

        if (!is.numeric(seed)) {
                stop("ERROR::Provide an interger for seed.")
        }

        if (!is.numeric(niter)) {
                stop("ERROR::Provide an interger for niter (number of iteration) to generate a network.")
        }

        if (is.null(searchMet)) {
                stop("ERROR::Please provide either greedy, leiden, betweeness, randomwalk, eigen, or louvain.")
        } else if (!searchMet %in% c("greedy", "leiden", "betweeness", "randomwalk", "eigen", "louvain")) {
                stop("ERROR::Please provide either greedy, leiden, betweeness, randomwalk, eigen, or louvain.")
        }

        if (!is.numeric(commWeight)) {
                stop("ERROR::Provide an interger for commWeight (Weight on community info) to update the layout.")
        }

        if (!is.numeric(cohortWeight)) {
                stop("ERROR::Provide an interger for cohortWeight (Weight on cohort info) to update the layout.")
        }

        if (isSymmetric(SOSet@corMat)) {
                mat <- SOSet@corMat
                mat[mat < coefThre] <- 0

                rSum <- apply(mat, 1, sum)
                cSum <- apply(mat, 2, sum)
                exclude <- rSum + cSum

                if (length(which(exclude == 0)) > 0) {
                        mat <- mat[-which(exclude == 0), -which(exclude == 0)]
                }

                if (nrow(mat) == 0) {
                        stop("ERROR::No coefficient values above the coefThre. Lower your coef treshold value.")
                } else {
                        # ---- Create a correlation network
                        message(paste0("Seed: ", seed))
                        set.seed(seed)
                        graph <- igraph::graph.adjacency(mat, weighted = TRUE, mode = "lower")
                        graph <- igraph::as.undirected(graph)
                        if (drop) {
                                cl <- igraph::components(graph)
                                graph <- igraph::delete_vertices(graph, igraph::V(graph)$name[cl$membership %in% which(cl$csize < mean(cl$csize))])
                        }
                        unweightedLayout <- igraph::layout.fruchterman.reingold(graph, niter = niter)
                        
                        # ---- Community search                        
                        nodes <- igraph::V(graph)$name
                        
                        message("\nCommunity search algorithm:")                        
                        set.seed(seed)
                        if (searchMet == "greedy") {
                                fc <- igraph::cluster_fast_greedy(graph)
                                message("\tFast Greedy")
                        } else if (searchMet == "leiden") {
                                fc <- igraph::cluster_leiden(graph)
                                message("\tLeiden")
                        } else if (searchMet == "betweeness") {
                                fc <- igraph::cluster_edge_betweenness(graph)
                                message("\tBetweeness")
                        } else if (searchMet == "randomwalk") {
                                fc <- igraph::cluster_walktrap(graph)
                                message("\tRandom Walk")
                        } else if (searchMet == "eigen") {
                                fc <- igraph::cluster_leading_eigen(graph)
                                message("\tEigen")
                        } else if (searchMet == "louvain") {
                                fc <- igraph::cluster_louvain(graph)
                                message("\tLouvain")
                        }

                        igraph::V(graph)$community <- fc$membership

                        # ---- New/Update layout
                        newGraph <- graph
                        igraph::E(newGraph)$weight <- 1

                        message("\nUpdating weights for")
                        for (comm in sort(unique(igraph::V(graph)$community))) {
                                message(paste0("\tCommunity #", comm))
                                eachComm <- igraph::V(graph)$name[which(igraph::V(graph)$community == comm)]
                                if (length(eachComm) > 1) {
                                        newGraph <- igraph::add_edges(newGraph, combn(eachComm, 2), attr = list(weight = commWeight))
                                }
                        }

                        message("\nUpdating weights for")
                        for (dName in unique(sapply(stringr::str_split(igraph::V(graph)$name, "\\$"), "[[", 1))) {
                                message(paste0("\tData: ", dName))
                                eachDset <- igraph::V(graph)$name[which(sapply(stringr::str_split(nodes, "\\$"), "[[", 1) == dName)]
                                if (length(eachDset) > 1) {
                                        newGraph <- igraph::add_edges(newGraph, combn(eachDset, 2), attr = list(weight = cohortWeight))
                                }
                        }

                        message("\nCalculating new layout based on new weights.")
                        set.seed(seed)
                        weightedLayout <- igraph::layout.fruchterman.reingold(newGraph, niter = niter)

                        # ---- Community colors
                        qualColPals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
                        colVector <- unlist(mapply(RColorBrewer::brewer.pal, qualColPals$maxcolors, rownames(qualColPals)))

                        if (max(fc$membership) > length(colVector)) {
                                commCols <- sample(colVector, max(fc$membership), replace = TRUE)
                        } else {
                                commCols <- sample(colVector, max(fc$membership), replace = FALSE)
                        }
                        names(commCols) <- c(1:max(fc$membership))

                        # ---- Sample2metagene or metagene2Sample
                        sample2metagene <- .assignSample(NMFobjL = SOSet@NMFobjL)

                        # ---- Community-level network
                        commInfo <- .generateCommunityNetwork(graph = graph, seed = seed, niter = niter)
                        commNetwork <- commInfo[["commNetwork"]]
                        commLay <- commInfo[["commLay"]]                        
                }

                SOTK <- new(
                        "SOTK",
                        SOSet = SOSet,
                        sample2metagene = sample2metagene,
                        corNetwork = graph,
                        unweightedLay = unweightedLayout,
                        weightedLay = weightedLayout,
                        commNetwork = commNetwork,
                        commLay = commLay, 
                        commCols = commCols,
                        parameters = list(coefThre = coefThre, seed = seed, niter = niter, drop = drop, searchMet = searchMet, commWeight = commWeight, cohortWeight = cohortWeight)
                )
        } else {
                stop("ERROR::Check correlation matrix (corMat) in SOSet.")
        }

        return(SOTK)
}
# TO-DOs: generic method for community search algorithm to include some parameters, e.g., resolution_parameter

.assignSample <- function(NMFobjL, perSample = FALSE) {
        if (is.null(NMFobjL)) {
                stop("ERROR::Check NMFobjL parameter.")
        }

        clusteringInfo <- lapply(seq_along(NMFobjL), function(x) {
                dName <- names(NMFobjL)[x]
                obj <- NMFobjL[[x]]

                clInfo <- c()
                for (k in as.numeric(names(obj$fit))) {
                        usage <- NMF::coef(get(as.character(k), obj$fit))
                        cl <- apply(usage, 2, function(x) {
                                return(which(x == max(x)))
                        })
                        rel <- data.frame(
                                sample = names(cl),
                                metagenes = paste0(dName, "$", sprintf("%02d", cl), "$", sprintf("%02d", k))
                        )
                        clInfo <- rbind(clInfo, rel)
                }

                if (!perSample) {
                        reForm <- list()
                        for (metagene in sort(unique(clInfo$metagenes))) {
                                samples <- clInfo$sample[which(clInfo$metagenes == metagene)]
                                if (is.null(reForm[[metagene]])) {
                                        reForm[[metagene]] <- samples
                                } else {
                                        reForm[[metagene]] <- append(reForm[[metagene]], samples)
                                }
                        }
                        clInfo <- reForm
                }

                return(clInfo)
        })
        names(clusteringInfo) <- names(NMFobjL)

        return(clusteringInfo)
}

.generateCommunityNetwork <- function(graph, seed, niter) {
        community <- igraph::V(graph)$community
        v <- aggregate(community, by=list(community), FUN=length)
        colnames(v) <- c("name", "size")
        names(community) <- igraph::V(graph)$name

        edgeList <- igraph::as_edgelist(graph, names = TRUE)
        edgeComm <- c()
        for(idx in c(1:nrow(edgeList))) {
                oneEdge <- edgeList[idx,]
                fromTo <- sort(community[oneEdge])
                if (fromTo[1] != fromTo[2]) {
                        edgeComm <- rbind(edgeComm, paste0(fromTo, collapse="--"))
                }
        }
        if (is.null(edgeComm)) {
                warning("WARNING::No edges among communities - lower your coefficient threshold value.")
        } else {
                edgeComm <- as.data.frame(edgeComm)
                colnames(edgeComm) <- "edge"
                e <- aggregate(edgeComm$edge, by=list(edgeComm$edge), FUN=length)
                colnames(e) <- c("edge", "weight")
        }

        nodeSize <- as.numeric(log2((v[,2]+1)^4)) # scaling... for better visualization
        nodeSize[nodeSize < 5] <- 5

        commNetwork <- igraph::graph.empty()
        commNetwork <- igraph::add.vertices(commNetwork, nrow(v), 
                name=as.character(v[,1]), 
                size=nodeSize
        )
        if (!is.null(edgeComm)) {
                commNetwork <- igraph::add.edges(commNetwork, t(matrix(unlist(stringr::str_split(e$edge, "--")), ncol = 2, byrow = TRUE)))
                commNetwork <- igraph::as.undirected(commNetwork)
                igraph::E(commNetwork)$weight <- log10(e[,2] + 1) * 3 # scaling... for better visualization
        }
        
        set.seed(seed)
        commLay <- igraph::layout.fruchterman.reingold(commNetwork, niter = niter)

        message("\nCommunity-level network generated.\n")

        return(list(commNetwork = commNetwork, commLay = commLay))
}

#' Show a SOTK
#'
#' @param object \code{SOTK}
#'
#' @return Prints the SOTK object to the output stream, and returns invisible NULL.
#' @export
#'
setMethod(
        "show",
        signature = signature(object = "SOTK"),
        function(object) {
                graph <- object@corNetwork
                nodes <- igraph::V(graph)$name
                cat("Correlation network:\n")
                cat("\tNodes       : ", length(nodes), "\n", sep="")
                cat("\tCommunities : ", length(object@commCols), " identified\n", sep="")
                cat("Parameters:\n")
                cat("\tcoefThre    : ", object@parameters[["coefThre"]], "\n", sep="")
                cat("\tseed        : ", object@parameters[["seed"]], "\n", sep="")
                cat("\tniter       : ", object@parameters[["niter"]], "\n", sep="")
                cat("\tdrop        : ", object@parameters[["drop"]], "\n", sep="")
                cat("\tsearchMet   : ", object@parameters[["searchMet"]], "\n", sep="")
                cat("\tcommWeight  : ", object@parameters[["commWeight"]], "\n", sep="")
                cat("\tcohortWeight: ", object@parameters[["cohortWeight"]], "\n", sep="")
        }
)
