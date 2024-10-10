plotCorrDensity.SOSet <- function(object,
                                    filename = "pairwiseCorrDensityPlot.pdf",
                                    breaks = NULL,
                                    width = 15, height = 15) {
        mat <- object@corMat
        rank <- object@NMFrankL

        if (is.null(mat)) {
                stop("ERROR::Provide correlation matrix.")
        }

        if (!isSymmetric(mat)) {
                stop("ERROR::Check your correlation matrix - it is not symmetric.")
        } else {
                mat[upper.tri(mat, diag = FALSE)] <- NA
        }

        if (!is.null(rank)) {
                rank <- rank[which(names(rank) %in% unique(sapply(stringr::str_split(colnames(mat), "\\$"), "[[", 1)))]
                if (length(unique(sapply(stringr::str_split(colnames(mat), "\\$"), "[[", 1))) != length(rank)) {
                        stop("ERROR::The number of data and rank are not identical in the lists.")
                }
        } else {
                buff <- do.call(rbind, stringr::str_split(colnames(mat), "\\$"))
                colnames(buff) <- c("Cohort", "k", "rank")
                maxRanks <- aggregate(rank ~ Cohort, data = buff, max)
                minRanks <- aggregate(rank ~ Cohort, data = buff, min)

                rank <- list()
                for (idx in seq_along(maxRanks$Cohort)) {
                        rank[[idx]] <- c(as.numeric(minRanks$rank[idx]):as.numeric(maxRanks$rank[idx]))
                }
                names(rank) <- maxRanks$Cohort
                message("WARNING::No rank info provided - all ranks will be included.")
        }

        if (is.null(breaks)) {
                breaks = seq(-1, 1, 0.1)
        }

        df <- data.frame(
                dName = sapply(stringr::str_split(colnames(mat), "\\$"), "[[", 1),
                k = as.numeric(sapply(stringr::str_split(colnames(mat), "\\$"), "[[", 2)),
                rank = as.numeric(sapply(stringr::str_split(colnames(mat), "\\$"), "[[", 3))
        )

        node2include <- Reduce(union, lapply(seq_along(rank), function(idx) {
                return(intersect(which(df$dName == names(rank)[idx]), which(df$rank %in% rank[[idx]])))
        }))

        mat <- mat[node2include, node2include]
        colNames <- sapply(stringr::str_split(colnames(mat), "\\$"), "[[", 1)
        rowNames <- colNames
        dNames <- unique(colNames)

        if (length(dNames) > 0) {
                if (!is.null(filename)) pdf(filename, width = width, height = height)
                par(mfrow = c(length(dNames), length(dNames)))

                # distribution
                for (idx in seq_along(dNames)) {
                        for (jdx in seq_along(dNames)) {
                                if (jdx == length(dNames) && idx == 1) {
                                        correlations <- mat
                                        if (any(is.na(correlations))) correlations <- correlations[!is.na(correlations)]                                        

                                        hist(correlations,
                                                prob = TRUE, main = "",
                                                xlab = paste0("N = ", nrow(mat) * (nrow(mat) - 1) / 2, " (", nrow(mat), " * ", nrow(mat) - 1, " / 2)"),
                                                sub = "All pairs", xlim = c(-1, 1), breaks = breaks
                                        )
                                        lines(density(correlations), lwd = 2, col = "black")
                                } else {
                                        rowIdx <- which(rowNames == dNames[idx])
                                        colIdx <- which(colNames == dNames[jdx])
                                        if (jdx == 1) {
                                                ylab <- "Density"
                                        } else {
                                                ylab <- ""
                                        }
                                        correlations <- mat[rowIdx, colIdx]
                                        if (length(is.na(correlations)) > 0) correlations <- correlations[!is.na(correlations)]
                                        if (length(correlations) > 0) {
                                                hist(correlations, prob = TRUE, main = "", xlab = paste0(dNames[idx], "-", dNames[jdx]), ylab = ylab, xlim = c(-1, 1), breaks = breaks)
                                                lines(density(correlations), lwd = 2, col = "black")
                                        } else {
                                                plot.new()
                                        }
                                }
                        }
                }

                if (!is.null(filename)) dev.off()

                if (!is.null(filename)) message(paste0(filename, " was generated."))
        } else {
                stop("ERROR::Parser error, please check your dataset name(s).")
        }
}

#' Generate correlation coefficient density plot
#'
#' @param object A \code{SOSet} object
#' @param filename \code{character} Output file name
#' @param width \code{numeric} Width size for output PDF
#' @param height \code{numeric} Height size for output PDF
#'
#' @return Coefficient density plot
#'
#' @docType methods
#' @name plotCorrDensity
#' @rdname plotCorrDensity
#'
#' @examples
#' plotCorrDensity(SOSet)
#'
#' @export
#'
setMethod(f = "plotCorrDensity", signature = "SOSet", definition = plotCorrDensity.SOSet)

plotNetwork.SOTK <- function(object, weighted = TRUE, label = FALSE, 
                           annot = "cohort", 
                           vertexSize = 5, vertexLabelCex = 1, edgeAlpha = 0.2,
                           filename = "network.pdf",
                           width = 10, height = 10) {

        if (!is.numeric(vertexSize) || vertexSize == 0) {
                stop("ERROR::vertexSize should be greater than 0.")
        }

        if (!is.numeric(vertexLabelCex) || vertexLabelCex == 0) {
                stop("ERROR::vertexLabelCex should be greater than 0.")
        }

        if (!is.numeric(edgeAlpha) || edgeAlpha < 0 || edgeAlpha > 1) {
                stop("ERROR:edgeAlpha should be greater or equal to 0 and less or equal to 1.")
        }

        if (!annot %in% c("cohort", "community")) {
                stop("ERROR::annot should be either cohort or community.")
        }

        if (class(object) != "SOTK") {
                stop("ERROR::Provide a SOTK object which was created by the SOTK function.")
        } else {                
                graph <- object@corNetwork
                nodes <- igraph::V(graph)$name
                if (weighted) {
                        layout <- object@weightedLay
                } else {
                        layout <- object@unweightedLay
                }
                col <- object@SOSet@dataCol
                existingNames <- names(object@SOSet@NMFobjL)
                col <- col[which(existingNames %in% names(col))]
        }

        if (!is.null(filename)) pdf(filename, width = width, height = height)

        if (annot == "cohort") {
                igraph::V(graph)$color <- sapply(nodes, .assignColor, col = col)
                if (label) {
                        vertexLabel <- paste(
                                sapply(stringr::str_split(nodes, "\\$"), "[[", 2),
                                sapply(stringr::str_split(nodes, "\\$"), "[[", 3),
                                sep = "/"
                        )
                } else {
                        vertexLabel <- NA
                }

                plot(
                        graph,
                        layout = layout,
                        vertex.label = vertexLabel,
                        vertex.color = igraph::V(graph)$color,
                        vertex.label.color = "black",
                        vertex.size = vertexSize,
                        vertex.label.cex = vertexLabelCex,
                        edge.color = scales::alpha("grey80", edgeAlpha)
                )
                legend("topleft", legend = names(col), fill = col, bty = "n")
                legend("bottomleft", legend = paste0(igraph::vcount(graph), " nodes"), bty = "n")
        } else if (annot == "community") {
                communityInfo <- igraph::V(graph)$community
                commCols <- object@commCols
                if (label) {
                        vertexLabel <- communityInfo
                } else {
                        vertexLabel <- NA
                }

                plot(
                        graph,
                        layout = layout,
                        vertex.label = vertexLabel,
                        vertex.label.color = "black",
                        vertex.label.cex = vertexLabelCex,
                        vertex.color = commCols[communityInfo],
                        vertex.shape = "circle",
                        vertex.size = vertexSize,
                        edge.color = scales::alpha("grey80", edgeAlpha)
                )
                title(stringr::str_to_title(object@parameters$searchMet), cex.main = 2, col.main = "black")
                legend("bottomleft", legend = paste0(max(communityInfo), " communities"), bty = "n")
        } else {
                plot(
                        graph,
                        layout = layout,
                        vertex.label = NA,
                        vertex.color = "white",
                        vertex.size = vertexSize,
                        edge.color = scales::alpha("grey80", edgeAlpha)
                )
                legend("topleft", legend = names(col), fill = col, bty = "n")
                legend("bottomleft", legend = paste0(igraph::vcount(graph), " nodes"), bty = "n")
        }

        if (!is.null(filename)) dev.off()
}

.assignColor <- function(x, col) {
        dName <- unlist(strsplit(x, "\\$"))[1]
        return(col[dName])
}

#' Generate correlation network plot
#'
#' @param object A \code{SOTK} object
#' @param weighted \code{logical} Whether generate a correlation network with weighted or unweighted layout. Default is TRUE.
#' @param label \code{logical} Whether annotate metagene label or not. Default is FALSE.
#' @param annot \code{character} Whether data or community annotation on the network. Default is data.
#' @param vertexSize \code{numeric} Vertex/node size on the network. Default is 5.
#' @param vertexLabelCex \code{numeric} Label font size on the nodes. Default is 1.
#' @param edgeAlpha \code{numeric} Alpha value for the edge transparency. Default is 0.2.
#' @param filename \code{character} Output file name
#' @param width \code{numeric} Width size for output PDF
#' @param height \code{numeric} Height size for output PDF
#'
#' @return Correlation network plot
#'
#' @docType methods
#' @name plotNetwork
#' @rdname plotNetwork
#'
#' @examples
#' plotNetwork(object = SOTK)
#'
#' @export
#'
setMethod(f = "plotNetwork", signature = "SOTK", definition = plotNetwork.SOTK)

plotCommNetwork.SOTK <- function(object, vertexInfo,
                           filename = "community.pdf",
                           width = 10, height = 10) {

        graph <- object@commNetwork
        layout <- object@commLay
        community <- c(1:length(object@commCols))

        if (is.null(vertexInfo)) {
                if (!is.null(filename)) pdf(filename, width = width, height = height)
                plot(
                        graph, 
                        layout = layout, 
                        vertex.label = igraph::V(graph)$name, 
                        vertex.label.color = "black", 
                        vertex.label.cex = 1,
                        vertex.color = "white", 
                        vertex.size = igraph::V(graph)$size,
                        edge.width = igraph::E(graph)$weight
                )
                legend("topleft", legend=c(
                        "Node size = Number of metagenes",
                        "Edge tickness = Number of edges"
                ), bty="n")

                if (!is.null(filename)) dev.off()
        } else {
                vertexLabel <- vertexInfo[["vertexLabel"]]
                vertexSize <- vertexInfo[["vertexSize"]]
                vertexPie <- vertexInfo[["vertexPie"]]
                legend <- vertexInfo[["legend"]]
                legendCol <- vertexInfo[["legendCol"]]

                vertexPieCol <- list()
                observed <- do.call(rbind, vertexPie)
                for (i in community) {
                        vertexPieCol[[names(vertexPie)[i]]] <- c("white", legendCol)
                }

                chisq <- chisq.test(observed)
                residual <- chisq$residuals
                residual[is.na(residual)] <- 0
                residual[residual < 0] <- 0

                wColSum <- apply(residual, 1, sum)

                for(i in seq_along(vertexPie)) {
                        dummy <- 0
                        if (wColSum[i] == 0) dummy <- 1
                        vertexPie[[i]] <- c(dummy, residual[i,])
                }

                if (!is.null(filename)) pdf(filename, width = width, height = height)

                plot(
                        graph,
                        layout = layout,
                        vertex.label = vertexLabel,
                        vertex.label.color = "black",
                        vertex.label.cex = 1,
                        vertex.size = as.numeric(vertexSize),
                        vertex.shape = "pie",
                        vertex.pie = vertexPie,
                        vertex.pie.color = vertexPieCol,
                        edge.width = igraph::E(graph)$weight
                )
                legend("topleft", legend = legend, fill = legendCol, bty="n")

                if (!is.null(filename)) dev.off()
        }
}

#' Generate community network plot
#'
#' @param object A \code{SOTK} object
#' @param vertexInfo \code{list} Pie annotation that includes node label, size, pie (proportion), and legend.
#' @param filename \code{character} Output file name
#' @param width \code{numeric} Width size for output PDF
#' @param height \code{numeric} Height size for output PDF
#'
#' @return Correlation network plot
#'
#' @docType methods
#' @name plotCommNetwork
#' @rdname plotCommNetwork
#'
#' @examples
#' plotCommNetwork(object = SOTK, vertexInfo = list(vertexLabel = vertexLabel, vertexSize = vertexSize, vertexPie = vertexPie, legend = legend, legendCol = legendCol))
#'
#' @export
#'
setMethod(f = "plotCommNetwork", signature = "SOTK", definition = plotCommNetwork.SOTK)

statComm.SOTK <- function(object, filename = "stats.pdf", width = 10, height = 10) {
        datasetNames <- names(object@SOSet@NMFobjL)
        graph <- object@corNetwork
        community <- igraph::V(graph)$community
        names(community) <- igraph::V(graph)$name

        buff_cnt <- c()
        buff_metagene <- c()
        for (dName in datasetNames) {
                subComm <- community[which(sapply(stringr::str_split(names(community), "\\$"), "[[", 1) == dName)]
                
                cnt <- c()
                metagene <- c()
                for (whichCom in c(1:max(community))) {
                        subCommunities <- names(subComm)[which(subComm == whichCom)]
                        subCommunities <- paste(sapply(stringr::str_split(subCommunities, "\\$"), "[[", 2), sapply(stringr::str_split(subCommunities, "\\$"), "[[", 3), sep="/")
                        cnt <- c(cnt, length(subCommunities))
                        metagene <- c(metagene, paste(subCommunities, collapse=","))
                }
                buff_cnt <- rbind(buff_cnt, cnt)
                buff_metagene <- rbind(buff_metagene, metagene)
        }

        numberOfMetagenes <- t(buff_cnt)
        rownames(numberOfMetagenes) <- paste0("Community_", c(1:max(community)))
        colnames(numberOfMetagenes) <- datasetNames

        commMetagenes <- t(buff_metagene)
        rownames(commMetagenes) <- paste0("Community_", c(1:max(community)))
        colnames(commMetagenes) <- datasetNames

        pdf(filename, width = width, height = height)
        for (dName in datasetNames) {

                ranks <- object@SOSet@NMFrankL[[dName]]
                subComm <- community[which(sapply(stringr::str_split(names(community), "\\$"), "[[", 1) == dName)]                

                row <- list()
                for (whichCom in c(1:max(community))) {
                        metageneNames <- names(subComm)[which(subComm == whichCom)]
                        row[[whichCom]] <- as.numeric(sapply(stringr::str_split(metageneNames, "\\$"), "[[", 3))
                }

                prop <- lapply(row, function(x) {
                        rankCnt <- rep(0, length(ranks))
                        names(rankCnt) <- ranks

                        buffMat <- table(x)
                        for(idx in seq_along(buffMat)) {
                                rankName <- names(buffMat)[idx]
                                rankValue <- buffMat[idx]
                                rankCnt[rankName] <- rankValue
                        }

                        return(rankCnt)
                })

                input <- do.call(rbind, prop)
                rowSum <- apply(input, 1, sum)

                output <- c(); rnames <- c()
                for (idx in c(1:max(community))) {
                        prop <- input[idx,]/rowSum[idx]
                        if (!is.na(prop[1])) {
                                output <- rbind(output, prop)
                                rnames <- c(rnames, idx)
                        }
                }
                colnames(output) <- paste0("rank=", ranks)
                rownames(output) <- paste0("C_", rnames)

                df <- reshape2::melt(as.matrix(output))
                colnames(df) <- c("Community", "rank", "Proportion")

                stackedBar <- ggplot2::ggplot(df, ggplot2::aes(x = Community, y = Proportion, fill = rank)) +
                        ggplot2::geom_bar(position="stack", stat="identity") + 
                        ggplot2::scale_fill_manual(values = rainbow(ncol(output))) +
                        ggplot2::labs(
                                title = dName,
                                subtitle = "Rank distribution in each community",
                                x = "Community", y = "Proportion"
                        ) +
                        ggplot2::theme_bw() +
                        ggplot2::theme(
                                axis.line = ggplot2::element_line(colour = "black"),
                                panel.grid.major = ggplot2::element_blank(),
                                panel.grid.minor = ggplot2::element_blank(),
                                panel.border = ggplot2::element_blank(),
                                panel.background = ggplot2::element_blank(),
                                axis.text.x = ggplot2::element_text(angle = 0, vjust = 0, hjust = 0.5)
                        )
                print(stackedBar)
        }
        dev.off()

        return(list(numberOfMetagenes, commMetagenes))
}

#' Generate community network stats
#'
#' @param object A \code{SOTK} object
#' @param filename \code{character} Output file name
#' @param width \code{numeric} Width size for output PDF
#' @param height \code{numeric} Height size for output PDF
#'
#' @return Stacked bar plot and number of Metagenes per community in each dataset
#'
#' @docType methods
#' @name statComm
#' @rdname statComm
#'
#' @examples
#' statComm(object = SOTK)
#'
#' @export
#'
setMethod(f = "statComm", signature = "SOTK", definition = statComm.SOTK)

selectRank.SOTK <- function(object,
                           vertexSize = 8, vertexLabelCex = 1, edgeAlpha = 0.2,
                           filename = "network.pdf",
                           width = 10, height = 10) {

        if (!is.numeric(vertexSize) || vertexSize == 0) {
                stop("ERROR::vertexSize should be greater than 0.")
        }

        if (!is.numeric(vertexLabelCex) || vertexLabelCex == 0) {
                stop("ERROR::vertexLabelCex should be greater than 0.")
        }

        if (!is.numeric(edgeAlpha) || edgeAlpha < 0 || edgeAlpha > 1) {
                stop("ERROR:edgeAlpha should be greater or equal to 0 and less or equal to 1.")
        }
        
        if (class(object) != "SOTK") {
                stop("ERROR::Provide a SOTK object.")
        } else {                
                ranks <- object@SOSet@NMFrankL[[1]]
                graph <- object@corNetwork
                layout <- object@weightedLay
                nodes <- igraph::V(graph)$name
                communityInfo <- igraph::V(graph)$community
                commCols <- object@commCols
        }

        if (!is.null(filename)) pdf(filename, width = width, height = height)

        for (rank in ranks) {
                vertexLabel <- paste(
                        sapply(stringr::str_split(nodes, "\\$"), "[[", 2),
                        sapply(stringr::str_split(nodes, "\\$"), "[[", 3),
                        sep = "/"
                )
                vertexLabel[which(as.numeric(sapply(stringr::str_split(nodes, "\\$"), "[[", 3)) != rank)] <- ""

                vertexColor <- commCols[communityInfo]
                
                liveCols <- unique(vertexColor[which(as.numeric(sapply(stringr::str_split(nodes, "\\$"), "[[", 3)) == rank)])
                if (length(liveCols) > 0) {
                        vertexColor[-which(vertexColor %in% liveCols)] <- "white"
                }
                
                vertexColor[which(vertexLabel != "")] <- "black"
                liveColNode <- length(which(vertexColor == "black"))

                plot(
                        graph,
                        layout = layout,
                        vertex.label = vertexLabel,
                        vertex.color = vertexColor,
                        vertex.label.color = "white",
                        vertex.size = vertexSize,
                        vertex.label.cex = vertexLabelCex,
                        edge.color = scales::alpha("grey80", edgeAlpha)
                )
                title(paste0("Rank = ", rank), cex.main = 2, col.main = "black")
                legend("bottomleft",
                        legend = c(
                                paste0(liveColNode, " Metagenes"),
                                paste0(length(liveCols), " communities"),
                                paste0("Ratio = ", round(liveColNode / length(liveCols), 2))
                        ),
                        bty = "n"
                )
        }

        if (!is.null(filename)) dev.off()
}

#' Generate community network stats
#'
#' @param object A \code{SOTK} object
#' @param vertexSize \code{numeric} Vertex/node size on the network. Default is 8.
#' @param vertexLabelCex \code{numeric} Label font size on the nodes. Default is 0.75.
#' @param edgeAlpha \code{numeric} Alpha value for the edge transparency. Default is 0.2.
#' @param filename \code{character} Output file name
#' @param width \code{numeric} Width size for output PDF
#' @param height \code{numeric} Height size for output PDF
#'
#' @return Metagene marking at each rank on the community-weighted layout
#'
#' @docType methods
#' @name selectRank
#' @rdname selectRank
#'
#' @examples
#' selectRank(object = SOTK)
#'
#' @export
#'
setMethod(f = "selectRank", signature = "SOTK", definition = selectRank.SOTK)

# --- TO-DOs ---
# geoMean
# get signature genes
