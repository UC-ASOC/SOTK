setGeneric(
        "plotCorrDensity", 
        function(
                object, 
                filename = "pairwiseCorrDensityPlot.pdf", 
                breaks = NULL,
                width = 15, 
                height = 15
        )
        standardGeneric("plotCorrDensity")
)

setGeneric(
        "plotNetwork", 
        function(
                object, 
                weighted = TRUE, 
                label = FALSE,
                annot = c("data", "community", "layout"),
                vertexSize = 5, 
                vertexLabelCex = 1, 
                edgeAlpha = 0.2,
                filename = "network.pdf",
                width = 10, 
                height = 10
        ) 
        standardGeneric("plotNetwork")
)

setGeneric(
        "selectRank", 
        function(
                object,
                vertexSize = 8, 
                vertexLabelCex = 0.75, 
                edgeAlpha = 0.2,
                filename = "Coverage.pdf",
                width = 10, 
                height = 10
        ) 
        standardGeneric("selectRank")
)

setGeneric(
        "plotCommNetwork", 
        function(
                object, 
                vertexInfo, 
                filename = "community.pdf",
                width = 10, 
                height = 10
        ) 
        standardGeneric("plotCommNetwork")
)

setGeneric(
        "statComm", 
        function(
                object, 
                filename = "stats.pdf",
                width = 10, 
                height = 10
        ) 
        standardGeneric("statComm")
)
