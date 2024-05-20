#' @title Cell identity characterization (CoGPS score) of scRNA-seq data
#'
#' @description
#' \code{ident_sc_score} quantifies the activity of context-specific gene program depend on \link{net_sc} for individual cells.
#'
#' @param x A data frame or a Seurat object containing single-cell gene expression data (library size normalization and log2 transformation) in a distinct cellular context. For a data frame, each row represents a cell, and each column represents a gene.
#' @param net A context-dependent functional network obtained by \link{net_sc}.
#' @param seed Set a random seed. Setting NULL will not set a seed.
#' @return A data frame with two columns (the column names are "cell_ID" and "CoGPS_score").
#' @import Seurat
#' @import Rmagic
#' @export

ident_sc_score <- function(x, net, seed = 1) {

	if (class(x) == "Seurat") {
		x <- x@assays$RNA@data
		x <- data.frame(t(as.matrix(x)))
	}

	x_smooth <- magic(as(as.matrix(x), "dgCMatrix"), genes = net$Node$Hub_gene, t = 2, seed = seed)$result
	x_smooth <- t(x_smooth)
	x_smooth <- data.frame(x_smooth[match(net$Node$Hub_gene, rownames(x_smooth)),])
	CoGPS_score <- data.frame(cell_ID = colnames(x_smooth), CoGPS_score = colMeans(x_smooth * net$Node$Weight))
	rownames(CoGPS_score) = NULL
	return(CoGPS_score)

}



#' @title Cell identity characterization (CoGPS state) of scRNA-seq data
#'
#' @description
#' \code{ident_sc_state} assigns the activity of context-specific gene program depend on \link{ident_sc_score} for individual cells.
#'
#' @param score A data frame with two columns (the column names are "cell_ID" and "CoGPS_score") obtained by \link{ident_sc_score}.
#' @return A data frame with three columns (the column names are "cell_ID", "CoGPS_score", and "CoGPS_state").
#' @import classInt
#' @export

ident_sc_state <- function(score) {
	Value = as.numeric(score$CoGPS_score)
	Intervals = classIntervals(var = Value, n = 3, style = 'fisher', warnLargeN = FALSE)
	pos = Intervals$brks
	state = data.frame(score, CoGPS_state = "")
	state[state$CoGPS_score >= pos[3], "CoGPS_state"] = "High"
	state[state$CoGPS_score >= pos[2] & state$CoGPS_score < pos[3], "CoGPS_state"] = "Medium"
	state[state$CoGPS_score < pos[2], "CoGPS_state"] = "Low"
	return(state)

}
