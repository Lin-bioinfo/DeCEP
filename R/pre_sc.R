#' @title Data preprocessing of single-cell RNA sequencing (scRNA-seq) data
#'
#' @description
#' \code{pre_sc} performs data preprocessing of scRNA-seq data across different cellular contexts.
#'
#' @param x A data frame or a Seurat object containing single-cell gene expression data (library size normalization and log2 transformation). For a data frame, each row represents a cell, and each column represents a gene.
#' @param label A data frame containing two columns (the column names are "cell_ID" and "cell_label"). It records the cellular context labels of individual cells.
#' @param modes The running mode of CoGPS, "discrete" or "continuous".
#' @param transition_point Numerical. The transition point of a specific gene program. It needs to be specified when modes = "continuous".
#' @param genelist A data frame containing one column (the column name is "gene"). It records the functional gene list that reflects a specific gene program.
#' @return Processed data. If modes = "discrete" it returns a list containing a series of data frames, if modes = "continuous" it returns a data frame via the time-lag method.
#' @import Seurat
#' @export

pre_sc <- function(x, label, modes = c("discrete", "continuous"), transition_point = NULL, genelist) {

	modes <- match.arg(modes)

	if (class(x) == "Seurat") {
		x <- x@assays$RNA@data
		x <- data.frame(t(as.matrix(x)))
	}
	cellID <- intersect(rownames(x), label$cell_ID)
	x <- x[match(cellID, rownames(x)), ]
	label <- label[match(cellID, label$cell_ID), ]
	genes <- genelist$gene

	if (sum(colnames(x) %in% genes) < 0.7 * length(genes)) {
		stop("The number of functional genes is not enough to characterize the gene program!")
	}

	if (modes == "discrete") {
		x_sc <- list()
		contexts <- unique(label$cell_label)
		for (i in 1:length(contexts)) {
			x_sc[[i]] <- x[rownames(x) %in% label[label$cell_label == contexts[i], "cell_ID"], ]
			x_sc[[i]] <- x_sc[[i]][ ,match(genes, colnames(x_sc[[i]]))[!is.na(match(genes, colnames(x_sc[[i]])))]]
		}
		names(x_sc) <- contexts
		return(x_sc)
	}

	if (modes == "continuous") {
		label <- label[order(label$cell_label, decreasing = FALSE), ]
		x <- x[match(label$cell_ID, rownames(x)), match(genes, colnames(x))[!is.na(match(genes, colnames(x)))]]
		pos <- sum(label$cell_label < transition_point)
		x_sc <- cbind(x[-c(1:pos), ], x[1:(dim(x)[1] - pos), ])
		colnames(x_sc)[1:(dim(x)[2])] <- paste0(colnames(x_sc)[1:(dim(x)[2])], "_T2")
		colnames(x_sc)[(dim(x)[2] + 1):dim(x_sc)[2]] <- paste0(colnames(x_sc)[(dim(x)[2] + 1):dim(x_sc)[2]], "_T1")
		x_sc <- as.data.frame(x_sc)
		return(x_sc)
	}

}
