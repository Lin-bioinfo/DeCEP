#' @title Data preprocessing and cell identity characterization (DeCEP state) of spatial transcriptomics (ST) data
#'
#' @description
#' \code{pre_st} performs data preprocessing of ST data across different spatial contexts.
#'
#' @param x A data frame or a Seurat object containing spatial gene expression data (raw count). For a data frame, each row represents a gene, and each column represents a spatial location.
#' @param coords A data frame with two columns (the column names are "x" and "y", corresponding to the two-dimensional coordinates of spatial locations).
#' @param genelist A data frame containing one column (the column name is "gene"). It records the functional gene list that reflects a specific gene program.
#' @param ref A reference scRNA-seq dataset. A data frame or a Seurat object containing single-cell gene expression data (raw count). For a data frame, each row represents a gene, and each column represents a cell.
#' @param cell_state A data frame containing the DeCEP states of a specific gene program corresponding to individual cells in the reference scRNA-seq dataset.
#' @param core_number The number of cores to use for parallel processing. If set to 1, no parallel processing is used.
#' @return A data frame with two columns (the column names are "spot_ID" and "DeCEP_state").
#' @import Seurat
#' @import spacexr
#' @import classInt
#' @export

pre_st <- function(x, coords, genelist, ref, cell_state, core_number = 2) {

	if (class(x) == "Seurat") {
		x <- x@assays$Spatial@counts
	}
	if (class(ref) == "Seurat") {
		ref <- ref@assays$RNA@counts
	}

	genes <- genelist$gene
	if (sum(rownames(x) %in% genes) < 0.7 * length(genes)) {
		stop("The number of functional genes is not enough to characterize the gene program!")
	}

	nUMI_st <- colSums(x) 
	puck <- SpatialRNA(coords, x, nUMI_st)

	ref_state <- as.factor(cell_state$DeCEP_state)
	names(ref_state) <- cell_state$cell_ID
	nUMI_sc <- colSums(ref) 
	reference <- Reference(ref, ref_state, nUMI_sc)

	myRCTD <- create.RCTD(puck, reference, max_cores = core_number, counts_MIN = 0, UMI_min = 0)
	myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

	weights <- myRCTD@results$weights
	norm_weights <- normalize_weights(weights)
	dat <- as.data.frame(as.matrix(norm_weights))

	Value <- as.numeric(dat$High)
	Intervals <- classIntervals(var = Value, n = 3, style = 'fisher', warnLargeN = FALSE)
	pos <- Intervals$brks
	dat$High_state <- 0
	dat[dat$High >= pos[3] & dat$High > 0.4, "High_state"] = 1

	Value <- as.numeric(dat$Low)
	Intervals <- classIntervals(var = Value, n = 3, style = 'fisher', warnLargeN = FALSE)
	pos <- Intervals$brks
	dat$Low_state <- 0
	dat[dat$Low >= pos[3] & dat$Low > 0.4, "Low_state"] = 1

	Value <- as.numeric(dat$Medium)
	Intervals <- classIntervals(var = Value, n = 3, style = 'fisher', warnLargeN = FALSE)
	pos <- Intervals$brks
	dat$Medium_state <- 0
	dat[dat$Medium >= pos[3] & dat$Medium > 0.4, "Medium_state"] <- 1

	dat$state <- ""
	dat[dat[,c("High_state")] == 1, "state"] <- "High"
	dat[dat[,c("Low_state")] == 1, "state"] <- "Low"
	dat[dat[,c("Medium_state")] == 1, "state"] <- "Medium"
	dat[apply(dat[ ,c("High_state", "Medium_state", "Low_state")], 1, sum) == 0, "state"] <- "Mixed"
	dat[apply(dat[ ,c("High_state", "Medium_state", "Low_state")], 1, sum) == 2, "state"] <- "Mixed"
	dat[apply(dat[ ,c("High_state", "Medium_state", "Low_state")], 1, sum) == 3, "state"] <- "Mixed"
	DeCEP_state <- data.frame(spot_ID = rownames(dat), DeCEP_state = dat$state)
	return(DeCEP_state)

}
