#' @title Cell identity characterization (CoGPS score) of ST data
#'
#' @description
#' \code{ident_st_score} quantifies the activity of context-specific gene program depend on \link{net_st} for individual spatial locations.
#'
#' @param x_impute A data frame containing spatial gene expression data (after imputation) in a distinct spatial context. For a data frame, each row represents a spatial location, and each column represents a gene.
#' @param involvedID A character vector containing the names of spatial locations in the region of interest (ROI) and corresponding neighborhood.
#' @param net A context-dependent functional network obtained by \link{net_st}.
#' @return A data frame with two columns (the column names are "spot_ID" and "CoGPS_score").
#' @export

ident_st_score <- function(x_impute, involvedID, net) {

	x_impute <- data.frame(x_impute[match(net$Node$Hub_gene, rownames(x_impute)), match(involvedID, colnames(x_impute))])
	CoGPS_score <- data.frame(spot_ID = colnames(x_impute), CoGPS_score = colMeans(x_impute * net$Node$Weight))
	rownames(CoGPS_score) = NULL
	return(CoGPS_score)

}
