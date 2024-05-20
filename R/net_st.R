#' @title Preparation for context-dependent functional network construction of ST data
#'
#' @description
#' \code{net_st_pre} constructs the functional network and pinpoints hub genes that reflect a specific gene program in a distinct spatial context.
#'
#' @param x A data frame containing spatial gene expression data (library size normalization and log2 transformation). For a data frame, each row represents a spatial location, and each column represents a gene.
#' @param coords A data frame with two columns (the column names are "x" and "y", corresponding to the two-dimensional coordinates of spatial locations).
#' @param genelist A data frame containing one column (the column name is "gene"). It records the functional gene list that reflects a specific gene program.
#' @param roiID A character vector containing the names of spatial locations in the region of interest (ROI).
#' @param r The distance used to identify the neighborhood of the user-specified ROI.
#' @return Processed data. A list containing two elements ("dat" and "neighborID"). The "dat" element contains a data frame with paired data comprising actual- and pseudo-spots. The "neighborID" element contains a character vector with the names of spatial locations in the corresponding neighborhood of the ROI.
#' @import bnlearn
#' @import parallel
#' @export

net_st_pre <- function(x, coords, genelist, roiID, r) {

	genes <- genelist$gene

	dis <- as.matrix(dist(coords))
	dis_list <- split(dis, rep(1:ncol(dis), each = nrow(dis)))
	names(dis_list) <- colnames(dis)
	roi_neighbor <- lapply(dis_list, function(x) {
		samp0 = rownames(dis)[x < r & x > 0]
		samp = samp0[!samp0 %in% roiID]
		return(samp)
	})
	names(roi_neighbor) <- rownames(dis)

	x1 <- x[ ,match(genes, colnames(x))[!is.na(match(genes, colnames(x)))]]
	roi_neighbor <- roi_neighbor[match(roiID, names(roi_neighbor))]
	roi_neighbor <- roi_neighbor[as.numeric(which(unlist(lapply(roi_neighbor, function(x) {(length(x))}[[1]][1])) != 0))]
	x1_neighbor <- lapply(roi_neighbor, function(x) {colMeans(x1[x, , drop = FALSE])})
	x1_neighbor <- do.call(rbind, x1_neighbor)
	x1 <- x1[names(roi_neighbor), ]
	colnames(x1_neighbor) <- paste0(colnames(x1_neighbor), "_neighbor")
	x_st <- as.data.frame(cbind(x1_neighbor, x1))
	neighborID <- unique(do.call(c, roi_neighbor))
	neighborID <- neighborID[!neighborID %in% roiID]
	x_st <- list(dat = x_st, neighborID = neighborID)
	return(x_st)

}


#' @title Context-dependent functional network construction of ST data
#'
#' @description
#' \code{net_st} constructs the functional network and pinpoints hub genes that reflect a specific gene program in a distinct spatial context.
#'
#' @param x_st A data frame containing paired gene expression data after \link{net_st_pre} data preprocessing.
#' @param genelist A data frame containing one column (the column name is "gene"). It records the functional gene list that reflects a specific gene program.
#' @param select A numerical value between 0 and 1 that specifies the top percentage of important spatial locations involved in network construction.
#' @param cluster The number of cores to use for parallel processing. If set to 1, no parallel processing is used.
#' @param seed Set a random seed. Setting NULL will not set a seed.
#' @return A list with two elements ("Edge" and "Node"). The "Edge" element containing edge attributes represents transcriptional co-regulation of functional genes. The "Node" element containing node attributes represents the weight of functional genes.
#' @import bnlearn
#' @import parallel
#' @export

net_st <- function(x_st, genelist, select = 0.2, cluster = NULL, seed = NULL) {

	genes <- genelist$gene
	select <- 1 - select

	x_sd <- apply(x_st, 1, sd)
	x_mean <- apply(x_st, 1, mean)
	pos_sd <- quantile(x_sd, probs = select)
	pos_mean <- quantile(x_mean, probs = select)
	x_select <- x_st[x_sd >= pos_sd & x_mean >= pos_mean, ]
	x_select <- x_select[rowSums(x_select) != 0, ]
	gene_filterd <- colnames(x_select)[colSums(x_select) == 0]
	gene_filterd <- unique(unlist(lapply(gene_filterd, function(x) {strsplit(x, "_", fixed = TRUE)[[1]][1]})))
	x_select <- x_select[ ,!colnames(x_select) %in% c(paste0(gene_filterd, "_neighbor"), gene_filterd)]

	if (dim(x_select)[2] * 0.5 < length(genes) * 0.5) {
		stop("The number of functional genes is not enough to characterize the gene program!")
	}

	gene_neighbor <- names(x_select)[grep("_neighbor", names(x_select))]
	gene <- names(x_select)[grep("_neighbor", names(x_select), invert = TRUE)]
	from <- as.character(sapply(gene, rep, times = length(gene)))
	to = rep(gene_neighbor, times = length(gene))
	blacklist_inter <- cbind(from, to)
	from_intra <- c(as.character(sapply(gene, rep, times = length(gene))), as.character(sapply(gene_neighbor, rep, times = length(gene_neighbor))))
	to_intra <- c(rep(gene, times = length(gene)), rep(gene_neighbor, times = length(gene_neighbor)))
	blacklist_intra <- cbind(from_intra, to_intra)
	blacklist <- rbind(blacklist_inter, blacklist_intra)

	x_select <- as.data.frame(scale(x_select))

	if (!is.null(cluster)) {
		clusterSetRNGStream(cluster, iseed = seed)
		net <- bn.boot(data = x_select, statistic = function(x) x, algorithm = "mmhc", algorithm.args = list(blacklist = blacklist), R = 100, m = nrow(x_select), cluster = cluster)
		stopCluster(cluster)
	} else {
		net <- bn.boot(data = x_select, statistic = function(x) x, algorithm = "mmhc", algorithm.args = list(blacklist = blacklist), R = 100, m = nrow(x_select))
	}

	directed_arcs <- lapply(net, directed.arcs)
	directed_arcs <- lapply(directed_arcs, function(x) {
		x = as.data.frame(x)
		x$Sign = paste0(x$from, "->", x$to)
		return(x)
	})

	directed_arcs <- do.call(rbind, directed_arcs)
	directed_arcs <- as.data.frame(table(directed_arcs$Sign))
	colnames(directed_arcs) <- c("Sign", "Freq")
	directed_arcs1 <- directed_arcs[directed_arcs$Freq > 20 * log10(sqrt(dim(x_select)[1])), ]
	directed_arcs1$from <- unlist(lapply(as.character(directed_arcs1$Sign), function(x) {strsplit(x, "->", fixed = TRUE)[[1]][1]}))
	directed_arcs1$to <- unlist(lapply(as.character(directed_arcs1$Sign), function(x) {strsplit(x, "->", fixed = TRUE)[[1]][2]}))

	arcs <- cbind(directed_arcs1, direction = "sign")
	arcs$Sign_rev <- paste0(arcs$to, "->", arcs$from)
	arcs$Sign_rev_match <- arcs$Sign_rev %in% arcs$Sign
	arcs$Freq_rev <- arcs[match(arcs$Sign_rev, arcs$Sign), "Freq"]
	arcs$minus_Freq <- (arcs$Freq - arcs$Freq_rev)/(arcs$Freq + arcs$Freq_rev)
	arcs$weight <- arcs$Freq / 100
	arcs$weight_rev <- arcs$Freq_rev / 100
	arcs$Sign_1 <- ""
	arcs$weight_1 <- ""
	arcs[is.na(arcs$minus_Freq), "Sign_1"] <- as.character(arcs[is.na(arcs$minus_Freq),]$Sign)
	arcs[is.na(arcs$minus_Freq), "weight_1"] <- arcs[is.na(arcs$minus_Freq),]$weight
	arcs[!is.na(arcs$minus_Freq) & arcs$minus_Freq > 0.1, "Sign_1"] <- as.character(arcs[!is.na(arcs$minus_Freq) & arcs$minus_Freq > 0.1,]$Sign)
	arcs[!is.na(arcs$minus_Freq) & arcs$minus_Freq > 0.1, "weight_1"] <- arcs[!is.na(arcs$minus_Freq) & arcs$minus_Freq > 0.1,]$weight
	arcs[!is.na(arcs$minus_Freq) & arcs$minus_Freq < -0.1, "Sign_1"] <- as.character(arcs[!is.na(arcs$minus_Freq) & arcs$minus_Freq < -0.1,]$Sign_rev)
	arcs[!is.na(arcs$minus_Freq) & arcs$minus_Freq < -0.1, "weight_1"] <- arcs[!is.na(arcs$minus_Freq) & arcs$minus_Freq < -0.1,]$weight_rev
	arcs[arcs$Sign_1 == "", "Sign_1"] <- unlist(lapply(as.character(arcs$Sign)[arcs$Sign_1 == ""], function(x) {
								gene = sort(c(strsplit(x, "->")[[1]][1], strsplit(x, "->")[[1]][2]))
								return(paste0(gene[1], "->", gene[2], " | ", gene[2], "->", gene[1]))
						 }))
	arcs[arcs$weight_1 == "", "weight_1"] <- (arcs[arcs$weight_1 == "", ]$weight + arcs[arcs$weight_1 == "", ]$weight_rev) / 2
	arcs <- arcs[!duplicated(arcs$Sign_1),]

	arcs_all <- arcs[c("Sign_1", "weight_1")]

	arcs_new0 = arcs_new1 = arcs_new2 = arcs_all[grep(" | ", arcs_all$Sign_1),]
	arcs_new1$Sign_1 <- unlist(lapply(arcs_new0$Sign_1, function(x) {strsplit(x, " | ")[[1]][1]}))
	arcs_new2$Sign_1 <- unlist(lapply(arcs_new0$Sign_1, function(x) {strsplit(x, " | ")[[1]][3]}))

	arcs_final <- rbind(arcs_all[!c(1:dim(arcs_all)[1]) %in% grep(" | ", arcs_all$Sign_1),], arcs_new1, arcs_new2)
	arcs_final$from <- unlist(lapply(arcs_final$Sign_1, function(x) {strsplit(x, "->")[[1]][1]}))
	arcs_final$to <- unlist(lapply(arcs_final$Sign_1, function(x) {strsplit(x, "->")[[1]][2]}))
	colnames(arcs_final) <- c("Edge", "Edge_strength", "From", "To")

	strength_score <- aggregate(as.numeric(arcs_final$Edge_strength), by = list(arcs_final$From), FUN = sum)
	colnames(strength_score) <- c("Hub_gene", "Weight")
	strength_score$Hub_gene <- gsub("_neighbor", "", strength_score$Hub_gene)
	net_result <- list(Edge = arcs_final, Node = strength_score)
	return(net_result)

}
