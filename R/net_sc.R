#' @title Context-dependent functional network construction of scRNA-seq data
#'
#' @description
#' \code{net_sc} constructs the functional network and pinpoints hub genes that reflect a specific gene program in a distinct cellular context.
#'
#' @param x A data frame containing single-cell gene expression data after \link{pre_sc} data preprocessing.
#' @param modes The running mode of DeCEP, "discrete" or "continuous".
#' @param genelist A data frame containing one column (the column name is "gene"). It records the functional gene list that reflects a specific gene program.
#' @param select A numerical value between 0 and 1 that specifies the top percentage of important cells involved in network construction.
#' @param cluster The number of cores to use for parallel processing. If set to 1, no parallel processing is used.
#' @param seed Set a random seed. Setting NULL will not set a seed.
#' @return A list with two elements ("Edge" and "Node"). The "Edge" element containing edge attributes represents transcriptional co-regulation of functional genes. The "Node" element containing node attributes represents the weight of functional genes.
#' @import bnlearn
#' @import parallel
#' @export

net_sc <- function(x, modes = c("discrete", "continuous"), genelist, select = 0.2, cluster = NULL, seed = NULL) {

	modes <- match.arg(modes)
	genes <- genelist$gene
	select <- 1 - select

	if (modes == "discrete") {
		x_sd <- apply(x, 1, sd)
		x_mean <- apply(x, 1, mean)
		pos_sd <- quantile(x_sd, probs = select)
		pos_mean <- quantile(x_mean, probs = select)
		x_select <- x[x_sd >= pos_sd & x_mean >= pos_mean, ]
		x_select <- x_select[rowSums(x_select) != 0, colSums(x_select) != 0]
		x_select <- x_select[, colSums(x_select) > dim(x_select)[1] * 0.01]

		if (dim(x_select)[2] < length(genes) * 0.5) {
			stop("The number of functional genes is not enough to characterize the gene program!")
		}
		x_select <- as.data.frame(scale(x_select))

		if (!is.null(cluster)) {
			clusterSetRNGStream(cluster, iseed = seed)
			net <- bn.boot(data = x_select, statistic = function(x) x, algorithm = "mmhc", R = 100, m = nrow(x_select), cluster = cluster)
			stopCluster(cluster)
		} else {
			net <- bn.boot(data = x_select, statistic = function(x) x, algorithm = "mmhc", R = 100, m = nrow(x_select))
		}
	}

	if (modes == "continuous") {
		x_sd <- apply(x, 1, sd)
		x_mean <- apply(x, 1, mean)
		pos_sd <- quantile(x_sd, probs = select)
		pos_mean <- quantile(x_mean, probs = select)
		x_select <- x[x_sd >= pos_sd & x_mean >= pos_mean, ]
		x_select <- x_select[rowSums(x_select) != 0, ]
		gene_filterd <- colnames(x_select)[colSums(x_select) == 0]
		gene_filterd <- unique(unlist(lapply(gene_filterd, function(x) {strsplit(x, "_", fixed = TRUE)[[1]][1]})))
		x_select <- x_select[ ,!colnames(x_select) %in% c(paste0(gene_filterd, "_T1"), paste0(gene_filterd, "_T2"))]

		if (dim(x_select)[2] * 0.5 < length(genes) * 0.5) {
			stop("The number of functional genes is not enough to characterize the gene program!")
		}

		gene_T1 <- names(x_select)[grep("_T1", names(x_select))]
		gene_T2 <- names(x_select)[grep("_T2", names(x_select))]
		from <- as.character(sapply(gene_T2, rep, times = length(gene_T2)))
		to <- rep(gene_T1, times = length(gene_T2))
		blacklist_inter <- cbind(from, to)
		from_intra <- c(as.character(sapply(gene_T2, rep, times = length(gene_T2))), as.character(sapply(gene_T1, rep, times = length(gene_T1))))
		to_intra <- c(rep(gene_T2, times = length(gene_T2)), rep(gene_T1, times = length(gene_T1)))
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
	if (modes == "continuous") {
		strength_score$Hub_gene <- gsub("_T1", "", strength_score$Hub_gene)
	}
	net_result <- list(Edge = arcs_final, Node = strength_score)
	return(net_result)

}
