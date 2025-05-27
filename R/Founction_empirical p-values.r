# After installing the gdGSE R package first, then run the following functions
# permutations
#' gdGSE: a novel pathway enrichment evaluation algorithm based on discretization of gene expression values
#'
#' @param exp_matrix The expression matrix of the samples or single cells for which enrichment scoring needs to be performed
#' @param condition The conditions of bulk "Tumor" or "Normal" or cell annotation of single-cell transcriptomes.
#' @param Signature A list of signatures,for example:the Gene Ontology Biological Process(GOBP) gene sets, cell type-specific signatures or a gene list correlated with interested features.
#' @param data_type "Bulk" or "SingleCell".
#' @return Description of the return value.
#' @export
#'
#' @importFrom stats na.omit sd
#' @export gdGSE
#'
gdGSE <- function(exp_matrix, condition, Signature, data_type) {
  if (missing(data_type)) {
    stop("The data_type is lack. The user must select Bulk or SingleCell")
  }
  if (data_type == "SingleCell") {
    cell_type = unique(condition[,2])
    result = c()
    for (p in 1:length(cell_type)) {
      celltype <- cell_type[[p]]
      aimCell <- condition[condition[,2] == celltype,]
      aimcell_matrix <- exp_matrix[, colnames(exp_matrix) %in% aimCell[,1]]

      Another <- condition[condition[,2] != celltype,]
      another_matrix <- exp_matrix[, colnames(exp_matrix) %in% Another[,1]]

      normal_rna_seq <- another_matrix
      row_mean = rowMeans(another_matrix, na.rm = TRUE)
      row_mean = as.data.frame(row_mean)
      GeneID = rownames(row_mean)
      row_mean1 = cbind(GeneID, row_mean)

      SD = apply(another_matrix, 1, function(x) sd(x, na.rm = TRUE))
      SD = as.data.frame(SD)
      GeneID = rownames(SD)
      SD1 = cbind(GeneID, SD)
      colnames(SD1) = c("GeneID", "nor_SD")

      # Prepare tumor RNA sequence
      GeneID = rownames(aimcell_matrix)
      GeneID = as.data.frame(GeneID)
      aimcell_matrix1 = cbind(GeneID, aimcell_matrix)
      normal_mean = row_mean1

      # Merge normal and tumor data
      cordata = merge(normal_mean, aimcell_matrix1, by = "GeneID")
      kb <- cordata$row_mean
      seq_data <- cordata[,-c(1:2)]
      rpk <- seq_data - kb
      rownames(rpk) = cordata$GeneID

      # Create binary matrix based on standard deviation
      SD = SD1$nor_SD
      y = rpk
      y3 <- ifelse(is.na(y), NA, ifelse(y > SD, 1, 0))
      
      # Genes in pathway > 0
      Signature <- lapply(Signature, function(gene_set) {
      gene_set <- gene_set[gene_set %in% rownames(aimcell_matrix)]
        return(gene_set)
      })
      Signature <- Signature[sapply(Signature, length) > 1]
      data <- Signature
      genesets <- as.data.frame(sapply(data, "[", i = 1:max(sapply(data, length))))

      # Calculate ratio for each gene set
      T_exp_dis = y3
      y1 <- matrix(0, nrow = dim(genesets)[2], ncol = dim(T_exp_dis)[2])  # Pre-allocate result storage space
      for (i in 1:dim(genesets)[2]) {
        geneset <- genesets[, i]
        geneset1 <- na.omit(geneset)
        geneset_exp <- T_exp_dis[rownames(T_exp_dis) %in% geneset1, ]
        geneset_exp1 <- na.omit(geneset_exp)

        for (j in 1:dim(geneset_exp1)[2]) {
          tab_exp1 <- geneset_exp1[, j]
          rat <- length(which(tab_exp1 == 1)) / length(geneset_exp1[, 1])
          y1[i, j] <- rat
        }
      }
      colnames(y1) <- colnames(geneset_exp1)
      # Check if it's single-cell data
      result = cbind(result, y1)
    }
   rownames(result) = colnames(genesets)
    return(result)  # Return the result for single-cell data processing
  } else {
    sample_type = unique(condition[,2])
    result = c()
    for (p in 1:length(sample_type)) {
      sample <- sample_type[[p]]
      aimSample <- condition[condition[,2] == sample,]
      aimSample_matrix <- exp_matrix[, colnames(exp_matrix) %in% aimSample[,1]]

      Another <- condition[condition[,2] != sample,]
      another_matrix <- exp_matrix[, colnames(exp_matrix) %in% Another[,1]]

      normal_rna_seq <- another_matrix
      row_mean = rowMeans(another_matrix, na.rm = TRUE)
      row_mean = as.data.frame(row_mean)
      GeneID = rownames(row_mean)
      row_mean1 = cbind(GeneID, row_mean)

      SD = apply(another_matrix, 1, function(x) sd(x, na.rm = TRUE))
      SD = as.data.frame(SD)
      GeneID = rownames(SD)
      SD1 = cbind(GeneID, SD)
      colnames(SD1) = c("GeneID", "nor_SD")

      # Prepare tumor RNA sequence
      GeneID = rownames(aimSample_matrix)
      GeneID = as.data.frame(GeneID)
      aimSample_matrix1 = cbind(GeneID, aimSample_matrix)
      normal_mean = row_mean1

      # Merge normal and tumor data
      cordata = merge(normal_mean, aimSample_matrix1, by = "GeneID")
      kb <- cordata$row_mean
      seq_data <- cordata[,-c(1:2)]
      rpk <- seq_data - kb
      rownames(rpk) = cordata$GeneID

      # Create binary matrix based on standard deviation
      SD = SD1$nor_SD
      y = rpk
      y3 <- ifelse(is.na(y), NA, ifelse(y > SD, 1, 0))

      # Genes in pathway > 0
      Signature <- lapply(Signature, function(gene_set) {
        gene_set <- gene_set[gene_set %in% rownames(aimSample_matrix)]
        return(gene_set)
      })
      Signature <- Signature[sapply(Signature, length) > 1]
      data <- Signature
      genesets <- as.data.frame(sapply(data, "[", i = 1:max(sapply(data, length))))

      # Calculate ratio for each gene set
      T_exp_dis = y3
      y1 <- matrix(0, nrow = dim(genesets)[2], ncol = dim(T_exp_dis)[2])  # Pre-allocate result storage space
      for (i in 1:dim(genesets)[2]) {
        geneset <- genesets[, i]
        geneset1 <- na.omit(geneset)
        geneset_exp <- T_exp_dis[rownames(T_exp_dis) %in% geneset1, ]
        geneset_exp1 <- na.omit(geneset_exp)

        for (j in 1:dim(geneset_exp1)[2]) {
          tab_exp1 <- geneset_exp1[, j]
          rat <- length(which(tab_exp1 == 1)) / length(geneset_exp1[, 1])
          y1[i, j] <- rat
        }
      }
      colnames(y1) <- colnames(geneset_exp1)
      rownames(y1) = colnames(genesets)
      result = cbind(result,y1)   
  # Check if it's single-cell data
  }
  return(result)
 }
}

#--------------------------------#
gdGSE1 <- function(exp_matrix, condition, gene_sets, data_type = "Bulk", n_perm = 1000) {
  
    # Step 1: Compute pathway enrichment scores across samples
    observed_es <- gdGSE(exp_matrix, condition, gene_sets, data_type = "Bulk")
    
    # Step 2:  Initialize null distributions as 3D arrays [pathways × samples × permutations]
    null_es <- array(NA, dim = c(nrow(observed_es), ncol(observed_es), n_perm),
                     dimnames = list(rownames(observed_es), colnames(observed_es), NULL))
    
    # Step3: Permutation
    for (k in seq_len(n_perm)) {
        perm_labels <- sample(condition[,2])                          # Shuffle labels
        perm_condition <- data.frame(sample = colnames(exp_matrix),
                                     State = perm_labels)
        null_es[,,k] <- gdGSE(exp_matrix, perm_condition, gene_sets, data_type = "Bulk") # Store full matrix
    }
    
    # Step4: Calculate p-values for each pathway 
    p_values <- sapply(seq_len(nrow(observed_es)), function(i) {
        obs_mean <- mean(observed_es[i, ])                    
        null_mean <- apply(null_es[i,,], MARGIN=2, FUN=mean)      
        (sum(null_mean >= obs_mean) +1 ) / (n_perm +1 )         # Empirical p-value
    })
    
    names(p_values) <- rownames(observed_es)
    return(p_values)
}

# run gdGSE1 founction
library(tidyverse)
set.seed(123)
p_values <- gdGSE1(exp_matrix, condition, gene_sets, data_type = "Bulk", n_perm = 1000) # Larger values of n_perm may significantly increase computation time.
