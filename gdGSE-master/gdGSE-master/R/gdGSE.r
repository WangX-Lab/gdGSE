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
    cell_type = unique(condition$CellType)

    result = c()
    for (p in 1:length(cell_type)) {
      celltype <- cell_type[[p]]

      aimCell <- condition[condition$CellType == celltype,]
      aimcell_count <- exp_matrix[, colnames(exp_matrix) %in% aimCell$CellID]

      Another <- condition[condition$CellType != celltype,]
      another_count <- exp_matrix[, colnames(exp_matrix) %in% Another$CellID]

      normal_rna_seq <- another_count
      row_mean <- rowMeans(normal_rna_seq, na.rm = TRUE)
      SD <- apply(normal_rna_seq, 1, sd, na.rm = TRUE)

      tumor_rna_seq2 <- aimcell_count
      rpk <- sweep(tumor_rna_seq2, 1, row_mean, "-")
      y3 <- ifelse(is.na(rpk), NA, ifelse(rpk > SD, 1, 0))

      y1 <- sapply(1:ncol(Signature), function(i) {
        geneset <- na.omit(Signature[, i])
        geneset_exp <- y3[rownames(y3) %in% geneset, , drop = FALSE]
        apply(geneset_exp, 2, function(tab_exp1) {
          sum(tab_exp1 == 1) / nrow(geneset_exp)
        })
      })

      colnames(y1) <- colnames(Signature)
      result = rbind(result, y1)
    }

    return(result)  # Return the result for single-cell data processing
  } else {

    exp_matrix <- as.data.frame(exp_matrix)
    condition <- as.data.frame(condition)

    pos_nor <- which(colnames(exp_matrix) %in% condition[which(condition[, 2] == "Normal"), 1])
    pos_dis <- which(colnames(exp_matrix) %in% condition[which(condition[, 2] == "Tumor"), 1])
    exp_nor <- exp_matrix[, pos_nor]
    exp_dis <- exp_matrix[, pos_dis]

    # Calculate row means for normal RNA sequence
    row_mean = rowMeans(exp_nor, na.rm = TRUE)
    row_mean = as.data.frame(row_mean)
    GeneID = rownames(row_mean)
    row_mean1 = cbind(GeneID, row_mean)

    # Calculate standard deviation for normal RNA sequence
    SD = apply(exp_nor, 1, function(x) sd(x, na.rm = TRUE))
    SD = as.data.frame(SD)
    GeneID = rownames(SD)
    SD1 = cbind(GeneID, SD)
    colnames(SD1) = c("GeneID", "nor_SD")

    # Prepare tumor RNA sequence
    GeneID = rownames(exp_dis)
    GeneID = as.data.frame(GeneID)
    exp_dis1 = cbind(GeneID, exp_dis)
    normal_mean = row_mean1

    # Merge normal and tumor data
    cordata = merge(normal_mean, exp_dis1, by = "GeneID")
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
      gene_set <- gene_set[gene_set %in% rownames(exp_dis)]
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
    return(y1)
  }
  # Check if it's single-cell data
}

# exp_matrix = read.csv("D:/E/Pathway_project/Rpackage/TCGA-CHOL_exp.csv", stringsAsFactors = FALSE, check.names = FALSE,row.names = 1)
# condition = read.csv("D:/E/Pathway_project/Rpackage/TCGA-CHOL_state.csv", stringsAsFactors = FALSE, check.names = FALSE)
# Signature = read.csv("D:/E/Pathway_project/Rpackage/Rpackage_pathway.csv", stringsAsFactors = FALSE, check.names = FALSE)

# gdGSE_scores <- gdGSE(exp_matrix, condition, Signature)
# usethis::use_data(exp_matrix, condition, Signature)

