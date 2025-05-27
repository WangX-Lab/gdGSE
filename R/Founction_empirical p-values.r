# After installing the gdGSE R package first, then run the following functions
# permutations
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
# library(gdGSE) 
library(tidyverse)
set.seed(123)
p_values <- gdGSE1(exp_matrix, condition, gene_sets, data_type = "Bulk", n_perm = 1000) # Larger values of n_perm may significantly increase computation time.
