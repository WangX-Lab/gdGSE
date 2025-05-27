# gdGSE

The gdGSE algorithm processes either bulk or single-cell transcriptomic data, and outputs a gene set enrichment matrix by discretizing the original gene expression matrix. 

&nbsp;

# Description

This vignette shows an example of how to use the gdGSE algorithm to calculate the gene signatures or pathways enrichment score in R. gdGSE is an algorithm to enrichment score at bulk or single-cell transcriptomic data. 



# Details

+ The function gdGSE()` is used to calculate the gene signatures or pathways enrichment score.
  + **exp_matrix:** gene expression profiles in tumor and/or normal samples (Bulk RNA-Seq or Single-cell  RNA-Seq). **Note**: when reading the input file, please set "row.names=1. **condition**: sample type ("Tumor" or "Normal") in bulk or cell annotation of single cells in "exp_matrix". **Note**: when in bulk, the "State" column is "Normal" and "Tumor". **Signature**: A list  of signatures, for example: the Gene Ontology Biological Process (GOBP) gene sets, cell type-specific signatures, or a gene list correlated with interested features. **Note**: The gene name type in the Signature must be consistent with the gene type in the exp_matrix (for example, if the gene type in the exp_matrix is "Entrez Gene ID", then the gene name type in the Signature should also be "Entrez Gene ID"). 
  + "data_type" is a character vector indicating the type of data (Bulk RNA-Seq or Single-cell  RNA-Seq). 

&nbsp;&nbsp;

# Installation

- You can install the released version of **gdGSE** with:
  &nbsp;

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("WangX-Lab/gdGSE")
```

&nbsp;
&nbsp;

# Examples

### **Example for "Bulk RNA-Seq"**

```R
library(gdGSE)
file_path <- system.file("extdata", "Bulk_example.RData", package = "gdGSE")
load(file_path)
ls()
#[1] "Bulk_condition"  "Bulk_exp_matrix" "Bulk_Signature"  "file_path" 

gdGSE_Score <- gdGSE(Bulk_exp_matrix, Bulk_condition, Bulk_Signature,data_type = "Bulk") #If Users calculate enrichment score in single-cell, please change "Bulk" into "SingleCell".
```

**Bulk_exp_matrix**

```R
Bulk_exp_matrix[1:5,1:5]
```

| row.names | TCGA-3X-AAV9-01A-72R-A41I-07 | TCGA-3X-AAVA-01A-11R-A41I-07 | TCGA-3X-AAVB-01A-31R-A41I-07 | TCGA-3X-AAVC-01A-21R-A41I-07 | TCGA-3X-AAVE-01A-11R-A41I-07 |
| :-------: | :--------------------------: | :--------------------------: | :--------------------------: | :--------------------------: | :--------------------------: |
| 100130426 |              0               |              0               |              0               |              0               |              0               |
| 100133144 |           1.68997            |           2.343863           |           1.271963           |              0               |           3.066434           |
| 100134869 |           4.506246           |           4.537551           |           2.438186           |           2.438186           |           3.105058           |
|   10357   |           6.956025           |           6.373866           |           7.202411           |           7.924523           |           6.919972           |
|   10431   |           9.750945           |           9.955643           |           9.092849           |           10.09296           |           9.827594           |

**Bulk_condition**

```R
Bulk_condition[1:45,]
```

|            sample            |  State   |
| :--------------------------: | :------: |
| TCGA-3X-AAV9-01A-72R-A41I-07 |  Tumor   |
| TCGA-3X-AAVA-01A-11R-A41I-07 | 14.16773 |
|             ....             |   ....   |
| TCGA-W5-AA34-11A-11R-A41I-07 |  Normal  |
| TCGA-ZU-A8S4-11A-11R-A41I-07 |  Normal  |

**Bulk_signature**

```R
Bulk_Signature[1:5,1:5]
```

| GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE | GOBP_REPRODUCTION | GOBP_SINGLE_STRAND_BREAK_REPAIR | GOBP_REGULATION_OF_DNA_RECOMBINATION | GOBP_REGULATION_OF_MITOTIC_RECOMBINATION |
| :-----------------------------------: | :---------------: | :-----------------------------: | :----------------------------------: | :--------------------------------------: |
|                 10000                 |        100        |            100133315            |                10039                 |                  10111                   |
|                 10891                 |       10007       |              1161               |                10097                 |                  126549                  |
|                 11232                 |     100125288     |               142               |                10111                 |                  201516                  |
|                  142                  |     100130958     |             200558              |                10189                 |                   2068                   |
|                 1763                  |     100130988     |              2074               |                10459                 |                   4292                   |

&nbsp;

### **Example for "Single cell RNA-Seq"**

```
## Example for "SingleCell"
file_path <- system.file("extdata", "SingleCell_example.RData", package = "gdGSE")
load(file_path)
ls()
# [1] "SC_condition"  "SC_exp_matrix" "SC_Signature"  "file_path" 
gdGSE_Score <- gdGSE(SC_exp_matrix, SC_condition, SC_Signature,data_type = "SingleCell")  #If Users calculate enrichment score in Bulk, please change "SingleCell" into "Bulk".
```

**SC_exp_matrix**

```R
SC_exp_matrix[1:5,1:5]
```

|   row.names   | HCC08T_TCTGGAACACAGACAG | HCC02T_ATTACTCAGCGCCTTG | HCC08T_CCACGGAGTGTGCGTC | HCC07T_AGTTGGTCACAAGTAA | HCC06T_CACAGGCTCTTTAGTC |
| :-----------: | :---------------------: | :---------------------: | :---------------------: | :---------------------: | :---------------------: |
| RP11-34P13.7  |            0            |            0            |            0            |            0            |            0            |
|  FO538757.2   |            1            |            1            |            0            |            0            |            0            |
|  AP006222.2   |            0            |            0            |            1            |            0            |            0            |
| RP4-669L17.10 |            0            |            0            |            0            |            0            |            0            |
| RP5-857K21.4  |            0            |            0            |            0            |            0            |            0            |

**SC_condition**

```
SC_condition[1:5,]
```

|         CellID          |  CellType   |
| :---------------------: | :---------: |
| HCC01T_AAACGGGTCTGACCTC |   Myeloid   |
| HCC01T_AAATGCCGTATAGGGC |      T      |
| HCC01T_AACCGCGCACGGATAG |   Myeloid   |
| HCC01T_AAGGTTCGTCATTAGC |      T      |
| HCC01T_AAGTCTGAGATGTGGC | Endothelial |

**SC_Signature**

```
SC_Signature[1:5,1:5]
```

| Hepatocyte | Myeloid |   B   | Endothelial | Fibroblast |
| :--------: | :-----: | :---: | :---------: | :--------: |
|   UBE2C    |  MARCO  | IGHG4 |    FABP4    |   COL1A1   |
|    TUBB    |  CD5L   | IGLC2 |    RBP7     |    LUM     |
|   VCX3A    |  C1QB   | IGHG1 |   IGFBP7    |    DCN     |
|  AKR1B10   |  C1QA   | IGKC  |   TM4SF1    |   COL3A1   |
|    NQO1    | SLC40A1 | IGLC3 |    CLDN5    |   COL1A2   |


**Calculate empirical p-values**
```R
library(tidyverse)
set.seed(123)

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

p_values <- gdGSE1(exp_matrix, condition, gene_sets, data_type = "Bulk", n_perm = 1000) # Larger values of n_perm may significantly increase computation time.
```

# Contact

E-mail any questions to Xiaosheng Wang (xiaosheng.wang@cpu.edu.cn)
