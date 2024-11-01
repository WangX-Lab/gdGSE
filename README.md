# gdGSE

The gdGSE algorithm processes either bulk or single-cell transcriptomic data, and outputs a gene set enrichment matrix by discretizing the original gene expression matrix. 

&nbsp;

# Description

This vignette shows an example of how to use the gdGSE algorithm to calculate the gene signature or pathway enrichment score in R. gdGSE is an algorithm to evaluate gene set enrichment score in bulk or single-cell transcriptomic data. 

&nbsp;

# Details

+ ### Input: three files needed for the gdGSE function:

  1. **exp_matrix:** gene expression profiles in tumor and/or normal samples (Bulk RNA-Seq or Single cell  RNA-Seq). **Note**: when reading the input file, please set "row.names=1". 

     ```
                           Table 1. mRNA expression of input data(ENTREZID)
     ```

     |           | TCGA-3X-AAV9-01A-72R-A41I-07 | TCGA-3X-AAVA-01A-11R-A41I-07 | TCGA-3X-AAVB-01A-31R-A41I-07 | TCGA-3X-AAVC-01A-21R-A41I-07 | TCGA-3X-AAVE-01A-11R-A41I-07 |
     | :-------: | :--------------------------: | :--------------------------: | :--------------------------: | :--------------------------: | :--------------------------: |
     | 100130426 |           0.000000           |           0.000000           |           0.000000           |           0.000000           |           0.000000           |
     | 100133144 |           1.689970           |           2.343863           |           1.271963           |           0.000000           |           3.066434           |
     | 100134869 |           4.506246           |           4.537551           |           2.438186           |           2.267596           |           3.105058           |
     |   10357   |           6.956025           |           6.373866           |           7.202411           |           7.924523           |           6.919972           |

     

  2. **condition**: sample type ("Tumor" or "Normal") in bulk or cell annotation of single cells in "exp_matrix". **Note**: when in bulk, the "State" column are "Normal" and "Tumor".

     ```
                       Table2. Identification of tumor and normal samples
     ```

     |            sample            | State  |
     | :--------------------------: | :----: |
     | TCGA-3X-AAV9-01A-72R-A41I-07 | Tumor  |
     | TCGA-3X-AAVA-01A-11R-A41I-07 | Tumor  |
     | TCGA-3X-AAVB-01A-31R-A41I-07 | Tumor  |
     | TCGA-3X-AAVC-01A-21R-A41I-07 | Tumor  |
     | TCGA-3X-AAVE-01A-11R-A41I-07 | Tumor  |
     | TCGA-W5-AA2I-11A-11R-A41I-07 | Normal |
     | TCGA-W5-AA2Q-11A-11R-A41I-07 | Normal |

     

  3. **Signature**: A list  of signatures, for example: the Gene Ontology Biological Process (GOBP) gene sets, cell type-specific signatures or a gene list correlated with interested features . 

     ```
                      Table3. GOBP gene sets(ENTREZID)
     ```

     | GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE | GOBP_REPRODUCTION | GOBP_SINGLE_STRAND_BREAK_REPAIR |
     | :-----------------------------------: | :---------------: | :-----------------------------: |
     |                 10000                 |        100        |            100133315            |
     |                 10891                 |       10007       |              1161               |
     |                 11232                 |     100125288     |               142               |
     |                  142                  |     100130958     |             200558              |
     |                 .....                 |       .....       |              .....              |

     Example **"exp_matrix"**: gene expression profiles in cholangiocarcinoma from the genomic data commons data portal (<https://portal.gdc.cancer.gov/>). There are 45 samples (36 tumor and 9 normal samples) and 20,531 genes in **"exp"**. The gdGSE function will output the pathway enrichment score for each of the 36 tumor samples as shown in **Table 4**. 

     ```
              Table 4. pathway enrichment score of tumor samples in output data
     ```

     |                                          | TCGA-3X-AAV9-01A-72R-A41I-07 | TCGA-3X-AAVA-01A-11R-A41I-07 |
     | ---------------------------------------- | :--------------------------: | :--------------------------: |
     | GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE    |           0.419354           |          0.4193548           |
     | GOBP_REPRODUCTION                        |           0.527494           |          0.4249830           |
     | GOBP_SINGLE_STRAND_BREAK_REPAIR          |          0.3636364           |          0.3636364           |
     | GOBP_REGULATION_OF_DNA_RECOMBINATION     |          0.5789474           |          0.4812030           |
     | GOBP_REGULATION_OF_MITOTIC_RECOMBINATION |          0.4285714           |          0.4285714           |
     | GOBP_MITOTIC_SPINDLE_ELONGATION          |          0.9166667           |          0.8333333           |

     

# Installation

- Users can install the released version of **gdGSE** with:
  &nbsp;

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("WangX-Lab/gdGSE")
```

&nbsp;

# Examples

```
library(gdGSE)

## Example for "Bulk"
file_path <- system.file("extdata", "Bulk_example.RData", package = "gdGSE")
load(file_path)
ls()
# [1] "Bulk_condition"  "Bulk_exp_matrix" "Bulk_Signature"  "file_path" 
gdGSE_Score <- gdGSE(Bulk_exp_matrix,Bulk_condition, Bulk_Signature,data_type = "Bulk") #If Users calculate enrichment score in Bulk, please set data_type into "Bulk".

## Example for "SingleCell"
file_path <- system.file("extdata", "SingleCell_example.RData", package = "gdGSE")
load(file_path)
ls()
# [1] "file_path"  "SC_condition"  "SC_exp_matrix" "SC_Signature" 
gdGSE_Score <- gdGSE(SC_exp_matrix,SC_condition, SC_Signature,data_type = "SingleCell")  #If Users calculate enrichment score in SingleCell, please set data_type into "SingleCell".
```

&nbsp;

# Contact

E-mail any questions to Xiaosheng Wang (xiaosheng.wang@cpu.edu.cn)
