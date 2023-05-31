# Research Project 
# Load Data

############ INSTALL PACKAGES ##################################################


install.packages('Seurat')

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("atakanekiz/CIPR-Package", build_vignettes = F)



# Install package scCATCH
install.packages("scCATCH")


# Install package rlang
install.packages("rlang")

#Install package to convert ensembled genes to gene symbols
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt", force = TRUE)

# Convert ensembled ID using annotables 
install.packages("devtools")
devtools::install_github("stephenturner/annotables")


#################### LOAD LIBRARIES  #############################################
library(Matrix)
#################################################################################

gc()
memory.limit(9999999999)

############################ PRE PROCESSING DATA ##############################

# CIPR and ScCATCH common part START
################################################################################
# File was pre-processed in python and ensembled genes were converted to 
# symbols 

main_data_vals = read.delim(file="python_pre_processed_5celltypes_2500each.txt",
                            header = TRUE, sep = "\t")


rows_main_data = main_data_vals[,1]
rownames(main_data_vals) = rows_main_data
main_data_vals = as.matrix(main_data_vals[,-1])
cols_main_data = colnames(main_data_vals)



# Loads subset of data
# Situation There are various proportions of cell types in single cell data
# 2000 t (1000 CD4 and 1000 CD8), cells 500 nk cells, 100 b cells, 800 Monocytes

# cols_main_data[1:2500] -> CD4 (T)
# cols_main_data[2501:5000] -> CD14 (Monocytes)
# cols_main_data[5001:7500] -> CD8 (T)
# cols_main_data[7501:10000] -> CD19 (B)
# cols_main_data[10001:12500] -> CD56 (NK)



set.seed(123)
CD4_T = sample(seq(from = 1, to = 2500, by = 1), size = 1000)
CD8_T = sample(seq(from = 5001, to = 7500, by = 1), size = 1000)
CD14_Monocytes = sample(seq(from = 2501, to = 5000, by = 1), size = 800)
CD19_B = sample(seq(from = 7501, to = 10000, by = 1), size = 100)
CD56_NK = sample(seq(from = 10001, to = 12500, by = 1), size = 500)


main_data_vals_diff_proportions = main_data_vals[,c(CD4_T,CD8_T,CD14_Monocytes,
                                                    CD19_B,CD56_NK)]




########################### RUN DEPENDING ON WHICH SUBSET OF DATA #############

# Loads full data set 
#pbmc.data <- Matrix(main_data_vals, sparse = TRUE)


# Loads subset with different proportions of cells  
#pbmc.data <- Matrix(main_data_vals_diff_proportions, sparse = TRUE)


