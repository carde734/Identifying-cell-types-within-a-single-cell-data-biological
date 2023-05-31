# Research Project 

# Correlation based method - CIPR 

#################### LOAD LIBRARIES  #############################################

library(dplyr)
library(Seurat)
library(CIPR)
library(biomaRt)
library(Matrix)
library(caret)
library(ggplot2)

#################################################################################



########################### RUN HERE FOR STEP BY STEP ANALYSIS #################


########################### RUN DEPENDING ON WHICH SUBSET OF DATA ##############

# Loads full data set 
pbmc.data <- Matrix(main_data_vals, sparse = TRUE)


# Loads subset with different proportions of cells  
pbmc.data <- Matrix(main_data_vals_diff_proportions, sparse = TRUE)



# CIPR and ScCATCH common part END
#######################################################################

##### CIPR ONLY


# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data)


# PRE PROCESSING 

# NORMALIZING DATA
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)


# Variable gene detection and scaling
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Perform PCA - To reduce dimensionality and identify possible clusters
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:10, verbose = F)

# Plots the standard deviations of the principle components for easy 
#identification of an elbow in the graph. This elbow often corresponds well with
#the significant dims 
ElbowPlot(pbmc)

# The elbow plot identified 5 significant dimensions

# Cluster cells

pbmc <- FindNeighbors(pbmc, dims=1:5)


### ADJUST RESOLUTION IN ORDER TO GET DIFFERENT COMBINATIONS OF # OF CLUSTERS ###


# Run this command and the others below for 13 clusters - We used the default value 
# When using the data with diff proportions the output is 11 clusters
pbmc = FindClusters(pbmc, resolution = 0.8)


# Run this command and the others below for 5 clusters - Ground truth of number
# of clusters
pbmc = FindClusters(pbmc, resolution = 0.05)




###############################################################################

# Run non-linear dimensionality reduction (tSNE)
pbmc <- RunTSNE(pbmc, dims = 1:5)
pbmc$unnamed_clusters <- Idents(pbmc)


saveRDS(pbmc, "pbmc_4_diff_proportions_clusters.rds")

# Find differentially expressed genes

# This is the step where we generate the input for CIPR's log fold change
# (logFC) comparison methods.

allmarkers <- FindAllMarkers(pbmc)
saveRDS(allmarkers, "allmarkers_4_diff_proportions_clusters.rds")

#########################################################################


# pbmc and all markers saved for 13 and 5 clusters and 11 and 4 (different proportions)
# in order to compare! Load pbmc and all markers

pbmc <- readRDS("pbmc_5_clusters.rds")
allmarkers <- readRDS("allmarkers_5_clusters.rds")

pbmc<- readRDS("pbmc_13_dims_clusters.rds")
allmarkers<- readRDS("allmarkers_13_dims_clusters.rds")

pbmc <- readRDS("pbmc_4_diff_proportions_clusters.rds")
allmarkers <- readRDS("allmarkers_4_diff_proportions_clusters.rds")

pbmc <- readRDS("pbmc_11_diff_proportions_clusters.rds")
allmarkers <- readRDS("allmarkers_11_diff_proportions_clusters.rds")


############################### ANALYSIS ###############################

# Create labels with cells (B,T, etc)
cd_cells = pbmc@meta.data$orig.ident
cells = as.character(cd_cells)
cells=replace(cells,cells=="CD14","Monocyte") 
cells=replace(cells,cells=="CD19","B cell") 
cells=replace(cells,cells=="CD56","NK cell") 
cells=replace(cells,cells=="CD4","CD4+ T cell") 
cells=replace(cells,cells=="CD8","CD8+ T cell") 
cells = as.factor(cells)


# Insert column ori_cells into the metadata
pbmc@meta.data$'ori_cells' = cells


# Limiting analysis to the select subsets of reference data

# Method 1 -  The method contains exactly same cell types as cell types in 
# the single cell data

CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "hsrnaseq", 
     plot_ind = T,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T,
     select_ref_subsets = c("B cell","CD4+ T cell", "CD8+ T cell", "Monocyte", "NK cell"))






# Plot all identity scores per cluster-reference cell pairs (hsrnaseq)
# Method 2 and 4 (but method 4 uses a different subset) - The method contains
# more cell types than in single cell data

CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "hsrnaseq", 
     plot_ind = T,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T,
)





# Method 3 - Some cell types that are in single cell data are not present in 
#the method.
# Removing B cells, for example - How B cells will be classified?

# Load hsrnaseq_samples.rda from CIPR git hub package and get unique cells

all_reference_cell_type_hsrnaseq = c("B cell","Basophil","CD4+ T cell","CD8+ T cell","Dendritic cell","Monocyte"
,"Neutrophil","NK cell","Progenitor","MAIT-gdT")   

# Removing B cell from the list 
No_Bcell_reference_cell_type_hsrnaseq = all_reference_cell_type_hsrnaseq[-1]

CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "hsrnaseq", 
     plot_ind = T,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T,
     select_ref_subsets = No_Bcell_reference_cell_type_hsrnaseq)



# Plot top scoring refernce subsets for each cluster
  
CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "hsrnaseq", 
     plot_ind = F,
     plot_top = T, 
     global_results_obj = T, 
     global_plot_obj = T)




# Manually identify clusters through CIPR plots:
ind_clu_plots$cluster0
ind_clu_plots$cluster1
ind_clu_plots$cluster2
ind_clu_plots$cluster3
ind_clu_plots$cluster4
ind_clu_plots$cluster6
ind_clu_plots$cluster8
ind_clu_plots$cluster10
ind_clu_plots$cluster9



################# VISUALIZATION ###########################################

# Visualize Seurat pbject --- BEFORE RUNNING FINDALLMARKERS
DimPlot(pbmc) + ggtitle('My title')
visualize_clusters=DimPlot(pbmc, group.by = 'seurat_clusters', label = TRUE) +ggtitle("Seurat Clustered Data")
# Perfect for Report! - ORIGINAL DATA
plot_pca_original_data=DimPlot(pbmc, group.by = "orig.ident", label = TRUE)+ggtitle("True Data")
plot_pca_original_data 


plot_pca_original_data_cell=DimPlot(pbmc, group.by = "ori_cells", label = TRUE)+ggtitle("True Data")
plot_pca_original_data_cell 

visualize_clusters 
#visualize_clusters_splitori

# Perfect to compare! 
plot_pca_original_data_cell|visualize_clusters
###################### END VISUALIZATION #######################################




################### RUN HERE FOR QUICK ANALYSIS  ##############################


#################### FUNCTIONS FOR FULL ANALYSIS ###############################


####### FUNC 1 - CALCULATE ACCURACY #############################################

calc_accuracy_plot_results= function(pbmc, allmarkers)
{
  
  visualize_clusters=DimPlot(pbmc, group.by = 'seurat_clusters', label = TRUE) +ggtitle("Seurat Clustered Data")
  
  
  clusters = as.numeric(unique(pbmc@meta.data$seurat_clusters))-1
  results_top_genes = matrix(data=NA,nrow = length(clusters), ncol = 2)
  
  for (c in clusters)
  {
    res=head(subset(CIPR_top_results, cluster==c),1)
    cell=as.character(res$reference_cell_type)
    results_top_genes[c+1,]=c(c,cell) 
  }
  colnames(results_top_genes)=c("Cluster", "Classified Cell")
  
  
  
  
  # Identifiying cells into cluster number 
  cell_ident_cluster = pbmc$seurat_clusters
  
  cells_classification = as.character(cell_ident_cluster)
  
  
  for (c in clusters)
  {
    cell_name = results_top_genes[,2][c+1]
    cells_classification=replace(cells_classification,cells_classification==paste(c),cell_name) 
  }
  
  cells_classification = as.factor(cells_classification)
  pbmc@meta.data$'cells_classification' = cells_classification
  
  
  # Visual Countings per cluster 
  counts_cell_cluster = table(cell_ident_cluster)
  
  # Visual Countings per class of cell identified 
  counts_cell_class = table(cells_classification)
  
  plot_pca_original_data_cell = DimPlot(pbmc, group.by = "ori_cells", label = TRUE)+ggtitle("True Data")
  
  visualize_class_cell =DimPlot(pbmc, group.by = 'cells_classification', label = TRUE) +ggtitle("Classified Clusters Into Cells")
  # Compare true data with the classified classes trough plot 
  
  
  
  # Metadata - Cluster results 
  # Nice confusion matrix that can be built! just need to replace the cluster
  # numbers with the identified cells so we can get a vector from 1-12500 with
  # the obtained predictions!!
  Cluster_res = pbmc@meta.data[,c('orig.ident', "seurat_clusters")]
  # How to define the clusters with cells !
  
  # Interesting that the method actually successed in clustering some of the cells 
  # pretty well. Cluster 0 has only CD14 cells (Monotocytes) , etc ,etc 
  
  
  # Get true labels 
  true_cell=pbmc@meta.data$ori_cells
  predicted_cell=pbmc@meta.data$cells_classification

  
  tot_cells_correct = sum(as.character(predicted_cell)==as.character(true_cell))
  accuracy =tot_cells_correct /length(true_cell)
  
  
  conf_matrix = table(true_cell,predicted_cell)
  
  
  results = list("pbmc_object" = pbmc,
                 "all_markers_object" = allmarkers,
                 "seurat_clusters_ident" = cell_ident_cluster,
                 "objects per cluster"=  table(Cluster_res),
                 "cluster_classification" = Cluster_res,
                 "true_cell" = true_cell,
                 "predicted_cell"= predicted_cell,
                 "conf matrix"=conf_matrix,
                 "Number of Cells correctly identified" = tot_cells_correct,
                 "Number of Cells misclassifed" = (length(true_cell)-tot_cells_correct),
                 "Accuracy"=accuracy,
                 "Error"=(1-accuracy), 
                 "plot_seurat_clusters" = visualize_clusters,
                 "plot_TrueData"= plot_pca_original_data_cell,
                 "plot_predicted_cell" =visualize_class_cell,
                 "plot_DataVSPredicted" = plot_pca_original_data_cell|visualize_class_cell)
  return(results)
  
}



# Runs CIPR analysis and takes pbmc file, allmarkers file and a comparison method 
# Comparison methods can be: 


# Method_1: The method contains exactly same cell types as cell types in the 
# single cell data

# Method_2: The method contains more cell types than in single cell data

# Method_3: Some cell types that are in single cell data are not present in the 
# method.

# Method_4: There are various proportions of cell types in single cell data
# (also method has all cells)

# Top_Score: More like a general overview over the genes that are the most higly
# expressed among each cluster



####### FUNC 2 - RUNS ANALYSIS PER PBMC AND METHOD ##############################


runs_CIPR_analysis = function( pbmc_file, allmarkers_file,comparison_method)
{
  
  pbmc = readRDS(paste0(pbmc_file))
  allmarkers = readRDS(paste0(allmarkers_file))

  # Create labels with cells (B,T, etc)
  cd_cells = pbmc@meta.data$orig.ident
  cells = as.character(cd_cells)
  cells=replace(cells,cells=="CD14","Monocyte") 
  cells=replace(cells,cells=="CD19","B cell") 
  cells=replace(cells,cells=="CD56","NK cell") 
  cells=replace(cells,cells=="CD4","CD4+ T cell") 
  cells=replace(cells,cells=="CD8","CD8+ T cell") 
  cells = as.factor(cells)
  
  
  # Insert column ori_cells into the metadata
  pbmc@meta.data$'ori_cells' = cells
  
  
  if (comparison_method == 'Method_1'){
    
    # Limiting analysis to the select subsets of reference data
    # Method 1 -  The method contains exactly same cell types as cell types in 
    # the single cell data
    
    
     suppressWarnings(
       
       CIPR(input_dat = allmarkers,
            comp_method = "logfc_dot_product", 
            reference = "hsrnaseq", 
            plot_ind = T,
            plot_top = F, 
            global_results_obj = T, 
            global_plot_obj = T,
            select_ref_subsets = c("B cell","CD4+ T cell", "CD8+ T cell", "Monocyte", "NK cell"))
       
       
       
     )
    

    
    results_list = 
      list("Cluster ind plots"=ind_clu_plots,
   "Accuracy and general plots" =calc_accuracy_plot_results(pbmc = pbmc, 
                                                      allmarkers = allmarkers))
    
    return(results_list)
    
    
  }
  
  
  if (comparison_method == 'Method_2' || comparison_method == 'Method_4' ){
    
    # Plot all identity scores per cluster-reference cell pairs (hsrnaseq)
    # Method 2 - The method contains more cell types than in single cell data
    CIPR(input_dat = allmarkers,
         comp_method = "logfc_dot_product", 
         reference = "hsrnaseq", 
         plot_ind = T,
         plot_top = F, 
         global_results_obj = T, 
         global_plot_obj = T,
    )
    
    results_list = 
      list("Cluster ind plots"=ind_clu_plots,
           "Accuracy and general plots" =calc_accuracy_plot_results(pbmc = pbmc, 
                                                                    allmarkers = allmarkers))
    
    return(results_list)    
  }
  
  if (comparison_method == 'Method_3'){
    
    # Method 3 - Some cell types that are in single cell data are not present in 
    #the method.
    # Removing B cells, for example - How B cells will be classified?
    
    CIPR(input_dat = allmarkers,
         comp_method = "logfc_dot_product", 
         reference = "hsrnaseq", 
         plot_ind = T,
         plot_top = F, 
         global_results_obj = T, 
         global_plot_obj = T,
         select_ref_subsets = c("Basophil","CD4+ T cell","CD8+ T cell","Dendritic cell","Monocyte"
                                ,"Neutrophil","NK cell","Progenitor","MAIT-gdT") )
    
    results_list = 
      list("Cluster ind plots"=ind_clu_plots,
           "Accuracy and general plots" =calc_accuracy_plot_results(pbmc = pbmc, 
                                                                    allmarkers = allmarkers))
    
    return(results_list)      
    
  }
  
  
  if (comparison_method == 'Top_Scoring'){
    
    # Plot top scoring reference subsets for each cluster
    
    CIPR(input_dat = allmarkers,
         comp_method = "logfc_dot_product", 
         reference = "hsrnaseq", 
         plot_ind = F,
         plot_top = T, 
         global_results_obj = T, 
         global_plot_obj = T)
    
    
    
  }

  
}


######################### END FUNCTIONS CREATION ##############################


####################### ANALYSIS RESULTS #######################################

start.time <- Sys.time()


# Calling the functions in order to get the analysis objects for each situation


# Larger number of clusters

analysis_13c_method_1= runs_CIPR_analysis(pbmc_file ="pbmc_13_dims_clusters.rds",
                                 allmarkers_file ="allmarkers_13_dims_clusters.rds",
                                 comparison_method = "Method_1")


analysis_13c_method_2= runs_CIPR_analysis(pbmc_file ="pbmc_13_dims_clusters.rds",
                                                allmarkers_file ="allmarkers_13_dims_clusters.rds",
                                                comparison_method = "Method_2")



analysis_13c_method_3= runs_CIPR_analysis(pbmc_file ="pbmc_13_dims_clusters.rds",
                                          allmarkers_file ="allmarkers_13_dims_clusters.rds",
                                          comparison_method = "Method_3")





analysis_13c_top_scoring= runs_CIPR_analysis(pbmc_file ="pbmc_13_dims_clusters.rds",
                                      allmarkers_file ="allmarkers_13_dims_clusters.rds",
                                      comparison_method = "Top_Scoring")




# Method 4 - There are various proportions of cell types in single 
# cell data + labeling method with all available cell types 
analysis_11c_Method_4= runs_CIPR_analysis(pbmc_file ="pbmc_11_diff_proportions_clusters.rds",
                                         allmarkers_file ="allmarkers_11_diff_proportions_clusters.rds",
                                         comparison_method = "Method_4")





# Smaller number of clusters



analysis_5c_method_1= runs_CIPR_analysis(pbmc_file ="pbmc_5_clusters.rds",
                                          allmarkers_file ="allmarkers_5_clusters.rds",
                                          comparison_method = "Method_1")


analysis_5c_method_2= runs_CIPR_analysis(pbmc_file ="pbmc_5_clusters.rds",
                                          allmarkers_file ="allmarkers_5_clusters.rds",
                                          comparison_method = "Method_2")



analysis_5c_method_3= runs_CIPR_analysis(pbmc_file ="pbmc_5_clusters.rds",
                                          allmarkers_file ="allmarkers_5_clusters.rds",
                                          comparison_method = "Method_3")




# Method 4 - There are various proportions of cell types in single 
# cell data + labeling method with all available cell types 
analysis_4c_Method_4= runs_CIPR_analysis(pbmc_file ="pbmc_4_diff_proportions_clusters.rds",
                                          allmarkers_file ="allmarkers_4_diff_proportions_clusters.rds",
                                          comparison_method = "Method_4")





end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

# Time difference of 4.73 mins


################## CALS OBJECTS (EXAMPLES) ###################################

analysis_13c_method_1$`Accuracy and general plots`

# Sort of Confusion Matrix 
analysis_13c_method_1$`Accuracy and general plots`$`conf matrix`

# Accuracy 
analysis_13c_method_1$`Accuracy and general plots`$Accuracy

# Other info 
analysis_13c_method_1$`Accuracy and general plots`$`Number of Cells correctly identified`
analysis_13c_method_1$`Accuracy and general plots`$`Number of Cells misclassifed`

# Plots  
analysis_13c_method_1$`Accuracy and general plots`$plot_seurat_clusters
analysis_13c_method_1$`Accuracy and general plots`$plot_TrueData
analysis_13c_method_1$`Accuracy and general plots`$plot_predicted_cell
analysis_13c_method_1$`Accuracy and general plots`$plot_DataVSPredicted


#################### SOME EXTRA ANALYSIS ####################################

#################### Manually identify clusters through CIPR plots ############

analysis_13c_top_scoring

analysis_13c_method_2$`Cluster ind plots`$cluster0
analysis_13c_method_2$`Cluster ind plots`$cluster6
  
################ TOP FEATURE ANALYSIS - 13 CLUSTERS - DEFAULT METHOD ###########

analysis_13c_method_2$`Cluster ind plots`$cluster0
analysis_13c_method_2$`Cluster ind plots`$cluster6


# Testing feature plotting to analyse cluster 6, since it's a cluster 
# with verified presence (very high reference idenitity score compared to the 
# other cell markers) of B cell markers.

# Cluster 6 - Marked as B cell 

head(subset(analysis_13c_method_2$`Accuracy and general plots`$all_markers_object, cluster==6))

# Observe first line for gene CD79A (Searched for the gene in PanglaoDB and this
# gene is clearly more present in B cells than other cells)

# This gene is detected in 96% percent of the cells in cluster 6 compared to 13% 
# of the cells in all the other clusters combined.

# Let's visualize this top feature
topfeature_CD79A=FeaturePlot(analysis_13c_method_2$`Accuracy and general plots`$pbmc_object, features = c('CD79A'), min.cutoff = 'q5') +
  ggtitle('Visualize Top Feature for Gene CD79A ')

topfeature_CD79A

analysis_13c_plot_TrueData = analysis_13c_method_2$`Accuracy and general plots`$plot_TrueData

analysis_13c_plot_TrueData|topfeature_CD79A


