# Research Project 

# Marker gene method Analysis - scCATCH 

######################## LOAD LIBRARIES  #########################################
library(Seurat)
library(rlang)
library(scCATCH)
library(ggplot2)



########################### RUN DEPENDING ON WHICH SUBSET OF DATA #############

# Loads full data set 
pbmc_full_data <- Matrix(main_data_vals, sparse = TRUE)


# Loads subset with different proportions of cells  
pbmc_diff_prop_data <- Matrix(main_data_vals_diff_proportions, sparse = TRUE)





##################### REFERENCE DATABASES ####################################
load(file = "cellmatch.rda")
# Transform database in cell match in order to only return the celltype (no subtype)
# , with the exception of CD8+ T Cell and CD4+ T Cell 
cellmatch_withoutSubtype = cellmatch
cellmatch_withoutSubtype$subtype1='NA'
cellmatch_withoutSubtype$subtype2[!(cellmatch_withoutSubtype$subtype2 %in% c('CD8+','CD4+') )] = 'NA'
cellmatch_withoutSubtype$subtype3='NA'






# Method 1 - Subset Cells to The method contains exactly same cell types as cell types in the single cell data

cellmatch_exact_celltype = cellmatch_withoutSubtype
cellmatch_exact_celltype= subset(cellmatch_exact_celltype,
                                 cellmatch_exact_celltype$celltype %in% c('Monocyte',
                                                                          'T Cell',
                                                                          'B Cell',
                                                                          'Natural Killer T (NKT) Cell',
                                                                          'Killer Cell'))



# Method 2 and 4 - All available Cells
cellmatch_allcell = cellmatch_withoutSubtype

# Method 3
# Some cell types that are in single cell data are not present in the method. 
#Subset for all cells, besides B cells. How will B cells be classified? 

cellmatch_noBcells = subset(cellmatch_allcell,cellmatch_allcell$celltype != 'B Cell')







################### RUN FOR STEP BY STEP ANALYSIS  ##############################
# ScCATCH ONLY 

#################################################################################

# Use as clusters original clusters used in CIPR: 
# 13 clusters (main case)
# 10 clusters - Different proportions of cells 


# CLUSTER OBJECTS CREATED from pbmc objects created and saved from 
# research_2_proj_CorrelationCIPR_pythonPP.R

pbmc_13c_object =readRDS("pbmc_13_dims_clusters.rds")
pbmc_10c_object =readRDS("pbmc_10_diff_proportions_clusters.rds")


# 13 clusters - Make sure that this pbmc object corresponds to the 13 clusters
clusters_13clusters = pbmc_13c_object$seurat_clusters


clusters_10clusters = pbmc_10c_object$seurat_clusters

# All pbmc objects have the original identification = Ground truth clusters
clusters_ground_truth = pbmc_13c_object$orig.ident


true_cell = clusters_ground_truth
true_cell = as.character(true_cell)
true_cell=replace(true_cell,true_cell=="CD14","Monocyte") 
true_cell=replace(true_cell,true_cell=="CD19","B Cell") 
true_cell=replace(true_cell,true_cell=="CD56","NK Cell") 
true_cell=replace(true_cell,true_cell=="CD4","CD4+ T Cell") 
true_cell=replace(true_cell,true_cell=="CD8","CD8+ T Cell") 
true_cell = as.factor(true_cell)




###############################################################################
data_sccatch = pbmc_full_data

# CHANGE HERE FOR DIFFERENT CLUSTER FORMATIONS:
#cluster_data_sccatch = as.character(clusters_13clusters)
#cluster_data_sccatch = as.character(clusters_5clusters)
#cluster_data_sccatch = as.character(clusters_10clusters)
cluster_data_sccatch = as.character(clusters_ground_truth)




# NORMALIZING DATA
data_sccatch <- NormalizeData(data_sccatch, normalization.method = "LogNormalize", scale.factor = 10000)
data_sccatch <- NormalizeData(data_sccatch)

all_genes = rownames(data_sccatch)

data_sccatch <- ScaleData(data_sccatch, features = all_genes)

################################## ANALYSIS ##############################################################
# revise gene symbols
data_sccatch = rev_gene(data = data_sccatch , data_type = "data",
                        species = "Human", geneinfo = geneinfo)



my_obj <- createscCATCH(data = data_sccatch , cluster = cluster_data_sccatch)




#################################################################################
my_obj <- findmarkergene(object = my_obj,
                         species = "Human",
                         marker = cellmatch,
                         tissue = "Peripheral blood")


my_obj <- findcelltype(object = my_obj)

my_obj@celltype

 
##################### Analysis ######################################

# Put this part inside a function to print plots and confusion matrix
# Do a manual calculation of accuracy through confusion matrix?? 


  
cluster_numbers = my_obj@celltype$cluster
cell_classified =my_obj@celltype$cell_type

cluster_classified = cluster_data_sccatch

for (i in 1:length(cluster_numbers)){
  cluster_classified=replace(cluster_classified,
                             cluster_classified==cluster_numbers[i],cell_classified[i]) 
}


predicted_cell = cluster_classified

conf_matrix = table(true_cell,predicted_cell)
conf_matrix

pbmc_13c_object@meta.data$'ori_cells_sccatch' = true_cell
pbmc_13c_object@meta.data$'cluster_classification_sccatch' = cluster_classified

plot_pca_original_data_cell = DimPlot(pbmc_13c_object, group.by = "ori_cells_sccatch", label = TRUE)+ggtitle("True Data")

visualize_class_cell =DimPlot(pbmc_13c_object, group.by = 'cluster_classification_sccatch', label = TRUE) +ggtitle("Classified Clusters Into Cells")

# Compare true data with the classified classes trough plot 
plot_pca_original_data_cell|visualize_class_cell



########################### RUN HERE FOR QUICK ANALYSIS #######################


############################ FUNCTION FOR ANALYSIS ############################

# Methods : 

# Method_1: The method contains exactly same cell types as cell types in the 
# single cell data

# Method_2: The method contains more cell types than in single cell data

# Method_3: Some cell types that are in single cell data are not present in the 
# method.

# Method_4: There are various proportions of cell types in single cell data
# (also method has all cells)


runs_ScCATCH_analysis = function(pbmc_file,comparison_method){
  
  
  if (comparison_method =='Method_1')
  {
    data_sccatch = pbmc_full_data
    reference_database = cellmatch_exact_celltype
    
  }
  
  
  if (comparison_method =='Method_2')
  {
    data_sccatch = pbmc_full_data
    reference_database = cellmatch_allcell
  }
  
  
  
  if (comparison_method =='Method_3')
  {
    data_sccatch = pbmc_full_data
    reference_database = cellmatch_noBcells
    
  }
  
  
  
  if (comparison_method =='Method_4')
  {
    data_sccatch = pbmc_diff_prop_data
    reference_database = cellmatch_allcell
    
  }
  
  
  pbmc_object = readRDS(paste0(pbmc_file))    
  cluster_data_sccatch = as.character(pbmc_object$seurat_clusters)
  
  
  # All pbmc objects have the original identification = Ground truth clusters
  clusters_ground_truth = pbmc_object$orig.ident
  
  
  true_cell = clusters_ground_truth
  true_cell = as.character(true_cell)
  true_cell=replace(true_cell,true_cell=="CD14","Monocyte") 
  true_cell=replace(true_cell,true_cell=="CD19","B Cell") 
  true_cell=replace(true_cell,true_cell=="CD56","NK Cell") 
  true_cell=replace(true_cell,true_cell=="CD4","CD4+ T Cell") 
  true_cell=replace(true_cell,true_cell=="CD8","CD8+ T Cell") 
  true_cell = as.factor(true_cell)
  
  
  
  # NORMALIZING DATA
  data_sccatch = NormalizeData(data_sccatch, normalization.method = "LogNormalize", scale.factor = 10000)
  data_sccatch = NormalizeData(data_sccatch)
  
  all_genes = rownames(data_sccatch)
  
  data_sccatch = ScaleData(data_sccatch, features = all_genes)
  
  # revise gene symbols
  data_sccatch = rev_gene(data = data_sccatch , data_type = "data",
                          species = "Human", geneinfo = geneinfo)
  
  
  
  my_obj = createscCATCH(data = data_sccatch , cluster = cluster_data_sccatch)
  
  
  my_obj = findmarkergene(object = my_obj,
                           species = "Human",
                           marker = reference_database,
                           tissue = "Peripheral blood")
  
  
  my_obj = findcelltype(object = my_obj)
  
  result_celltype = my_obj@celltype
  
  
  cluster_numbers = my_obj@celltype$cluster
  cell_classified =my_obj@celltype$cell_type
  
  cluster_classified = cluster_data_sccatch
  
  for (i in 1:length(cluster_numbers)){
    cluster_classified=replace(cluster_classified,
                               cluster_classified==cluster_numbers[i],cell_classified[i]) 
  }
  
  
  predicted_cell = cluster_classified
  
  conf_matrix = table(true_cell,predicted_cell)
  
  true_cell_groupTcell = as.character(true_cell)
  true_cell_groupTcell[agrep('+ T', true_cell_groupTcell)] ='T Cell'
  true_cell_groupTcell = as.factor(true_cell_groupTcell)
  
  
  predicted_cell_groupTcell = as.character(predicted_cell)
  predicted_cell_groupTcell[agrep('+ T', predicted_cell_groupTcell)] ='T Cell'
  predicted_cell_groupTcell = as.factor(predicted_cell_groupTcell)
  
  
  tot_cells_correct = sum(as.character(predicted_cell_groupTcell)==as.character(true_cell_groupTcell))
  accuracy =tot_cells_correct /length(true_cell_groupTcell)
  
  
  
  pbmc_object@meta.data$'ori_cells_sccatch' = true_cell
  pbmc_object@meta.data$'cluster_classification_sccatch' = cluster_classified
  
  plot_pca_original_data_cell = DimPlot(pbmc_object, group.by = "ori_cells_sccatch", label = TRUE)+ggtitle("True Data")
  
  visualize_class_cell =DimPlot(pbmc_object, group.by = 'cluster_classification_sccatch', label = TRUE) +ggtitle("Classified Clusters Into Cells")
  
  # Compare true data with the classified classes trough plot 
  
  
  results = list("Results ScCATCH" = result_celltype,
                 "Confusion Matrix" = conf_matrix,
                 "Number of Cells correctly identified" = tot_cells_correct,
                 "Number of Cells misclassifed" = (length(true_cell)-tot_cells_correct),
                 "Accuracy"=accuracy,
                 "Error"=(1-accuracy),
                 "plot_TrueData"= plot_pca_original_data_cell,
                 "plot_predicted_cell" =visualize_class_cell,
                 "plot_DataVSPredicted" =  plot_pca_original_data_cell|visualize_class_cell
                 )
  
  return(results)

}



####################### ANALYSIS RESULTS #######################################


# Larger number of clusters


analysis_ScCATCH_13c_method_1= runs_ScCATCH_analysis(pbmc_file ="pbmc_13_dims_clusters.rds",
                                          comparison_method = "Method_1")


analysis_ScCATCH_13c_method_2= runs_ScCATCH_analysis(pbmc_file ="pbmc_13_dims_clusters.rds",
                                          comparison_method = "Method_2")


analysis_ScCATCH_13c_method_3= runs_ScCATCH_analysis(pbmc_file ="pbmc_13_dims_clusters.rds",
                                                     comparison_method = "Method_3")


analysis_ScCATCH_11c_method_4= runs_ScCATCH_analysis(pbmc_file ="pbmc_11_diff_proportions_clusters.rds",
                                                     comparison_method = "Method_4")



# Smaller number of clusters


analysis_ScCATCH_5c_method_1= runs_ScCATCH_analysis(pbmc_file ="pbmc_5_clusters.rds",
                                                     comparison_method = "Method_1")


analysis_ScCATCH_5c_method_2= runs_ScCATCH_analysis(pbmc_file ="pbmc_5_clusters.rds",
                                                     comparison_method = "Method_2")


analysis_ScCATCH_5c_method_3= runs_ScCATCH_analysis(pbmc_file ="pbmc_5_clusters.rds",
                                                     comparison_method = "Method_3")




analysis_ScCATCH_4c_method_4= runs_ScCATCH_analysis(pbmc_file ="pbmc_4_diff_proportions_clusters.rds",
                                                    comparison_method = "Method_4")




################## CALS OBJECTS (EXAMPLES) ###################################

# Results of ScCATCH 
analysis_ScCATCH_5c_method_1$`Results ScCATCH`

# Sort of Confusion Matrix 
analysis_ScCATCH_5c_method_1$`Confusion Matrix`

# Accuracy 
analysis_ScCATCH_5c_method_1$Accuracy

# Other info 
analysis_ScCATCH_5c_method_1$`Number of Cells correctly identified`
analysis_ScCATCH_5c_method_1$`Number of Cells misclassifed`

# Plots  
analysis_ScCATCH_5c_method_1$plot_TrueData
analysis_ScCATCH_5c_method_1$plot_predicted_cell
analysis_ScCATCH_5c_method_1$plot_DataVSPredicted

