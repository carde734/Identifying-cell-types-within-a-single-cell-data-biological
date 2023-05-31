# Identifying-cell-types-within-a-single-cell-data-biological
Single cell data biological analysis includes a primary step that is investigating which type each cell belongs to CIPR (Correlation based method) and ScCATCH (Marker gene method) methods were analysed and numerically compared under different situations.


## Introduction
This project has the objective of supporting single cell data biological analysis by
investigating one important initial step of the process, that is, determine which type each
cell belongs to. To make such annotation, various machine learning methods were
proposed [Pasquini] and the most representative method from each of the two following
categories was chosen for this research. From the category “Marker gene method” the
chosen method was ScCATCH and from the category “Correlation based method” the
chosen method was CIPR. A detailed numerical comparison of these methods was made
in order to conclude which method obtained better results, under the following situations:\
\
• 1: The method contains exactly same cell types as cell types in the single cell data.\
• 2: The method contains more cell types than in single cell data.\
• 3: Some cell types that are in single cell data are not present in the method.\
• 4: There are various proportions of cell types in single cell data.

## Results 

### True Data 
![image](https://github.com/carde734/Identifying-cell-types-within-a-single-cell-data-biological/assets/90332007/cc2660a8-c895-4edd-b6eb-13c2644e1ed7)




# References
[Pasquini] Pasquini, G., Arias, J. E. R., Schäfer, P., & Busskamp, V. (2021). Automated
methods for cell type annotation on scRNA-seq data. Computational and Structural
Biotechnology Journal, 19, 961-969.\
[Ekiz] Ekiz, H.A., Conley, C.J., Stephens, W.Z. et al. CIPR: a web-based R/shiny app
and R package to annotate cell clusters in single cell RNA sequencing experiments. BMC
Bioinformatics 21, 191 (2020)\
[Shao] Shao, X., Liao, J., Lu, X., Xue, R., Ai, N., & Fan, X. (2020). scCATCH: Automatic
Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data. In
iScience (Vol. 23, Issue 3, p. 100882). Elsevier BV.
