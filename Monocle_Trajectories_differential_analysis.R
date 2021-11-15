


#link - https://cole-trapnell-lab.github.io/monocle3/docs/differential/

#Let's begin with a small set of genes that we know are important in ciliated neurons to demonstrate Monocle's capabilities:
ciliated_genes <- c("che-1","hlh-17","nhr-6", "dmd-6", "ceh-36","ham-1")




cds_subset <- cds_traj_order[rowData(cds_traj_order)$gene_short_name %in% ciliated_genes,]
#cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]


#identify genes that vary over time by fitting this model to each one, and then testing whether its 
#?? is significantly different from zero. To do so, we first call the fit_models() function:

gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time")

#Now let's see which of these genes have time-dependent expression. First, we extract a table of coefficients from each model using the coefficient_table() function:
fit_coefs <- coefficient_table(gene_fits)
fit_coefs

#Note that the table includes one row for each term of each gene's model. We generally don't care about the intercept term 
#??0, so we can easily just extract the time terms:
emb_time_terms <- fit_coefs %>% filter(term == "embryo.time")


#Now, let's pull out the genes that have a significant time component.
#coefficient_table() tests whether each coefficient differs significantly from zero under the Wald test
#adjusted values can be found in the q_value column. We can filter the results and control the false discovery rate as follows:


emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

#visualize the differences revealed by the tests above. One type of plot is a "violin" plot.
plot_genes_violin(cds_subset, group_cells_by="embryo.time.bin", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


#Controlling for batch effects and other factors
gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time + batch")
fit_coefs <- coefficient_table(gene_fits)
fit_coefs %>% filter(term != "(Intercept)") %>%
  select(gene_short_name, term, q_value, estimate)


#Evaluating models of gene expression. We can evaluate the fits of each model using the evaluate_fits() function:

evaluate_fits(gene_fits)


#Should we include the batch term in our model of gene expression or not? 
#Monocle provides a function compare_models() that can help you decide. 
#Compare models takes two models and returns the result of a likelihood ratio test between them. 
#Any time you add terms to a model, it will improve the fit. 
#But we should always to use the simplest model we can to explain our data. 
#The likelihood ratio test helps us decide whether the improvement in fit is large enough to justify the complexity our extra terms introduce. You run compare_models() like this:


#The first of the two models is called the full model. 
#This model is essentially a way of predicting the expression value of each gene in a given cell knowing both what time it was collected and which batch of cells it came from. 
time_batch_models <- fit_models(cds_subset,
                                model_formula_str = "~embryo.time + batch",
                                expression_family="negbinomial")
time_batch_models

#second model, called the reduced model, does the same thing, 
#but it only knows about the time each cell was collected.
time_models <- fit_models(cds_subset,
                          model_formula_str = "~embryo.time",
                          expression_family="negbinomial")
compare_models(time_batch_models, time_models) %>% select(gene_short_name, q_value)

#As we can see, all of the genes' likelihood ratio tests are significant, indicating that there are substantial batch effects in the data. 
#We are therefore justified in adding the batch term to our model.

#Choosing a distribution for modeling gene expression


#######################################################
#Finding genes that change as a function of pseudotime
#Identifying the genes that change as cells progress along a trajectory is a core objective of this type of analysis. 
#Knowing the order in which genes go on and off can inform new models of development



plot_cells(cds_traj_order, color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE)
#How do we find the genes that are differentially expressed on the different paths through the trajectory? How do we find the ones that are restricted to the beginning of the trajectory? Or excluded from it?

#Let's perform graph_test(), this time passing it neighbor_graph="principal_graph", which tells it to test whether cells at similar positions on the trajectory have correlated expression:
ciliated_cds_pr_test_res<- graph_test(cds_traj_order, neighbor_graph="principal_graph", cores=8)

pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
pr_deg_ids 
#Here are a couple of interesting genes that score as highly significant according to graph_test():
plot_cells(cds_traj_order, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

#we can collect the trajectory-variable genes into modules:

#gene_module_df <- find_gene_modules(cds_traj_order[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
gene_module_df <- find_gene_modules(cds_traj_order[pr_deg_ids,])

#Here we plot the aggregate module scores within each group of cell types as annotated by Packer & Zhu et al:

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_traj_order)), 
                                cell_group=colData(cds_traj_order)$cell.type)

agg_mat <- aggregate_gene_expression(cds_traj_order, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,scale="column", clustering_method="ward.D2")


#We can also pass gene_module_df to plot_cells() as we did when we compared clusters in the L2 data above.

plot_cells(cds_traj_order,
           genes=gene_module_df %>% filter(module %in% c(1, 2,5, 6, 7, 8)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


#Monocle offers another plotting function that can sometimes give a clearer view of a gene's dynamics along a single path. 
#You can select a path with choose_cells() or by subsetting the cell data set by cluster, cell type, or other annotation that's restricted to the path. Let's pick one such path, the AFD cells:

AFD_genes <- c("gcy-8", "dac-1", "oig-8")
AFD_lineage_cds <- cds_traj_order[rowData(cds_traj_order)$gene_short_name %in% AFD_genes,
                       colData(cds_traj_order)$cell.type %in% c("AFD")]


#The function plot_genes_in_pseudotime() takes a small set of genes 
#and shows you their dynamics as a function of pseudotime:

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="embryo.time.bin",
                         min_expr=0.5)


#Analyzing branches in single-cell trajectories
#Analyzing the genes that are regulated around trajectory branch points often provides insights into the genetic circuits that control cell fate decisions

cds_subset <- choose_cells(cds_traj_order)

#This will identify genes with interesting patterns of expression that fall only within the region of the trajectory you selected, giving you a more refined and relevant set of genes.

subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)

pr_deg_ids1 <- row.names(subset(subset_pr_test_res, q_value < 0.05))

#Grouping these genes into modules can reveal fate specific genes or those that are activate immediate prior to or following the branch point:

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids1,], resolution=0.001)

#We will organize the modules by their similarity (using hclust) over the trajectory so it's a little easier to see which ones come on before others:

agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
















#Graph-autocorrelation analysis for comparing clusters
#In the L2 worm data, we identified a number of clusters that were very distinct as neurons:

# reload and reprocess the data as described in the 'Clustering and classifying your cells' section
expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds_pro <- preprocess_cds(cds, num_dim = 50)
cds_dim <- reduce_dimension(cds_pro)
cds_clus <- cluster_cells(cds_dim, resolution=1e-5)

plot_cells(cds_clus, color_cells_by="partition")

colData(cds_clus)$assigned_cell_type <- as.character(partitions(cds_clus))

library(dplyr)
colData(cds_clus)$assigned_cell_type <- dplyr::recode(colData(cds_clus)$assigned_cell_type,
"1"="Germline", "2"="Body wall muscle","3"="Unclassified neurons",
"4"="Vulval precursors","5"="Failed QC","6"="Seam cells",
"7"="Pharyngeal epithelia", "8"="Coelomocytes","9"="Am/PH sheath cells",
"10"="Failed QC", "11"="Touch receptor neurons","12"="Intestinal/rectal muscle",
"13"="Pharyngeal neurons","14"="NA","15"="flp-1(+) interneurons",
"16"="Canal associated neurons","17"="Ciliated sensory neurons",
"18"="Other interneurons", "19"="Pharyngeal gland","20"="Failed QC",
 "21"="Ciliated sensory neurons","22"="Oxygen sensory neurons",
"23"="Ciliated sensory neurons", "24"="Ciliated sensory neurons",
"25"="Ciliated sensory neurons","26"="Ciliated sensory neurons",
"27"="Oxygen sensory neurons","28"="Ciliated sensory neurons",
"29"="Unclassified neurons","30"="Socket cells","31"="Failed QC",
"32"="Pharyngeal gland", "33"="Ciliated sensory neurons",
"34"="Ciliated sensory neurons","35"="Ciliated sensory neurons",
 "36"="Failed QC","37"="Ciliated sensory neurons","38"="Pharyngeal muscle")



plot_cells(cds_clus, group_cells_by="partition", color_cells_by="assigned_cell_type")

#Subset just the neurons:

neurons_cds <- cds_clus[,grep("neurons", colData(cds_clus)$assigned_cell_type, ignore.case=TRUE)]
plot_cells(neurons_cds, color_cells_by="partition")
plot_cells(neurons_cds, group_cells_by="partition", color_cells_by="assigned_cell_type")



#To investigate which genes are expressed differentially across the clusters, 
#we could use the regression analysis tools discussed above. 
#However, Monocle provides an alternative way of finding genes that vary between
#groups of cells in UMAP or t-SNE space. The function graph_test() uses a statistic from spatial autocorrelation analysis called Moran's I
pr_graph_test_res <- graph_test(neurons_cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))


# If you'd like to rank the genes by effect size, sort this table by the morans_Icolumn, which ranges from -1 to +1. A value of 0 indicates no effect, while +1 indicates perfect positive autocorrelation and suggests that nearby cells have very similar values of a gene's expression. Significant values much less than zero are generally rare.
#Positive values indicate a gene is expressed in a focal region of the UMAP space (e.g. specific to one or more clusters).

#Finding modules of co-regulated genes
#Once you have a set of genes that vary in some interesting way across the clusters, 
#Monocle provides a means of grouping them into modules.
#You can call find_gene_modules(), which essentially runs UMAP on the genes (as opposed to the cells) and then groups them into modules using Louvain community analysis:
  
gene_module_df <- find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)

#The data frame gene_module_df contains a row for each gene and identifies the module it belongs to
#Monocle provides a simple utility function called aggregate_gene_expression for this purpose:

cell_group_df <- tibble::tibble(cell=row.names(colData(neurons_cds)), 
                                cell_group=partitions(cds_clus)[colnames(neurons_cds)])
agg_mat <- aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
          scale="column", clustering_method="ward.D2",fontsize=6)

#The second way of looking at modules and their expression is to pass gene_module_df directly to plot_cells(). If there are many modules, it can be hard to see where each one is expressed, so we'll just look at a subset of them:

plot_cells(neurons_cds, 
           genes=gene_module_df %>% filter(module %in% c(8, 28, 33, 37)),
           group_cells_by="partition", color_cells_by="partition",
           show_trajectory_graph=FALSE)




