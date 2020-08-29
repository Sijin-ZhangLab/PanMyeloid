### NicheNet analysis for cDC1-derived and cDC2-derived LAMP3+ DC

## /data2/csj/tools/R-3.6.3/bin/R
.libPaths("/data2/csj/tools/Rlib3.6.3")

library(nichenetr)
library(tidyverse)

ligand_target_matrix = readRDS("ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns


h5ad <- readRDS("/data2/csj/Pan_Myeloid/A20191105/NicheNet_analysis/cDC_for_NicheNet.rds")

DC1_sig <- readRDS("/data2/csj/Pan_Myeloid/A20191105/NicheNet_analysis/cDC1_maturation_sig.rds")
DC2_sig <- readRDS("/data2/csj/Pan_Myeloid/A20191105/NicheNet_analysis/cDC2_maturation_sig.rds")

DC1_tmp <- DC1_sig[which(DC1_sig$Adj_pval < 0.05 & DC1_sig$log2FC > 0.2 ),]
DC2_tmp <- DC2_sig[which(DC2_sig$Adj_pval < 0.05 & DC2_sig$log2FC > 0.15),]
overlapped_genes <- intersect(DC1_tmp$gene,DC2_tmp$gene)
DC1_specific_vector <- DC1_tmp[!(DC1_tmp$gene %in% overlapped_genes),]
DC2_specific_vector <- DC2_tmp[!(DC2_tmp$gene %in% overlapped_genes),]
library(dplyr)
DC1_specific_df <- DC1_specific_vector %>%  top_n(n = 100, wt = log2FC)
DC2_specific_df <- DC2_specific_vector %>%  top_n(n = 100, wt = log2FC)
common_sig_df <- DC1_tmp[DC1_tmp$gene %in% overlapped_genes, ] %>%  top_n(n = 100, wt = log2FC)

DC1_specific <- as.vector(DC1_specific_df$gene)
DC2_specific <- as.vector(DC2_specific_df$gene)
common_sig <- as.vector(common_sig_df$gene)

df_genes <- data.frame("cDC1-specific" = DC1_specific, "cDC2-specific" = DC2_specific, "shared" = common_sig)
write.csv(df_genes,"cDC maturation sigantures.csv", quote = F, row.names = F, col.names = F)

DC2_specific <- DC2_specific[DC2_specific!='IGLV2-14']
DC2_specific <- DC2_specific[DC2_specific!="IGKV3-20"]
length(common_sig)
length(DC1_specific)
length(DC2_specific)



expression = h5ad$expression
sample_info = h5ad$metadata
sample_info$cell = rownames(sample_info)

### Define expressed genes in sender and receiver cell populations
# here we do not consider sender population as we only have DC cells

cDC3_cDC1_ids = sample_info %>% filter(MajorCluster == 'cDC3-cDC1') %>% pull(cell)
## selected expressed genes (We will consider a gene to be expressed when it is expressed in at least 10% of cells in one cluster.)
expressed_genes_receiver = expression[,cDC3_cDC1_ids] %>% apply(1,function(x){sum(x>0)/length(x)})  %>% .[. >= 0.1] %>% names()

# Check the number of expressed genes: should be a 'reasonable' number of total expressed genes in a cell type, e.g. between 5000-10000 (and not 500 or 20000)
length(expressed_genes_receiver)
## [1] 6351


### Define the gene set of interest and a background of genes
geneset_oi = DC1_specific %>% .[. %in% rownames(ligand_target_matrix)] 

length(geneset_oi)
## [1] "SERPINE1" "TGFBI"    "MMP10"    "LAMC2"    "P4HA2"    "PDPN"

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)
head(background_expressed_genes)
## [1] "RPS11"   "ELMO2"   "PNMA1"   "MMP2"    "TMEM216" "ERCC5"


### Define a set of potential ligands
lr_network = readRDS("lr_network.rds")

# If wanted, users can remove ligand-receptor interactions that were predicted based on protein-protein interactions and only keep ligand-receptor interactions that are described in curated databases. To do this: uncomment following line of code:
# lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = ligands ## here we used all ligands intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

### Perform NicheNet’s ligand activity analysis on the gene set of interest
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

### Infer target genes of top-ranked ligands and visualize in a heatmap
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)


active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

cDC1_specific = vis_ligand_target %>% make_heatmap_ggplot("Pro cDC1 maturation ligands","cDC1-specific maturation signature", color = "red",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "red", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

cDC1_specific








##########cDC2 maturation signature
cDC3_cDC2_ids = sample_info %>% filter(MajorCluster == 'cDC3-cDC2') %>% pull(cell)
expressed_genes_receiver = expression[,cDC3_cDC2_ids] %>% apply(1,function(x){sum(x>0)/length(x)})  %>% .[. >= 0.1] %>% names()

# Check the number of expressed genes: should be a 'reasonable' number of total expressed genes in a cell type, e.g. between 5000-10000 (and not 500 or 20000)
length(expressed_genes_receiver)
## [1] 6351


### Define the gene set of interest and a background of genes
geneset_oi = DC2_specific %>% .[. %in% rownames(ligand_target_matrix)] 
length(geneset_oi)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)

### Define a set of potential ligands
lr_network = readRDS("lr_network.rds")

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = ligands ## here we used all ligands intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

### Perform NicheNet’s ligand activity analysis on the gene set of interest
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

### Infer target genes of top-ranked ligands and visualize in a heatmap
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()


active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets[order_targets %in% rownames(active_ligand_target_links)],order_ligands] %>% t()


cDC2_specific = vis_ligand_target %>% make_heatmap_ggplot("Pro cDC2 maturation ligands","cDC2-specific maturation signature", color = "red",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "red", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

cDC2_specific






##########common maturation signature
cDC3_cDC2_ids = sample_info %>% filter(MajorCluster == 'cDC3-cDC2') %>% pull(cell)
cDC3_cDC1_ids = sample_info %>% filter(MajorCluster == 'cDC3-cDC1') %>% pull(cell)

expressed_genes_receiver = expression[,c(cDC3_cDC2_ids,cDC3_cDC1_ids)] %>% apply(1,function(x){sum(x>0)/length(x)})  %>% .[. >= 0.1] %>% names()

# Check the number of expressed genes: should be a 'reasonable' number of total expressed genes in a cell type, e.g. between 5000-10000 (and not 500 or 20000)
length(expressed_genes_receiver)
## [1] 6351


### Define the gene set of interest and a background of genes
geneset_oi = common_sig %>% .[. %in% rownames(ligand_target_matrix)] 
length(geneset_oi)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)

### Define a set of potential ligands
lr_network = readRDS("lr_network.rds")

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = ligands ## here we used all ligands intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

### Perform NicheNet’s ligand activity analysis on the gene set of interest
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

### Infer target genes of top-ranked ligands and visualize in a heatmap
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()


active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets[order_targets %in% rownames(active_ligand_target_links)],order_ligands] %>% t()

common_res = vis_ligand_target %>% make_heatmap_ggplot("Pro cDC maturation ligands","cDC common maturation signature", color = "red",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "red", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

common_res

pdf(file="Common_maturation.pdf",width=6.99, height=4.96)
common_res
dev.off()

pdf(file="cDC2_specific_maturation.pdf",width=4.28, height=4.96)
cDC2_specific
dev.off()

pdf(file="cDC1_specific_maturation.pdf",width=5.36, height=4.96)
cDC1_specific
dev.off()



