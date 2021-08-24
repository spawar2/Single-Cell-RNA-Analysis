library(tidyverse)
library(Seurat)
library(gridExtra)
library(phateR)
library(sctransform)
# library(slingshot)
library(pheatmap)
library(patchwork)
library(MAST)


set.seed(123)

get_cell_numbers <- function(dl){
  my_d <- Read10X(dl)
  lapply(colnames(my_d), function(x) unlist(strsplit(x, "_"))[1]) %>%
    unlist() %>% table() %>% print()
}

my.make.exp <- function(dl, project.name, species = "Human", downsample = NA, op_dir = "./"){
  my.d <- Read10X(dl) # dl is a named directory list
  if(!is.na(downsample) && downsample > 0){
    snames <- unlist(lapply(colnames(my.d), function(x) unlist(strsplit(x, "_"))[1]))
    inds <- split(1:length(snames), snames)
    set.seed(123)
    use_inds <- lapply(names(inds), function(i){
      s <- sample(inds[[i]], downsample, replace = F)
    })
    use_inds <- unlist(use_inds)
    my.d <- my.d[, use_inds]
  }
  
  my.exp <- CreateSeuratObject(my.d, project = project.name, min.cells = 3, min.features = 200)
  
  pdfs <- paste(file.path(op_dir, "QC_"), 1:2, ".pdf", sep = "")
  
  if(species == "Human"){
    my.exp[["percent.mito"]] <- PercentageFeatureSet(object = my.exp, pattern = "^MT-")
  } else {
    my.exp[["percent.mito"]] <- PercentageFeatureSet(object = my.exp, pattern = "^mt-")
  }
  
  VlnPlot(my.exp, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
  ggsave(pdfs[1])
  
  plot1 <- FeatureScatter(my.exp, feature1 = "nCount_RNA", feature2 = "percent.mito")
  plot2 <- FeatureScatter(my.exp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  ggsave(pdfs[2], width = 14, height = 7, units = "in")
  
  my.exp
}

filter.normalize.scale <- function(my.exp, cellcycle = FALSE, op_dir = "./"){ # }, mito.cutoff = 5, feature.cutoff = 2500){
  pdfs <- paste(file.path(op_dir, "QC_"), 3:4, ".pdf", sep = "")
  # mito.cutoff <- mito.cutoff
  # feature.cutoff <- feature.cutoff
  # my.exp <- subset(my.exp, subset = nFeature_RNA > 200 & nFeature_RNA < feature.cutoff & percent.mito < mito.cutoff)
  
  # my.exp <- NormalizeData(my.exp, normalization.method = "LogNormalize", scale.factor = 10000)
  if(cellcycle){
    s_genes <- cc.genes$s.genes
    g2m_genes <- cc.genes$g2m.genes
    my.exp <- CellCycleScoring(my.exp, s.features = s_genes, g2m.features = g2m_genes)
    my.exp <- SCTransform(my.exp, vars.to.regress = c("S.Score", "G2M.Score", "percent.mito", "nCount_RNA"))
  } else {
    my.exp <- SCTransform(my.exp, vars.to.regress = c("percent.mito", "nCount_RNA"))
  }
  
  my.exp <- FindVariableFeatures(my.exp, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(my.exp), 10)
  print(top10) # for debugging
  
  plot1 <- VariableFeaturePlot(my.exp)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
  CombinePlots(plots = list(plot1, plot2))
  ggsave(pdfs[1], width = 14, height = 7, units = "in")
  saveRDS(my.exp, file.path(op_dir, "experiment.RDS"))
  
  all.genes <- rownames(my.exp)
  print(length(all.genes)) # for debugging
  # my.exp <- ScaleData(my.exp, features = all.genes, vars.to.regress = c("percent.mito", "nCount_RNA"))
  # my.exp <- ScaleData(my.exp, features = all.genes, vars.to.regress = c("percent.mito", "nCount_RNA"))
  # my.exp <- ScaleData(my.exp)
  my.exp
}

reduce.dimensions <- function(my.exp, op_dir = "./", ...){
  pdfs <- file.path(op_dir, paste("plot_", 1:4, ".pdf", sep = ""))
  my.exp <- RunPCA(object = my.exp, features = VariableFeatures(my.exp))
  saveRDS(my.exp, file.path(op_dir, "experiment-clustered.RDS"))
  
  VizDimLoadings(my.exp, dims = 1:2, reduction = "pca")
  
  ggsave(pdfs[1])
  
  my.exp <- FindNeighbors(my.exp, dims = 1:30)
  saveRDS(my.exp, file.path(op_dir, "experiment-clustered.RDS"))
  my.exp <- FindClusters(my.exp, resolution = 0.5)
  saveRDS(my.exp, file.path(op_dir, "experiment-clustered.RDS"))
  my.exp <- RunTSNE(my.exp, dims = 1:30, seed.use = 123, ...)
  my.exp <- RunUMAP(my.exp, dims = 1:30, method = "uwot")
  saveRDS(my.exp, file.path(op_dir, "experiment-clustered.RDS"))
  print("TSNE worked")
  ##  tryCatch(my.exp <- RunUMAP(my.exp, dims = 1:30),
  ##           error = function(e) {
  ##                print("UMAP not installed")
  ##	           },
  ##           finally = {})
  p1 <- DimPlot(my.exp, label = T, label.size = 6)
  p2 <- DimPlot(my.exp, group.by = "orig.ident")
  p <- CombinePlots(plots = list(p1, p2), ncol = 2)
  ggsave(pdfs[2], width = 14, height = 7)
  my.exp
}

use.phate <- function(my.exp){
  ## Use PHATE on the data.
  my.exp.phate <- phateR::phate(t(GetAssayData(my.exp, slot = "scale.data")), knn = 20, gamma = 0)
  my.exp[["phate"]] <- CreateDimReducObject(embeddings = my.exp.phate$embedding, key = "PHATE_", assay = DefaultAssay(my.exp))
  my.exp
}

get.species.cells <- function(fn, name, sp = "hg38"){
  ## Select all the cells that belong to a single species from the gem.txt file.
  df <- read.delim(fn, sep = ",", stringsAsFactors = F, header = T)
  df %>% filter(call == sp) %>% select(barcode) %>% transmute(barcode = gsub("-1", "",  barcode)) %>%
    unlist() %>% as.vector() %>% paste(name, ., sep = "_")
}

make_analysis_from_table <- function(fn){
  df_exp <- read.delim(fn, header = T, stringsAsFactors = F, sep = "\t")
  # we are expecting a tab-separated-values file. We expect each line has five
  # columns, first column denotes name of the analysis folder, second column
  # denotes comma-separated list of directores (full paths) to be used in the
  # analysis, and third column is a comma separated list of labels to be used
  # for each of the analyses.  The fourth column is the species of the
  # reference.  The last column should be a unique project ID so we can merge
  # the analyses later if required.  Currently this analysis assumes that the
  # gene expression is going to be corrected for cell-cycle gene expression and
  # for marker genes we are going to use the repressed and overexpressed genes
  # for each cluster.
  for(i in 1:nrow(df_exp)){
    op_dir <- df_exp[i, 1]
    dl     <- unlist(strsplit(df_exp[i, 2], ","))
    names(dl) <- unlist(strsplit(df_exp[i, 3], ","))
    print(dl) # for debugging
    use_species <- df_exp[i, 4]
    project_name <- df_exp[i, 5]
    message(paste("Starting analysis ...", project_name))
    message(paste(op_dir, dl, use_species, project_name, sep = "\n"))
    message(paste("Creating output directory ", op_dir))
    system(paste("mkdir -p ", op_dir))
    my_exp <- make.analysis(dl, op_dir = op_dir, project.name = project_name,
                            species = use_species, cellcycle = TRUE, 
                            only.pos = FALSE)
    saveRDS(my_exp, paste(op_dir, "experiment-clustered.RDS", sep = ""))
    message(paste("Finished analysis ", project_name, "."))
  }
}

make.analysis <- function(dl, op_dir = "./", project.name = "SCRNAseq", species = "Human", 
                          downsample = NA, cellcycle = FALSE, only.pos = T){
  ## Complete the analysis using normal defaults, generate graphs and also generate differential gene expression profile.
    op_dir = dl
  my.exp <- my.make.exp(dl, project.name = project.name, species = species, downsample = downsample, op_dir = op_dir)
  my.exp <- filter.normalize.scale(my.exp, cellcycle = cellcycle, op_dir = op_dir)
  my.exp <- reduce.dimensions(my.exp, op_dir = op_dir)
  my.exp <- use.phate(my.exp)
  all.markers <- FindAllMarkers(my.exp, only.pos = only.pos, min.pct = 0.25, test.use = "MAST")
  write.table(all.markers, paste(op_dir, "all-cluster-markers.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")
  avg.exp <- AverageExpression(my.exp, add.ident = "orig.ident")
  write.table(avg.exp[["SCT"]], paste(op_dir, "average-expression-per-cluster-per-condition.txt", sep = ""), row.names = T, col.names = T, quote = F, sep = "\t")
  all_p <- make_all_graphs(my.exp)
  ## my.exp@meta.data %>% select(orig.ident, seurat_clusters) %>% table() %>%
  ##    barplot(beside = T, legend.text = T, args.legend = list(bty = "n"), main = "Cell number by Clusters", ylab = "Frequency", xlab = "Cluster ID")
  ggsave(device = "pdf", height = 21, width = 11, units = "in", dpi = 300, filename = paste(op_dir, "number-of-cells-by-cluster.pdf", sep = ""))
  my.exp
}

make_bp_sample_cluster <- function(my_exp){
  ## Make a barplot
  my_exp@meta.data %>%  dplyr::select("orig.ident", "seurat_clusters") %>% table() %>%
    as.data.frame() %>% ggplot(aes_string(x = "seurat_clusters", y = "Freq",
                                          fill = "orig.ident", group = "orig.ident")) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_bw() + scale_fill_discrete(name = "Sample") +
    ggtitle("Cell Numbers by cluster by sample") +
    xlab("Cluster ID") + ylab("Cell Numbers") -> bp_cellnumbers
}

make_bp_prop <- function(my_exp){
  my_exp@meta.data %>% dplyr::select("orig.ident", "seurat_clusters") %>%
    table() %>% prop.table(2) %>% as.data.frame() %>%
    ggplot(aes_string(x = "seurat_clusters", y = "Freq",
                      fill = "orig.ident", group = "orig.ident")) +
    geom_bar(stat = "identity") +
    theme_minimal() + scale_fill_discrete("Sample") +
    ggtitle("Cell Proportions by Cluster by Sample") +
    xlab("Cluster ID") + ylab("Proportion") -> bp_proportion
}

make_all_graphs <- function(my_exp){
  ## make all the graphs that we want to for the given data-set.  Start with the number of cells
  ## per sample per cluster.
  bpn <- make_bp_sample_cluster(my_exp)
  bpp <- make_bp_prop(my_exp)
  t_p1  <- DimPlot(my_exp, reduction = "tsne", label = T, label.size = 4) + NoLegend()
  t_p2  <- DimPlot(my_exp, reduction = "tsne", group.by = "orig.ident")
  u_p1 <- DimPlot(my_exp, reduction = "umap", label = T, label.size = 4) + NoLegend()
  u_p2 <- DimPlot(my_exp, reduction = "umap", group.by = "orig.ident")
  p_p1 <- DimPlot(my_exp, reduction = "phate", label = T, label.size = 4) + NoLegend()
  p_p2 <- DimPlot(my_exp, reduction = "phate", group.by = "orig.ident")
  r1 <- t_p1 + t_p2
  r2 <- u_p1 + u_p2
  r3 <- p_p1 + p_p2
  all_plot <- (bpn / bpp / r1 / r2 / r3) + plot_layout(guides = "collect")
  all_plot
}

make_dot_plot_top_n <- function(my_exp, n = 2){
  ## Create a dot plot for the top 2 genes per cluster.
  all_markers <- FindAllMarkers(my_exp, test.use = "MAST", only.pos = T,
                                random.seed = 1234, min.pct = 0.25)
  top_n <- all_markers %>% group_by(cluster) %>% top_n(n, avg_logFoldChange)
  p <- DotPlot(my_exp, unique(top_n$gene), cols = c("lightgrey", "red"))
  op.l <- list(all_markers = all_markers, dotplot = p)
  op.l
}

save_markers <- function(my.exp, fn){
  all_markers <- FindAllMarkers(my.exp, min.pct = 0.25, only.pos = T, test.use = "MAST")
  write.table(all_markers, fn, row.names = T, col.names = T, quote = F, sep = "\t")
  all_markers
}

save_per_cluster_expression <- function(my.exp, fn){
  ## save the average expression values per cluster per sample.
  avg_exp <- AverageExpression(my.exp, add.ident = "orig.ident")
  write.table(avg_exp["SCT"], fn, row.names = T, col.names = T, quote = F, sep = "\t")
  avg_exp
}

runSlignshot <- function(object, reduction = "phate", group.by = NULL,
                         start.clus = NULL, end.clus = NULL,
                         approx_points = FALSE, allow.breaks = TRUE){
  rd <- Embeddings(object, reduction)
  cl <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  cl <- FetchData(object, vars = group.by) %>% pull(`group.by`)
  
  object@misc[["sds"]] <- list("dr" = reduction, "data" = slingshot(rd, cl))
  ps <- slingPseudotime(object@misc[['sds']]$data)
  object@meta.data[, colnames(ps)] <- as.data.frame(ps)
  object <- LogSeuratCommand(object = object)
  return(object)
}

runPseudoTimeDGE <- function(object){
  var_genes <- VariableFeatures(object)
  
  DGE <- list()
  for(c_name in names(object@misc$sds$data@curves)){
    object@misc$sds$dge[[c_name]] <- FetchData(object, append(var_genes, c_name, 0)) %>%
      tidyr::gather(gene, signal, -one_of(c_name)) %>%
      dplyr::rename(curve = 1) %>%
      tidyr::nest(-gene) %>%
      mutate(
        fit = purrr::map(data, ~ gam::gam(signal ~ gam::lo(curve), data = .x)),
        tidied = purrr::map(fit, broom::tidy)
      ) %>%
      tidyr::unnest(tidied) %>%
      dplyr::filter(term != 'Residuals')
  }
  
  object <- LogSeuratCommand(object = object)
  return(object)
}

plotPseudoTime <- function(object, group.by = NULL, reduction = "phate", dims = 1:2){
  object[["ident"]] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  dims <- paste0(Key(object = object[[reduction]]), dims)
  
  curved <- bind_rows(lapply(names(object@misc$sds$data@curves),
                             function(x) {
                               my_c <- slingCurves(object@misc$sds$data)[[x]]
                               my_df <- as.data.frame(my_c$s[my_c$ord, dims])
                               my_df$curve <- x
                               return(my_df)
                             }))
  data <- FetchData(object = object, vars = c(dims, group.by))
  
  p <- ggplot(data, aes_string(x = dims[1], y = dims[2])) +
    geom_point(aes(color = !!sym(group.by))) + theme_bw() +
    theme(legend.position = "top")
  
  if(is.character(data[[group.by]]) || is.factor(data[[group.by]])){
    p <- p + guides(col = guide_legend(nrow = 2))
  } else {
    p <- p + scale_color_distiller(palette = "RdYlBu", na.value = 'grey90')
  }
  p <- p + geom_path(aes_string(dims[1], dims[2], linetype = "curve"), curved, size = 1)
  
  p
}

plotCurveHeatmaps <- function(object = NULL, curve = NULL, n = 25){
  cells <- object@meta.data %>% tibble::rownames_to_column(var = 'cellid') %>%
    dplyr::arrange(!!sym(curve)) %>% dplyr::filter(!is.na(!!sym(curve)))
  genes <- object@misc$sds$dge[[curve]] %>% dplyr::arange(p.value) %>%
    head(n) %>% pull(gene)
  FetchData(object = object, vars = genes, cells = cells$cellid) %>%
    t(.) %>%
    pheatmap(., scale = "row", cluster_cols = F, )}