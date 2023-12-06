plot_heatmap <- function(meta) {
  normalized_data <-meta$import_data(paste0('dds_',meta$data_parameters$NORMALIZED_FUNCTION,'.rds'))
  coldata <- meta$import_data('col_data.txt')
  rownames(normalized_data) <- rownames(normalized_data) %>%map_key_to_col(.,key=meta$data_parameters$origin_seq_id,column= meta$data_parameters$heatmap_id,db=meta$DB)
  # 4.2 heatmap
  # select top 30 expressed genes for heatmap
  pipeline_controller(list(meta,normalized_data),
                      meta$flow_controller$plots_to_draw$expression_heatmap,
                      timer(visualize_expression_heatmap))
  # 4.3 correlation matrix
  pipeline_controller(list(meta, normalized_data),meta$flow_controller$plots_to_draw$cor_plot,timer(visualize_cor_plot))
  # 4.3.1 Sample distance by euclidean
  pipeline_controller(list(meta, normalized_data),meta$flow_controller$plots_to_draw$sample_cluster,timer(visualize_sample_cluster))
}

visualize_expression_heatmap <- function(meta,normalized_data) {
  expression_heatmap_parameters <- get_heatmap_parameters( meta,normalized_data)
  plot_expression_heatmap(meta=meta,expression_heatmap_parameters=expression_heatmap_parameters)
  return(NULL)
}

get_heatmap_parameters <- function(meta,normalized_data ) {
  coldata <- meta$import_data('col_data.txt')
  selection <-
    order(rowMeans(normalized_data), decreasing = TRUE)[1:meta$plot_parameters$expression_heatmap$HEATMAP_GENE_NUM]
  expression_heatmap_parameters <- meta$plot_parameters$expression_heatmap
  expression_heatmap_parameters[['mat']] <- normalized_data[selection,]
  expression_heatmap_parameters[['annotation_col']] <- coldata
  return(expression_heatmap_parameters)
}



plot_expression_heatmap <- function(meta,expression_heatmap_parameters) {
  expression_heatmap_parameters <- expression_heatmap_parameters[-match(c('HEATMAP_GENE_NUM'),names(expression_heatmap_parameters))]
  heat <- do.call('pheatmap',expression_heatmap_parameters)
  heatmap_path <-
    meta$get_file_path( "Plot", paste0("top",meta$plot_parameters$heatmap$HEATMAP_GENE_NUM,"_heatmap_row.png"))
  ggsave(filename = heatmap_path ,
         heat,
         dpi=600,
         height = 14,
         width = 10)
  return(NULL)
}

visualize_cor_plot <- function(meta, normalized_data) {
  # 4.3.2 spearman correlation
  # spearman correlation is suitable for small sample sizes
  # pearson is for normal distributed samples
  cor_plot_parameters <- get_cor_plot_parameters(meta, normalized_data)
  plot_cor_plot(cor_plot_parameters, meta)
}
get_cor_plot_parameters <- function(meta, normalized_data) {
  coldata <- meta$import_data('col_data.txt')
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  cor_plot_parameters <- meta$plot_parameters$cor_plot
  sampleCor <- cor(normalized_data, method =meta$plot_parameters$cor_plot$cor_method)
  cor_plot_parameters[['mat']] <-   sampleCor
  cor_plot_parameters[['col']] <-   colors
  cor_plot_parameters[['annotation_col']] <-   coldata
  return(cor_plot_parameters)
}
plot_cor_plot <- function(cor_plot_parameters, meta) {
  corr <- do.call('pheatmap',cor_plot_parameters[-match(c('cor_method'),names(cor_plot_parameters))])
  correlation_path <-
    meta$get_file_path( "Plot", "correlation.png")
  ggsave(filename = correlation_path,
         corr,
         dpi=600,
         height = 8,
         width = 10)
}

visualize_sample_cluster <- function(meta,normalized_data) {
  sample_cluster_parameters <- get_sample_cluster_parameters(meta,normalized_data)
  plot_sample_distance(meta,sample_cluster_parameters)
}
get_sample_cluster_parameters <- function(meta, normalized_data) {
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  coldata <- meta$import_data('col_data.txt')
  sampleDistMatrix <- as.matrix(dist(t(normalized_data), method = meta$plot_parameters$sample_cluster$normalized_method))
  sample_cluster_parameters <- meta$plot_parameters$sample_cluster
  sample_cluster_parameters[['col']] <- colors
  sample_cluster_parameters[['annotation_col']] <- coldata
  sample_cluster_parameters[['mat']] <- sampleDistMatrix
  sample_cluster_parameters[['cellwidth']] <-  25
  sample_cluster_parameters[['cellheight']] <-  25
  sample_cluster_parameters[['display_numbers']] <-  TRUE
  return(sample_cluster_parameters)
}
plot_sample_distance <- function(meta,sample_cluster_parameters) {
  dist <- do.call(pheatmap,sample_cluster_parameters)
  sample_distance_path <-  meta$get_file_path( "Plot", "sample_distance.png")
  ggsave(filename = sample_distance_path,
         dist,
         dpi=600,
         height = 8,
         width = 10)
}
