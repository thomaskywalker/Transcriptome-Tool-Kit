plot_venn_heatmap <- function(meta) {
  # import datasets
  REF_GROUP <- meta$data_parameters$REF_GROUP[[meta$data_parameters$REF_GROUP %>% attributes() %>% unlist()]]
  data_normalized <- meta$import_data(paste0("dds_",meta$data_parameters$NORMALIZED_FUNCTION,".rds"))
  coldata <- meta$import_data('col_data.txt')
  resLFC <- meta$import_data('resLFC_all.rds')

  log_info('filter the DEGs with |log2FC| > {meta$data_parameters$lfcThreshold} and padj < {meta$data_parameters$pvalueThreshold}')
  pval_threshold <- meta$data_parameters$pvalueThreshold
  lfc_threshold <- meta$data_parameters$lfcThreshold
  groups_num <- length(resLFC)
  group_gene_id_set <- vector('list',groups_num)
  for(idx in seq_along(resLFC)){
    #group_gene_id_set[[idx]] <- row.names(resLFC[[idx]][which(resLFC[[idx]]$padj < pval_threshold, resLFC[[idx]]$log2FoldChange >= lfc_threshold), ])
    # since have add test to validate lfc>thresholds in null hypothesis, do not need to cut-off lfc as hard thresholds. 
    group_gene_id_set[[idx]] <- row.names(resLFC[[idx]][which(resLFC[[idx]]$padj < pval_threshold), ])
  }
  
  DEG_intersection <- Reduce(intersect, group_gene_id_set)
  
  DEG_unique <-  vector('list',groups_num)
  for (idx in seq_along(resLFC)) {
    unique_idx <- !group_gene_id_set[[idx]] %in% (group_gene_id_set[-idx] %>% unlist())
    DEG_unique[[idx]] <- group_gene_id_set[[idx]][unique_idx]
  }
  DEG_unique <- unlist(DEG_unique)

  pipeline_controller(list(meta, group_gene_id_set, groups_num),
                      (meta$flow_controller$plots_to_draw$venn_diagram),
                      visualize_venn)

  if(meta$flow_controller$plots_to_draw$deg_heatmap)
    deg_heat(meta=meta,data_normalized=data_normalized,genelist=DEG_intersection, coldata=coldata,file_name="intersect")
  # pipeline_controller (list(meta=meta,data_normalized=data_normalized,genelist=DEG_intersection, coldata=coldata,file_name="intersect"),
  #                      meta$flow_controller$plots_to_draw$deg_heatmap,
  #                      deg_heat)

  if(meta$flow_controller$plots_to_draw$deg_heatmap)
    deg_heat(meta=meta,data_normalized=data_normalized,genelist=DEG_unique ,coldata=coldata,file_name="unique")
  #break, spent too much time .cuased by pipeline_contoller , when use if condition everything is fine.
  # pipeline_controller(list(meta=meta,data_normalized=data_normalized,genelist=DEG_unique ,coldata=coldata,file_name="unique"),
  #                      meta$flow_controller$plots_to_draw$deg_heatmap,
  #                      deg_heat)
  
  ggggggggg <<- group_gene_id_set 
  hhhhhhhhh <<- Reduce(cbind,group_gene_id_set) 
  return(NULL)
}

  visualize_venn <- function(meta, group_gene_id_set, groups_num) {
    venn_parameters <- meta$plot_parameters$venn_diagram
    venn_parameters[['x']] <- group_gene_id_set
    venn_parameters[['resolution']] <- 600
    venn_parameters[['category.names']] <- meta$group_list
    venn_parameters[['cat.dist']] <- rep(0.1,groups_num)
    if(is.null(venn_parameters[['fill']])){
      venn_parameters[['fill']] = rainbow(groups_num)
    }
    if(is.null(venn_parameters[['filename']])){
      venn_parameters[['filename']] <-meta$get_file_path('Plot',"VennDiagram.png")
    }else{
      venn_parameters[['filename']] <-meta$get_file_path('Plot',venn_parameters[['filename']])
    }
    do.call('venn.diagram',venn_parameters)
    log_info("Venn diagram saved")
  }
  deg_heat <- function(meta,data_normalized, genelist, coldata, file_name) {
      deg_heatmap_parameters <- meta$plot_parameters$deg_heatmap
      deg_heatmap_parameters[['mat']] <- data_normalized[genelist, ]
      if (is.null(deg_heatmap_parameters[['main']]))
        deg_heatmap_parameters[['main']] <- file_name
      deg_heatmap_parameters[['annotation_col']] <- coldata
      heat_deg <- do.call('pheatmap',deg_heatmap_parameters)
      # save files
      deg_heat_path <- meta$get_file_path('Plot',paste0("heat_",file_name,"_DEG.png"))
      ggsave(filename = deg_heat_path, heat_deg, height = 10, width = 8,dpi=600)
      log_info('{deg_heat_path} saved')
  }
