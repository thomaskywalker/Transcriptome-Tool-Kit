plot_pvalue_LFC <- function(meta) {

  groups_list <- meta$import_data('deseq2_dds.rds') %>% resultsNames(.) %>% .[-1]
  resLFC_all <- meta$import_data('resLFC_all.rds') %>% tranform_LFC_id(meta,groups_list,.)
  pipeline_controller(list(groups_list,meta, resLFC_all),
                      meta$flow_controller$plots_to_draw$MA_plot,
                      plot_MA
                      )
  pipeline_controller(list(groups_list,meta, resLFC_all),
                      meta$flow_controller$plots_to_draw$volcano_plot,
                      plot_volcano
                      )

}

tranform_LFC_id <- function(meta, groups_list,resLFC_all) {
  for (idx in seq_along(groups_list)) {
    rownames(resLFC_all[[idx]]) <- map_key_to_col(rownames(resLFC_all[[idx]]),key=meta$data_parameters$origin_seq_id,column= meta$data_parameters$tranform_LFC_id,meta$DB)
  }
  return(resLFC_all)
}

get_volcano_plot_parameters <- function(meta,data) {
  # 3 Exploring and exporting results
  # 3.2 Volcano plot
  # use EnhancedVolcano function
  volcano_plot_parameters <- meta$plot_parameters$volcano_plot
  volcano_plot_parameters[['toptable']] <- data
  volcano_plot_parameters[['lab']] <- rownames(data)
  volcano_plot_parameters[['y']] <- "pvalue"
  volcano_plot_parameters[['x']] <- "log2FoldChange"
   volcano_plot_parameters[['col']] <-  c("grey30", "#2a2ac2b5", "royalblue", "#f1a71d")
  volcano_plot_parameters[['legendLabels']] <- c("Not sig.", "Log2FC", "p-value", "p-value & Log2FC")
  return(volcano_plot_parameters)
}

get_MA_parameters <- function(meta, data) {
  MA_plot_parameters <- meta$plot_parameters$MA_plot
  MA_plot_parameters[['ylim']] <-  c(-10, 10)
  MA_plot_parameters[['object']] <- data
  return(MA_plot_parameters)
}

visualize_MA_plot <- function(){

}
plot_MA <- function(groups_list, meta, resLFC_all) {
  for (idx in seq_along(groups_list)) {
    group_name <- groups_list[[idx]]
    png(width = 10, height = 10, units = "in",filename=meta$get_file_path('Plot', paste0(group_name,"_MA_plot.png")),res=600)
    MA_plot_parameters <- get_MA_parameters(meta, data=resLFC_all[[idx]])
    maplot <- do.call('plotMA',MA_plot_parameters)
    dev.off()
  }
}
plot_volcano <- function(groups_list, meta, data) {
  for (idx in seq_along(groups_list)) {
    group_name <- groups_list[[idx]]
    volcano_plot_parameters <- get_volcano_plot_parameters(meta,data=data[[idx]])
    volcano <- do.call('EnhancedVolcano',volcano_plot_parameters)
    ggsave(
      filename = meta$get_file_path( 'Plot', paste0(group_name,"_volcano_plot.png")),
      volcano,
      dpi=600,
      height = 10,
      width = 10
    )
  }
}
