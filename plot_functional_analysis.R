plot_functional_analysis <- function(meta) {
  # 0.2 load functional enrichment results
  coldata <- meta$import_data('col_data.txt')
  all_fe <- meta$import_data('functional_all.rds')
  subcat <- names(all_fe[[1]])

  for (idx in seq_along(meta$group_list)) {
    input <- all_fe[[idx]]
    fname <- meta$group_list[idx]
    log_debug('{fname} start')
    folder_name <- paste0(fname, "_fe_visual")
    dir.create(meta$get_file_path('Plot', folder_name), showWarnings = FALSE)
    chosen_plot <- meta$flow_controller$plots_to_draw
    for (category in subcat) {
      dat <- input[[category]][['results']]
      pipeline_controller(list(meta,dat, fname, category, folder_name),chosen_plot$dot_plot & (!is.enrichResult(dat)),plot_Normal_dotplot)
      pipeline_controller(list(meta,dat, fname, category,  folder_name),chosen_plot$dot_plot & (is.enrichResult(dat)),plot_enrich_dotplot)
      pipeline_controller(list(meta, dat, fname, category, folder_name),(chosen_plot$ridge_plot & (!is.enrichResult(dat))),plot_ridge_plot)
      #ORA ridge break whn use pipeline_controller, spent too much time (>20 mins).cuased by pipeline_contoller , when use if condition everything is fine.
      if((chosen_plot$ridge_plot & (is.enrichResult(dat))))
        plot_enrich_emapplot(dat, meta, folder_name, fname, category)
    }
  }
  return(NULL)
}

is.enrichResult <- function(obj) {obj %>% class(.)=="enrichResult"}

plot_Normal_dotplot <- function(meta,dat, fname, category, folder_name) {
  try({
    log_debug('{fname}:{category} dotplot start')
    dot <-
      enrichplot::dotplot(dat, showCategory = 10, split = ".sign") +
      facet_grid(. ~ ".sign") +
      ggplot2::ggtitle(paste(fname, "_", as.character(category))) +
      theme(plot.title = element_text(size = 15, hjust = 0.5))
    ggsave(
      filename = meta$get_file_path(
        'Plot',
        paste0(folder_name, '/' , fname, "_", category, "_dot.png")
      ),
      dot,
      dpi=600,
      height = 12,
      width = 10
    )
    log_info('{fname}:{category} dotplot saved')
  })
}

plot_enrich_dotplot <- function(meta,dat, fname, category,  folder_name) {
  try({
    log_debug('{fname}:{category} dotplot start')
    dot <- enrichplot::dotplot(dat, showCategory = 10) +
      ggplot2::ggtitle(paste(fname, "_", as.character(category))) +
      theme(plot.title = element_text(size = 15, hjust = 0.5))
    ggsave(
      filename = meta$get_file_path(
        'Plot',
        paste0(folder_name, '/' , fname, "_", category, "_dot.png")
      ),
      dot,
      dpi=600,
      height = 12,
      width = 10
    )
    log_info('{fname}:{category} dotplot saved')
  })
}

plot_ridge_plot <- function(meta, dat, fname, category, folder_name) {
  # ORA can't run ridge plot
  try({
    log_debug('{fname}:{category} ridgeplot start')
    ridge <- enrichplot::ridgeplot(dat, showCategory = 20) +
      ggplot2::ggtitle(paste(fname, "_", as.character(category))) +
      theme(plot.title = element_text(size = 15, hjust = 0.5))
    ggsave(
      filename = meta$get_file_path(
        'Plot',
        paste0(folder_name, '/' , fname, "_", category, "_ridge.png")
      ),
      ridge,
      dpi=600,
      height = 12,
      width = 10
    )
    log_info('{fname}:{category} ridgeplot saved')
  })
}

plot_enrich_emapplot <- function(dat, meta, folder_name, fname, category) {
  # ORA can't run ridge plot
  try({
    log_debug('{fname}:{category} emapplot start')
    p <- emapplot(pairwise_termsim(dat))
    ggsave(meta$get_file_path(
      'Plot',
      paste0(folder_name, '/' , fname, "_", category, "_emap.png"))
      ,p,dpi=600)
    log_info('{fname}:{category} emapplot saved')
  })

  log_info('{fname}:{category} plot saved')
}
