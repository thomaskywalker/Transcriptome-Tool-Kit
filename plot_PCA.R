plot_PCA <- function(meta) {
  log_info('start draw PCA plot')
  # 6 PCA analysis
  normalized_data <- meta$import_data(paste0("dds_",meta$data_parameters$NORMALIZED_FUNCTION,".rds"))
  coldata<- meta$import_data('col_data.txt')
  rownames(normalized_data) <- rownames(normalized_data) %>%map_key_to_col(.,key=meta$data_parameters$origin_seq_id,column= meta$data_parameters$PCA_id,db=meta$DB)
  ## removing the lower 10% of variables based on variance
  p <- PCAtools::pca(normalized_data, metadata = coldata, removeVar = 0.1)
  # 6.1 PCA biplot
  pipeline_controller(list(meta, p),
                      meta$flow_controller$plots_to_draw$biplot,
                      visualize_to_biplot
  )
  pipeline_controller(list(meta, p),
                      meta$flow_controller$plots_to_draw$scree_plot,
                      visualize_to_screeplot
  )

  # 6.3 pairs plot
  # The number of PCs to visualise is determined in scree plot
  pipeline_controller(list(meta, p),
                      meta$flow_controller$plots_to_draw$pairsplot,
                      visualize_to_pairsplot
  )

  # 6.4 loadings plot
  pipeline_controller(list(meta, p),
                      meta$flow_controller$plots_to_draw$loading_plot,
                      visualize_to_loading_plot
  )

  return(NULL)
}
visualize_to_loading_plot <-  function(meta, p) {
  loading_plot_parameters <- get_loading_plot_parameters(meta, p)
  plot_loading_plot(loading_plot_parameters, meta)
}
visualize_to_pairsplot <-   function(meta, p) {
  pairsplot_parameters <- get_pairsplot_parameters(meta, p)
  plot_pairsplot(pairsplot_parameters, meta)
}
visualize_to_biplot <- function(meta, p) {
  biplot_parameters <- get_biplot_parameters(meta, p)
  plot_biplot(biplot_parameters, meta)}

visualize_to_screeplot <- function(meta, p) {
  scree_plot_parameters <- get_scree_parameters(meta, p)
  plot_scree_plot(scree_plot_parameters, meta)
}
get_biplot_parameters <- function(meta, p) {
  biplot_parameters <- meta$plot_parameters$biplot
  biplot_parameters[['pcaobj']] <- p
  biplot_parameters[['vline']] <- c(-25, 0, 25)
  biplot_parameters[['colby']] <- 'group'
  return(biplot_parameters)
}

get_scree_parameters <- function( meta,p ) {
  horn <- parallelPCA(meta$import_data(paste0("dds_",meta$data_parameters$NORMALIZED_FUNCTION,".rds"))) # Horn’s parallel analysis
  elbow <- findElbowPoint(p$variance) # Elbow method
  scree_plot_parameters <- meta$plot_parameters$scree_plot
  scree_plot_parameters[['pcaobj']] <- p
  scree_plot_parameters[['components']] <- getComponents(p, 1:15)
  scree_plot_parameters[['vline']] <- c(horn$n, elbow)
  return(scree_plot_parameters)
}

get_pairsplot_parameters <- function(meta, p) {
  horn <- parallelPCA(meta$import_data(paste0("dds_",meta$data_parameters$NORMALIZED_FUNCTION,".rds"))) # Horn’s parallel analysis
  pairsplot_parameters <- meta$plot_parameters$pairsplot
  pairsplot_parameters[['pcaobj']] <- p
  pairsplot_parameters[['num_components']] <- horn$n
  pairsplot_parameters[['components']] <- getComponents(p, seq(horn$n))
  pairsplot_parameters[['margingaps']] <- unit(c(0.1, 0.1, 0.1, 0.1), "in")
  pairsplot_parameters[['colby']] <- 'group'
  return(pairsplot_parameters)
}

get_loading_plot_parameters <- function(meta, p) {
  loading_plot_parameters <- meta$plot_parameters$loading_plot
  loading_plot_parameters[['pcaobj']] <- p
  return(loading_plot_parameters)
}

plot_biplot <- function(biplot_parameters, meta) {
  bi  <- do.call('biplot',biplot_parameters)
  ggsave(
    filename = meta$get_file_path("Plot", "PCA_biplot.png"),
    bi,
    dpi=600,
    height = 10,
    width = 10
  )
  log_info("PCA_biplot.png saved")
}

plot_scree_plot <- function(scree_plot_parameters, meta) {
  scree <- do.call('screeplot',scree_plot_parameters) +
    geom_label(aes(
      x = scree_plot_parameters['vline'][[1]][[1]] - 0.5,
      y = 40,
      label = "Horn's method",
      vjust = -1,
      size = 8
    )) +
    geom_label(aes(
      x = scree_plot_parameters['vline'][[1]][[2]] + 0.5,
      y = 50,
      label = "Elbow method",
      vjust = -1,
      size = 8
    ))
  ggsave(
    filename = meta$get_file_path( "Plot", "PCA_screeplot.png"),
    scree,
    dpi=600,
    height = 10,
    width = 10
  )
  log_info("PCA_screeplot.png saved")
}

plot_loading_plot <- function(loading_plot_parameters, meta) {
  load <- do.call('plotloadings',loading_plot_parameters)
  ggsave(
    filename = meta$get_file_path( "Plot", "PCA_loadingsplot.png"),
    load,
    dpi=600,
    height = 10,
    width = 10
  )
  log_info("PCA_loadingsplot.png saved")
}

plot_pairsplot <- function(pairsplot_parameters, meta) {
  pair <- do.call('pairsplot',pairsplot_parameters[-match(c('num_components'),names(pairsplot_parameters))])
  ggsave(
    filename = meta$get_file_path( "Plot", "PCA_pairsplot.png"),
    pair,
    dpi=600,
    height = 10,
    width = 10,
    bg = "white"
  )
  log_info("PCA_pairsplot.png saved")
}
