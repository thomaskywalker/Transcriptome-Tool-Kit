# the script load all of function and libraries we need.
# load libraries ----
library(renv)
library(logger)
log_layout(layout_glue_colors)
log_info('Start loading Packages')
library(tidyverse)
library(magrittr)
library(AnnotationHub)
library(AnnotationDbi)
library(clusterProfiler) # GO, KEGG and GSEA analyses
library(msigdbr) # gene set database for GSEA
library(DOSE) # DO, NCG analyses
library(ReactomePA) # reactome analysis
library(pasilla)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggalt)
library(RColorBrewer)
library(pheatmap)
library(VennDiagram)
library(EnhancedVolcano)
library(PCAtools)
library(glue)
library(assertthat)
library(BiocParallel)
library(enrichplot)
library(pathview)
library(fgsea)
library(ggridges)
library(yaml)
library(foreach)
library(doParallel)


# load function in other R scripts ----
Rfile_pre <- list.files(pattern='\\.R$')
Rfile <- Rfile_pre[-match(c('RUN.R','dependency.R'),table=Rfile_pre)]
lapply(Rfile,source)

# load function ----
map_key_to_col <- function(gene_list,key='ENSEMBL',column= "SYMBOL",db) {
  log_debug('Map the {key} gene id to {column} for visualisation')
  col <- gene_list %>%
    mapIds(db,
           keys = .,
           column = column,
           keytype = key) %>%
    ifelse(is.na(.) | duplicated(.), names(.), .)
  return(col)
}


create_folder <- function(path,folder_name) {
  if(folder_name %in% dir(path)){
    return(NULL)
  }else{
    dir.create(file.path(path,folder_name), showWarnings = FALSE)
    log_debug('Folder {folder_name} Created.')
  }
}

read_data <- function(file_path){
  identify_file_type <- function(file_path) {
    # 獲取文件擴展名並轉換為小寫
    file_extension <- tolower(tools::file_ext(file_path))

    if (file_extension == "csv") {
      return("csv")
    } else if (file_extension == "txt") {
      return("txt")
    } else if (file_extension == "rds") {
      return("rds")
    } else {
      # 如果擴展名無法確定，嘗試讀取文件
      tryCatch({
        readRDS(file_path)
        return("rds")
      }, error = function(e) {
        return("unknown")
      })
    }
  }
  type <- identify_file_type(file_path)
  if(type == 'csv')
    return(read.csv(file_path))
  if(type == 'txt')
    return(read.table(file_path))
  if(type == 'rds')
    return(readRDS(file_path))
}



flow_control <- function(meta){
    pipeline_controller(list(meta, data_dds_transformation, data_normalized, data_functional_analyses),
                        meta$flow_controller$data_pipe_line$Ctrl,
                        timer(data_transform_steps))


    pipeline_controller(list(meta, plot_heatmap, plot_pvalue_LFC, plot_PCA, plot_venn_heatmap, plot_functional_analysis),
                        meta$flow_controller$plots_to_draw$Ctrl,
                        timer(plot_steps))
}

plot_steps <- function(meta, plot_heatmap, plot_pvalue_LFC, plot_PCA, plot_venn_heatmap, plot_functional_analysis) {
  parameters <- meta$flow_controller$plots_to_draw
  meta <- list(meta)
  pipeline_controller(meta,c(parameters$heatmap,parameters$cor_plot,parameters$sample_cluster),timer(plot_heatmap))
  pipeline_controller(meta,c(parameters$MA_plot,parameters$volcano_plot),timer(plot_pvalue_LFC))
  pipeline_controller(meta,c(parameters$loading_plot,parameters$pairsplot,parameters$scree_plot,parameters$biplot),timer(plot_PCA))
  pipeline_controller(meta,c(parameters$deg_heatmap,parameters$venn_diagram),timer(plot_venn_heatmap))
  pipeline_controller(meta,c(parameters$ridge_plot,parameters$dot_plot),timer(plot_functional_analysis))
}

data_transform_steps <- function(meta, plot_heatmap, plot_pvalue_LFC, plot_PCA, plot_venn_heatmap, plot_functional_analysis) {
  parameters <- meta$flow_controller$data_pipe_line
  meta <- list(meta)
  pipeline_controller(meta,parameters$dds_transformation,timer(data_dds_transformation))
  pipeline_controller(meta,parameters$normalized_data_generation,timer(data_normalized))
  pipeline_controller(meta,parameters$enrichment_analyses$Ctrl,timer(data_functional_analyses))
}

pipeline_controller <- function(FUN_input, decision_paras, FUN, fun_name = deparse(substitute(FUN))) {
  if (any(decision_paras)) {
    log_info('flow_controller set TRUE, Start {fun_name}.')
    do.call(FUN, FUN_input)
  } else {
    log_info('flow_controller set FALSE, Skip {fun_name}.')
  }
}

timer <- function(FUN, fun_name = deparse(substitute(FUN)), level = 0) {
  if (level == 0) {
    fun <- function(...) {
      timer_name <- paste(fun_name, "_level_", level, sep = "")
      start_time <- proc.time()
      result <- FUN(...)
      end_time <- proc.time()
      duration <- end_time - start_time
      cat('Function', timer_name, 'took', round(duration["elapsed"], 2), 'seconds to execute.\n')
      return(result)
    }
  } else {
    fun <- timer(FUN, fun_name, level = level - 1)
  }

  return(fun)
}

rds_to_csv <- function(meta){
  normalized_rds_to_csv <- function(meta){
    deseq2_dds <- meta$import_data('dds_vst.rds')
    write.csv(deseq2_dds,meta$get_file_path('data','dds_vst.csv'))
  }
 

  resLFC_rds_to_csv <- function(meta){
    resLFC_all <- meta$import_data('resLFC_all.rds')
    write.csv(resLFC_all,meta$get_file_path('data','resLFC_all.csv'))
    resLFC_all <- read.csv(meta$get_file_path('data','resLFC_all.csv'))
  }

  functional_all_rds_to_csv <- function(meta){
  functional_all <- meta$import_data('functional_all.rds')
  for (i in seq_along(functional_all)) {
    for (j in seq_along(functional_all[[i]])) {
      dat <- functional_all[[i]][[j]]
      write.csv(dat$results,meta$get_file_path('data',paste0('functional_',meta$group_list[i],'_',dat$analysis_name,'.csv')))
    }
  }
  }
  normalized_rds_to_csv(meta)
  resLFC_rds_to_csv(meta)
  functional_all_rds_to_csv(meta)
}
