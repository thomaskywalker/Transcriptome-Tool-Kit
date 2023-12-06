data_normalized <- function(meta) {
  dds <- meta$import_data('deseq2_dds.rds')
  # 4.0 Data trasformation for downstream analyses
  normalized_data <- assay(do.call(meta$data_parameters$NORMALIZED_FUNCTION,list(dds)))
  rds_file_path <- meta$get_file_path( 'data', paste0("dds_",meta$data_parameters$NORMALIZED_FUNCTION,".rds"))
  saveRDS(normalized_data, rds_file_path)
  log_info('{rds_file_path} saved')

  log_info('The input for downstream function contains :{resultsNames(dds)}')
  groups_list <- resultsNames(dds)[-1]
  groups_num <- length(groups_list)
  resLFC_all <- vector('list', groups_num)
  for (idx in 1:groups_num) {
    contra <- groups_list[[idx]]
    log_info('Preparing data of {contra},Shrink Log Fold Change by apeglm')
    resLFC <- lfcShrink(dds, coef = contra, type = 'apeglm')# you can try "normal" or "ashr" instead
    resOrdered <- resLFC[order(resLFC$padj),]
    resLFC_all[[idx]] <- resOrdered
  }

  rds_file_path <- meta$get_file_path( 'data', "resLFC_all.rds")
  saveRDS(resLFC_all, rds_file_path)
  log_info('{rds_file_path} saved')
}
