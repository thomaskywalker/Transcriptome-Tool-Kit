data_dds_transformation <- function(meta) {
  # import data
  cts <- meta$import_data(meta$path$COUNT_DATA_PATH, sep = ",", header = TRUE, row.names = 1)  # cts (raw count)
  log_debug('Count data read from  :{meta$path$COUNT_DATA_PATH} .Dimensions are nCol:{ncol(cts)},nRow:{nrow(cts)}.')

  coldata <- meta$import_data(meta$path$COLUMN_DATA_PATH) #coldata (metadata)
  log_debug('Column data read from  :{meta$path$COLUMN_DATA_PATH} .Dimensions are nCol:{ncol(coldata)}, nRow:{nrow(coldata)}.')

  # reorder cts's columns based on row order of metadata (coldata)
  if (!(all(rownames(coldata) == colnames(cts)))){
    cts <- cts[, rownames(coldata)]
  }

  dds <- DESeqDataSetFromMatrix(
      countData = cts,
      colData = coldata,
      design = meta$data_parameters$EXP_DESIGN_FORMULA
  )
  log_info('Secessfully created DESeq Data by countData:{meta$path$COUNT_DATA_PATH},column data:{meta$path$COLUMN_DATA_PATH} ,design:{meta$data_parameters$EXP_DESIGN_FORMULA }.')


  # 1.1 Pre-filtering
  # perform a minimal pre-filtering to keep only rows that have at least 10 reads total
  dds <- (rowSums(counts(dds)) >= 10) %>% dds[., ]
  log_info('Perform a minimal pre-filtering to keep only rows that have at least 10 reads. {nrow(cts)-nrow(dds)} rows has been filtered.')

  # 1.2 Note on factor levels
  factors <- meta$data_parameters$REF_GROUP %>% attributes(.) %>% .[['names']]
  for ( name in factors){
    dds[[name]]  <- dds[[name]]  %>% as.factor() %>% relevel(.,ref=meta$data_parameters$REF_GROUP[[name]])
  }


  # 2.0 run DESeq function
  log_debug('Start estimating dispersion  by DESeq2, it may takes some times.')
  ## probably the most important command :)
  dds <- DESeq(dds,parallel =FALSE)

  fname <- meta$get_file_path('data','deseq2_dds.rds')
  saveRDS(dds, file = fname)
  log_info('Data has been saved as {fname}')

  return(dds)
  }

