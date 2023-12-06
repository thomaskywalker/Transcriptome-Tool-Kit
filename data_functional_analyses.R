# 8.0 enrichment analyses
# We use GSEA methods to perform emrichment analyses i.e. gseGO, gseKEGG and GSEA
# contruct geneList object i.e. gene names as index and decreasing foldchange as value
data_functional_analyses <- function(meta) {

  msigdata <- function(species,cat, subcat) {
      db <- msigdbr(species, category = cat, subcategory = subcat) %>%
          dplyr::select(gs_name, entrez_gene) %>%
          as.data.frame()
      return(db)
  }

  Normal_analyses <- function(meta,dds_result, fname,chosen_analyses) {
      # prepare geneList object i.e.
      # gene names as index and decreasing foldchange as value
      analyses_results <- list()
      # 8.1 Gene Ontology
      # perform gseGO with minimum size of 10 and maximum size of 300
      if(chosen_analyses$GO){
      meta$efficiency_parameters$seed %>% set.seed()
      log_debug('start GO: {fname}')
      names(dds_result$log2FoldChange) <- map_key_to_col(rownames(dds_result),key=meta$data_parameters$origin_seq_id,column= meta$data_parameters$GO_id,db = meta$DB)
      genelist <- sort(dds_result$log2FoldChange, decreasing = TRUE)
      gores <- clusterProfiler::gseGO(
          geneList = genelist,
          OrgDb = meta$DB,
          eps = 0,
          ont = "ALL",
          minGSSize = 50,
          maxGSSize = 500,
          pvalueCutoff = 1,
          nPermSimple = 10000,
          verbose = FALSE
      )
      gores <- DOSE::setReadable(gores, OrgDb = meta$DB, keyType = meta$data_parameters$GO_id)
      analyses_results[['GO']] <- list(analysis_name = 'GO',results = gores)
      }

      # 8.2 KEGG
      # perform kegg with "cmt" database
      if(chosen_analyses$KEGG){
      names(dds_result$log2FoldChange) <- map_key_to_col(rownames(dds_result),key=meta$data_parameters$origin_seq_id,column= meta$data_parameters$KEGG_id,db = meta$DB)
      genelist <- sort(dds_result$log2FoldChange, decreasing = TRUE)
        meta$efficiency_parameters$seed %>% set.seed()
      log_debug('start KEGG: {fname}')
      kkres <- clusterProfiler::gseKEGG(
          geneList = genelist,
          organism = meta$data_parameters$KEGG_organism,
          pvalueCutoff = 1,
          eps = 0
      )
      analyses_results[['KEGG']] <- list(analysis_name = 'KEGG',results = kkres)
      }

      # 8.3 KEGG module
      # A KEGG Module is a group of genes and proteins working together
      # in a specific biological function or metabolic subsystem
      if(chosen_analyses$MKEGG){
      names(dds_result$log2FoldChange) <- map_key_to_col(rownames(dds_result),key=meta$data_parameters$origin_seq_id,column= meta$data_parameters$MKEGG_id,db = meta$DB)
      genelist <- sort(dds_result$log2FoldChange, decreasing = TRUE)
      meta$efficiency_parameters$seed %>% set.seed()
      log_debug('start MKEGG: {fname}')
      mkkres <- clusterProfiler::gseMKEGG(
          gene = genelist,
          organism = meta$data_parameters$KEGG_organism,
          pvalueCutoff = 1
      )
      analyses_results[['MKEGG']] <- list(analysis_name = 'MKEGG',results = mkkres)
      }

      # 8.4 Disease enrichment analysis
      # 8.4.1 diseaes ontology
      if(chosen_analyses$DO){
      names(dds_result$log2FoldChange) <- map_key_to_col(rownames(dds_result),key=meta$data_parameters$origin_seq_id,column= meta$data_parameters$DO_id,db = meta$DB)
      genelist <- sort(dds_result$log2FoldChange, decreasing = TRUE)
        meta$efficiency_parameters$seed %>% set.seed()
      log_debug('start DO: {fname}')
      dores <- DOSE::gseDO(
          geneList = genelist,
          minGSSize = 10,
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          verbose = FALSE
      )
      dores <- DOSE::setReadable(dores, meta$DB, keyType = meta$data_parameters$DO_id)
      analyses_results[['DO']] <- list(analysis_name = 'DO',results = dores)
      }

      # 8.5.2 network of cancer genes
      if( chosen_analyses$NCG){
      names(dds_result$log2FoldChange) <- map_key_to_col(rownames(dds_result),key=meta$data_parameters$origin_seq_id,column= meta$data_parameters$NCG_id,db = meta$DB)
      genelist <- sort(dds_result$log2FoldChange, decreasing = TRUE)
        meta$efficiency_parameters$seed %>% set.seed()
      log_debug('start NCG: {fname}')
      ncgres <- DOSE::gseNCG(
          geneList = genelist,
          pvalueCutoff =1, #0.5,
          pAdjustMethod = "BH",
          verbose = FALSE
      )
      ncgres <- DOSE::setReadable(ncgres, meta$DB, keyType = meta$data_parameters$NCG_id)
      analyses_results[['NCG']] <- list(analysis_name = 'NCG',results = ncgres)
      }

      # 8.5 reactome analysis
      # reactomePA database is used
      if(chosen_analyses$reactomePA){
      names(dds_result$log2FoldChange) <- map_key_to_col(rownames(dds_result),key=meta$data_parameters$origin_seq_id,column= meta$data_parameters$reactomePA_id,db = meta$DB)
      genelist <- sort(dds_result$log2FoldChange, decreasing = TRUE)
        meta$efficiency_parameters$seed %>% set.seed()
      log_debug('start reactomePA: {fname}')
      reactres <- ReactomePA::gsePathway(
          geneList = genelist,
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          verbose = FALSE
      )
      reactres <- DOSE::setReadable(reactres, meta$DB, keyType = meta$data_parameters$reactomePA_id)
      analyses_results[['reactomePA']] <- list(analysis_name = 'reactomePA',results = reactres)
      }
      # 8.6 GSEA
      # load molecular signature database
      if(chosen_analyses$GSEA){
        meta$efficiency_parameters$seed %>% set.seed()
      log_debug('start GSEA: {fname}')
        category_num <- length(meta$data_parameters$category)
        gsea <- vector('list',category_num)
        molecular_signature_db <- vector('list',category_num)
        for (i in 1:category_num) {
          molecular_signature_db[[i]] <- msigdata(meta$data_parameters$species,meta$data_parameters$category[[i]], NULL)
          analysis_name <- paste0('GSEA',meta$data_parameters$category[[i]])
          analyses_results[[analysis_name]] <- list(analysis_name = analysis_name,results = clusterProfiler::GSEA(geneList = genelist, TERM2GENE = molecular_signature_db[[i]], eps = 0))
        }
      }
      return(analyses_results)
  }

  ORA_analyses <- function(meta,dds_result, fname,chosen_analyses){
    symbol <- map_key_to_col(rownames(dds_result),key=meta$data_parameters$origin_seq_id,column= "SYMBOL",db=meta$DB)
    entrezid <- map_key_to_col(rownames(dds_result),key=meta$data_parameters$origin_seq_id,column= "ENTREZID",db=meta$DB)
    analyses_results <- list()
    ORA_id_symbol <- dds_result %>% subset(padj < 0.01, abs(log2FoldChange) >= 1) %>% rownames() %>% map_key_to_col(.,key=meta$data_parameters$origin_seq_id,column= "SYMBOL",db=meta$DB)
    ORA_id_entrezid <- dds_result %>% subset(padj < 0.01, abs(log2FoldChange) >= 1) %>% rownames() %>% map_key_to_col(.,key=meta$data_parameters$origin_seq_id,column= "ENTREZID",db=meta$DB)
      names(dds_result$log2FoldChange) <- entrezid
      genelist <- sort(dds_result$log2FoldChange, decreasing = TRUE)
    if(chosen_analyses$GO){
      meta$efficiency_parameters$seed %>% set.seed()
      log_debug('start GO-ORA: {fname}')
      gores <- clusterProfiler::enrichGO(
        gene = ORA_id_entrezid,
        universe = entrezid,
        OrgDb = meta$DB,
        ont = "ALL",
        minGSSize = 50,
        maxGSSize = 500,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        pAdjustMethod = "BH",

        readable = TRUE
      )
      gores <- DOSE::setReadable(gores, meta$DB, keyType =meta$data_parameters$GO_ORA_id)
      analyses_results[['GO_ORA']] <- list(analysis_name = 'GO_ORA',results = gores)
    }
    # 8.2 KEGG
    # perform kegg with "cmt" database
    if(chosen_analyses$KEGG){
      meta$efficiency_parameters$seed %>% set.seed()
      log_debug('start KEGG-ORA: {fname}')
      kkres <- clusterProfiler::enrichKEGG(
        gene = ORA_id_symbol,
        universe = symbol,
        organism = meta$data_parameters$KEGG_organism,
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
      analyses_results[['KEGG_ORA']] <- list(analysis_name = 'KEGG-ORA',results = kkres)
    }
    return(analyses_results)
  }


  analyses_flow <- function(meta) {
  chosen_analyses <- meta$flow_controller$data_pipe_line$enrichment_analyses

  resLFC_all <- meta$import_data('resLFC_all.rds')

    log_info("perform analyses")
    functional_all <- vector('list',length = length(meta$group_list))
    for(idx in 1:length(meta$group_list)){
        log_debug("{meta$group_list[idx]}")
        dun <- list()
      if(chosen_analyses$Normal$Ctrl){
        dun <- c(dun,Normal_analyses(meta,resLFC_all[[idx]], meta$group_list[idx],chosen_analyses=chosen_analyses$Normal))
      }
      if(chosen_analyses$ORA$Ctrl){
        dun <- c(dun,ORA_analyses(meta,resLFC_all[[idx]], meta$group_list[idx],chosen_analyses=chosen_analyses$ORA))
      }
      functional_all[[idx]] <- dun
    }
    return(functional_all)
  }

  functional_all <- analyses_flow(meta)
  saveRDS(functional_all,file.path(meta$path$OUTPUT_PATH,'data','functional_all.rds'))
  return(functional_all)
}
