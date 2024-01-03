meta_data_factory <- function(config,std)
{
  read_data <- function(file_path,...){
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
          readRDS(file_path,...)
          return("rds")
        }, error = function(e) {
          return("unknown")
        })
      }
    }
    type <- identify_file_type(file_path)
    if(type == 'csv')
      return(read.csv(file_path,...))
    if(type == 'txt')
      return(read.table(file_path,...))
    if(type == 'rds')
      return(readRDS(file_path,...))


  }
  get_file_path <- function(type,file_name){
    return(file.path(config$path$OUTPUT_PATH,type,file_name))
  }
  get_group_list <- function(){
    return(import_data(config$path$COLUMN_DATA_PATH)$group %>% factor(.) %>% relevel(. ,ref=config$data_parameters$REF_GROUP[[names(config$data_parameters$REF_GROUP)]]) %>% levels(.) %>% .[-1])
  }
  import_data <- function(data_name,...){
    data <- read_data(get_file_path('data',data_name),...)
    log_debug("import {data_name}")
    return(data)
  }
  calc_combinations <- function() {
    formula <- config$data_parameters$EXP_DESIGN_FORMULA
    factors <- config$data_parameters$REF_GROUP %>% attributes(.) %>% .[['names']]

    data <- import_data(config$path$COLUMN_DATA_PATH)
    for ( name in factors){
    data[[name]]  <- data[[name]]  %>% as.factor() %>% relevel(.,ref=config$data_parameters$REF_GROUP[[name]])
    }
    combinations <- model.matrix(as.formula(formula),data=data) %>% colnames()
    return(combinations)
  }


  retrieve_AnnotationHub <- function(DB_id){
    hub <- AnnotationHub::AnnotationHub()
    AH  <- hub[[DB_id]]
    return(AH)
  }
  load_orgdb <- function(db_name)  {
    # 尝试从本地路径加载数据库
    if (!is.null(db_name) && file.exists(db_name)) {
      message("Loading orgDb from local path: ", db_name)
      db <- loadDb(db_name)
      return(db)
    }

    # 尝试作为已安装的包来加载数据库
    if (requireNamespace(db_name, quietly = TRUE)) {
      library(db_name,character.only = TRUE)
      message("Loading orgDb from installed packages: ", db_name)
      db <- get(db_name)
      return(db)
    }
  return(retrieve_AnnotationHub(db_name))
}

  stopifnot('Fisrt class attributes were not set' = (names(config) == names(std)))
  stopifnot('data_parameters were missed. (compare to std.yaml)'=names(config$data_parameters) == names(std$data_parameters))
  stopifnot('path were missed. (compare to std.yaml)'=names(config$path) == names(std$path))
  stopifnot('OUTPUT_PATH do not exsits' = file.exists(config$path$OUTPUT_PATH))
  stopifnot('COUNT_DATA_PATH do not exsits' = file.exists(get_file_path('data',config$path$COUNT_DATA_PATH) ))
  stopifnot('COLUMN_DATA_PATH do not exsits' = file.exists(get_file_path('data',config$path$COLUMN_DATA_PATH) ))
  stopifnot('Column data file gene id do not contain all sample id in count data'= all(rownames(import_data(config$path$COLUMN_DATA_PATH)) %in% colnames(import_data(config$path$COUNT_DATA_PATH)) ))
  config$DB <- load_orgdb(config$data_parameters$DB_id)
  config$get_file_path <- get_file_path
  config$read_data <- read_data
  config$import_data <- import_data
  config$data_parameters$EXP_DESIGN_FORMULA <- as.formula(config$data_parameters$EXP_DESIGN_FORMULA)
  config$group_list <- get_group_list()
  config$data_parameters$combinations <- calc_combinations()
  config$data_parameters$combinations_num <-  length(config$data_parameters$combinations)
  file.copy(config$get_file_path('data',config$path$COUNT_DATA_PATH),config$get_file_path('data','count_data.csv'),overwrite=FALSE,copy.date = TRUE)
  file.copy(config$get_file_path('data',config$path$COLUMN_DATA_PATH),config$get_file_path('data','col_data.txt'),overwrite=FALSE,copy.date = TRUE)
  return(config)
  }

