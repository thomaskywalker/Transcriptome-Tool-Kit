start_analyse <- function(project_name){
  t1 <- Sys.time()
  source('dependency.R')
  # log_threshold('DEBUG')
  log_threshold('INFO')
  config <- yaml::read_yaml(file.path(project_name,"config.yaml"))
  std <- yaml::read_yaml('std.yaml')
  meta <- meta_data_factory(config,std,project_name)
  flow_control(meta)
  rds_to_csv(meta)
  t2 <- Sys.time()
  print(t2-t1)
}



start_analyse(project_name = 'example_data')
