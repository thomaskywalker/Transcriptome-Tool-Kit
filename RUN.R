start_analyse <- function(project_name){
  source('dependency.R')
  log_threshold('DEBUG')
  config <- yaml::read_yaml(file.path(project_name,"config.yaml"))
  std <- yaml::read_yaml('std.yaml')
  meta <- meta_data_factory(config,std)
  flow_control(meta)
}

# set working directory to location of RUN.R
# setwd('R')
# input project name
# start_analyse(project_name = 'setA')
# start_analyse(project_name = 'setB')
# start_analyse(project_name = 'setall')
# project_name = 'setall_AF'
# start_analyse(project_name = 'setall_AF')

start_analyse(project_name = '../gene_seq_analyses - 複製/setC')


