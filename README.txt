# User guide
- Add a folder for analyse job under this folder.
- create a folder named `data` under analyse job folder(in small case!)
- add read counts (.csv) and column data(.txt,which record the meta data) to `data` folder
- copy `std.yaml` to analyse job folder
- rename it to `config.yaml` 
- open `config.yaml` to adjust parameters.
- data_parameters and path must to be filled in, other parameters can change when needed.
- open `RUN.R`
- setwd to the folder `RUN.R` located (depends)
-  find the command `start_analyse(project_name = 'polyps_results_test_2')` ,and change 'polyps_results_test_2' to your analyse job folder name. 
- If you have multiple analyse job folder, you can copy it multiple times to different lines.
```
start_analyse(project_name = 'job1')
start_analyse(project_name = 'job2')
```
- copy the command `source(RUN.R)` to the console panel and press enter to start analyse. 


# Future develop route
- add function to analyse interaction term
- add upset plot to replace venn diagram for contrasts which more than five groups.    
- targets package integrate 
- a GUI interface (shiny app?)
- more plots for enrichment analyses
- PPI ?
- reduce redundant compute (especially functional analysis)
- renv or docker technique to reduce deploy difficulty (fix R version,package version)
- use .RData to save functions
- .Rmd-like report generation 
- log file
- complie function to speedup function 
- pararellized function 
- integrate fastp

# constraint
- only accept transcriptome count data
- can only accept 1 factor and no other effect at a time.
- need to fix one treatment for control group
- need Biological Replicates
- the category of treatments under 5 for vennplot

# sessionInfo
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22621)