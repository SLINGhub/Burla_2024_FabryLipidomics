## Plasma and Platelet Lipidome Changes in Fabry Disease

*Bo Burla, Jeongah Oh, Albina Nowak, Nathalie Piraud, Eduardo Meyer, Ding Mei, Anne K Bendt, Jan-Dirk Studt, Beat M Frey, Federico Torta, Markus R Wenk, and Pierre-Alexandre Krayenbuehl*

Affiliations of main authors: *Singapore Lipidomics Incubator, Life Sciences Institute; Department of Biochemistry Yong Loo Lin School of Medicine; National University of Singapore*

2024

## Summary

This repository contains the R used to process and QC filter the raw data (peak areas), and for the data analysis and generating of figures and tables


## Getting the Code

All R code and datasets are provided as an [RStudio](https://www.rstudio.com/products/RStudio) project. The easiest way run this code is to clone the repository within RStudio (<https://github.com/SLINGhub/Burla_2024_FabryLipidomics.git>). Alternatively, you can also download the Github Repository and open the Rstudio project.

## Setting up the Code

By default, the dependencies for this project are managed using [renv](https://rstudio.github.io/renv/) to improve reproducibility and facilitate installation of required packages. After cloning, the required packages can be automatically installed by running following command in the console:

``` r
renv::restore()
```

Should you prefer to use your local R library instead, you can turn `renv` off by running following command:

``` r
renv::deactivate() 
```

## Directory Structure

-   `analysis` `Quarto` notebooks with data processing and analysis workflows
-   `data/`
    -   `LipidomicsAnalyses/` Peak area files (\*.CSV in subfolder MH) and analysis metadata ( \*.XLSM)
    -   `processed/` Processed data (concentrations, QC filtered), output of `01_data-preprocessing.qmd`
    -   `StudyMetadata/` Subject metadata and clinical laboratory analysis values
-   `output/` Generated figures and tables 
-   `R/` R scripts with functions for plotting and managing of lipid names
-   `renv/` Used by `renv`

## Running the Analysis

-   Open the folder `analysis`
-   To process the raw data run all chunks in `01_data-preprocessing.qmd`
-   To generate figures and tables run all chunks in `02_manuscript-figures-tables.qmd`

## Releases

v1.0.0: First submitted version

## Contact

-   Code and Data Analysis: *Bo Burla* ([lsibjb\@nus.edu.sg](mailto:lsibjb@nus.edu.sg){.email}) 
-   Manuscript: *Bo Burla* ([lsibjb\@nus.edu.sg](mailto:lsibjb@nus.edu.sg){.email}) and *Pierre-Alexandre Krayenbuhl* ([pierrea.krayenbuehl\@usz.ch](mailto:pierrea.krayenbuehl@usz.ch){.email}) 

## License

The code in this analysis is covered by the `MIT` license

## sessionInfo

``` r
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.0

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Asia/Singapore
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] patchwork_1.1.3  ggpubr_0.6.0     ggrepel_0.9.3    ggsignif_0.6.4   gt_0.9.0         rlang_1.1.1      glue_1.6.2       ggpattern_1.0.1  lubridate_1.9.3  forcats_1.0.0    stringr_1.5.0    dplyr_1.1.3     
[13] purrr_1.0.2      readr_2.1.4      tidyr_1.3.0      tibble_3.2.1     tidyverse_2.0.0  here_1.0.1       midar_0.1.1.9001 ggplot2_3.4.3   

loaded via a namespace (and not attached):
 [1] polynom_1.4-1         readxl_1.4.3          magrittr_2.0.3        clue_0.3-65           GetoptLong_1.0.5      matrixStats_1.0.0     compiler_4.3.1        mgcv_1.9-0            systemfonts_1.0.4    
[10] png_0.1-8             vctrs_0.6.3           quantreg_5.97         pkgconfig_2.0.3       shape_1.4.6           crayon_1.5.2          fastmap_1.1.1         backports_1.4.1       labeling_0.4.3       
[19] utf8_1.2.3            tzdb_0.4.0            ragg_1.2.5            MatrixModels_0.5-2    bit_4.0.5             xfun_0.40             broom_1.0.5           parallel_4.3.1        cluster_2.1.4        
[28] R6_2.5.1              stringi_1.7.12        RColorBrewer_1.1-3    car_3.1-2             cellranger_1.1.0      Rcpp_1.0.11           iterators_1.0.14      knitr_1.44            usethis_2.2.2        
[37] IRanges_2.34.1        Matrix_1.6-1.1        splines_4.3.1         timechange_0.2.0      tidyselect_1.2.0      abind_1.4-5           rstudioapi_0.15.0     yaml_2.3.7            doParallel_1.0.17    
[46] codetools_0.2-19      lattice_0.21-9        withr_2.5.1           survival_3.5-7        ggpp_0.5.4            zip_2.3.0             xml2_1.3.5            circlize_0.4.15       pillar_1.9.0         
[55] BiocManager_1.30.22   carData_3.0-5         renv_1.0.3            foreach_1.5.2         stats4_4.3.1          generics_0.1.3        vroom_1.6.4           rprojroot_2.0.3       S4Vectors_0.38.2     
[64] hms_1.1.3             munsell_0.5.0         scales_1.2.1          tools_4.3.1           ggpmisc_0.5.4-1       SparseM_1.81          openxlsx_4.2.5.2      fs_1.6.3              grid_4.3.1           
[73] colorspace_2.1-0      nlme_3.1-163          cli_3.6.1             textshaping_0.3.6     fansi_1.0.4           ComplexHeatmap_2.16.0 gtable_0.3.4          ggh4x_0.2.6           rstatix_0.7.2        
[82] digest_0.6.33         BiocGenerics_0.46.0   rjson_0.2.21          farver_2.1.1          htmltools_0.5.6       lifecycle_1.0.3       GlobalOptions_0.1.2   bit64_4.0.5           MASS_7.3-60          
```
