# G2S3: a Sparse and Smooth Signal of Gene Graph-based imputation method for scRNA-seq data

G2S3 is an imputation method that applies graph signal processing to extract gene graph structure from scRNA-seq data and recover true expression levels by borrowing information from adjacent genes in the gene graph. It is computationally efficient for imputation in large-scale scRNA-seq datasets.


## Installation Instructions:

### Required software and packages
To use G2S3 in R, you will need to install both the R and Matlab toolboxes.

1. Matlab R2019b (https://au.mathworks.com)
  
2. R (http://www.r-project.org/)
  
3. UNLocBoX (https://epfl-lts2.github.io/unlocbox-html/)

4. GSPBOX (https://epfl-lts2.github.io/gspbox-html/download.html)

Install the required Matlab Toolboxes before running G2S3.m. 


## Usage instructions:
`g2s3_config.R`: R source code for configure running scripts for G2S3. It includes function to generate `run_G2S3.m` file in the folder which can be called in R.

`example_script.R`: an example to use G2S3 in R enviroment.

`run_G2S3.m`: example Matlab script to run `G2S3.m`.

`example1.rds`: example dataset.
