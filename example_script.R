unlocbox.dir = '~/Softwares/unlocbox/'
gspbox.dir = '~/Softwares/gspbox/;'
g2s3.dir ='.'
source('g2s3_config.R')
g2s3_config(unlocbox.dir, gspbox.dir, g2s3.dir) # this will configure the run_G2S3.m file for imputation
dat = readRDS('example1.rds')
# run default 1 diffusion step
dat_impute = run_g2s3(dat) 

# run user defined number of steps (t=3)
dat_impute = run_g2s3(dat, t = 3) 
