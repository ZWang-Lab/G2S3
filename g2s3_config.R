g2s3_config <- function(unlocbox.dir, gspbox.dir, g2s3.dir = ".", t = 1){
    cmd.out <- paste0('origPath = pwd;\n')
    cmd.out <- paste0(cmd.out,"cd ", unlocbox.dir, ";\n")
    cmd.out  <- paste0(cmd.out, "init_unlocbox();\n")
    cmd.out <- paste0(cmd.out, "cd ", gspbox.dir, ";\n")
    cmd.out <- paste0(cmd.out,"gsp_start;")
    cmd.out <- paste0(cmd.out,"cd(origPath);\n")
    cmd.out <- paste0(cmd.out, "addpath(\"", g2s3.dir, "\");\n")
    cmd.out <- paste0(cmd.out,"cd(origPath);\n")
    cmd.out <- paste0(cmd.out, "x = csvread('x.csv');\n")
    cmd.out  <- paste0(cmd.out, "[x_impute,network] = G2S3(x, 't', ",t, ");\n")
    cmd.out <- paste0(cmd.out,"csvwrite('x_impute.csv', x_impute);")
    cat(cmd.out,file=paste0("run_G2S3_", t, ".m"),append=F)
}



run_g2s3 <- function(df, t = 1){
    # suppose df is the scRNA-seq data that is Gene * Cell
    x  = t(as.matrix(df))
    write.table(x, file='x.csv', sep=",", row.names=FALSE, col.names=FALSE)
    system(paste0("matlab -nodisplay -r \"run('run_G2S3_", t, ".m'); exit\""))
    x_g2s3 = read.table("x_impute.csv", sep = ",");x_g2s3 = t(x_g2s3)
    file.remove('x.csv');file.remove('x_impute.csv')
    dimnames(x_g2s3) = dimnames(df)
    return(x_g2s3)
}