origPath = pwd;
cd ~/Softwares/unlocbox/;
init_unlocbox();
cd ~/Softwares/gspbox/;;
gsp_start;addpath(".");
cd(origPath);
x = csvread('x.csv');
[x_impute,network] = G2S3(x);
csvwrite('x_impute.csv', x_impute);