function res = niqe_metric(img)

cd niqe_release
load modelparameters.mat
 
blocksizerow    = 96;
blocksizecol    = 96;
blockrowoverlap = 0;
blockcoloverlap = 0;

res = computequality(img,blocksizerow,blocksizecol,blockrowoverlap,blockcoloverlap, ...
    mu_prisparam,cov_prisparam);

cd ..
end

