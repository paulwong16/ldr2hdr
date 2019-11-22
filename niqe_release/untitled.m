load modelparameters.mat
 
blocksizerow    = 96;
blocksizecol    = 96;
blockrowoverlap = 0;
blockcoloverlap = 0;

niqe_res = zeros(20,5);

for k = 1:20
    img_name = sprintf('/Users/wangzhijie/Desktop/ECE613/research_project/experiments/FWI/%d-hdr.jpg', k);
    img = imread(img_name);
    niqe_res(k,1) = computequality(img,blocksizerow,blocksizecol,blockrowoverlap,blockcoloverlap, ...
    mu_prisparam,cov_prisparam);
end

for k = 1:20
    img_name = sprintf('/Users/wangzhijie/Desktop/ECE613/research_project/experiments/HIS/%d-hdr.jpg', k);
    img = imread(img_name);
    niqe_res(k,2) = computequality(img,blocksizerow,blocksizecol,blockrowoverlap,blockcoloverlap, ...
    mu_prisparam,cov_prisparam);
end

for k = 1:20
    img_name = sprintf('/Users/wangzhijie/Desktop/ECE613/research_project/experiments/LRA/%d-hdr.jpg', k);
    img = imread(img_name);
    niqe_res(k,3) = computequality(img,blocksizerow,blocksizecol,blockrowoverlap,blockcoloverlap, ...
    mu_prisparam,cov_prisparam);
end

for k = 1:20
    img_name = sprintf('/Users/wangzhijie/Desktop/ECE613/research_project/experiments/SLH/%d-hdr.jpg', k);
    img = imread(img_name);
    niqe_res(k,4) = computequality(img,blocksizerow,blocksizecol,blockrowoverlap,blockcoloverlap, ...
    mu_prisparam,cov_prisparam);
end

for k = 1:20
    img_name = sprintf('/Users/wangzhijie/Desktop/ECE613/research_project/experiments/YING/%d-hdr.jpg', k);
    img = imread(img_name);
    niqe_res(k,5) = computequality(img,blocksizerow,blocksizecol,blockrowoverlap,blockcoloverlap, ...
    mu_prisparam,cov_prisparam);
end

csvwrite('/Users/wangzhijie/Desktop/ECE613/research_project/experiments/niqe.csv', niqe_res);