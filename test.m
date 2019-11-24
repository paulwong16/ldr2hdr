clear all;
img = imread('./test5-.png');
img = uint8(double(img) / 16383 * 255);
%% test computing time
tic
fwi_res = FWI(img);
toc

tic
his_res = histogram_hdr(img);
toc

tic
slh_res = SLH(img);
toc

tic
lra_res = lra_hdr(img);
toc


tic
ying_res = Ying_2017_CAIP(img);
toc


tic
cla_res = CLA(img);
toc

%% show result
hdr_res = imread('./test5_res.jpg');
imshow([img, fwi_res, his_res, slh_res; lra_res, ying_res, cla_res, hdr_res]);
%imshow([img, fwi_res, his_res; slh_res, lra_res, ying_res]);