clear all;
img = imread('./test5-.png');
img = uint8(double(img) / 16383 * 255);
hdr_res = imread('./test5_res.jpg');
%% test computing time
tic;
fwi_res = FWI(img);
toc;
text_str = cell(4,1);
text_str{1} = ['FU2016'];
text_str{2} = ['Time: ' num2str(toc) ' seconds'];
text_str{3} = ['GMSD: ' num2str(GMSD(rgb2gray(hdr_res), rgb2gray(fwi_res)))];
text_str{4} = ['NIQE: ' num2str(niqe_metric(fwi_res))];
position = [20 20;20 120;20 220;20 320]; 
fwi_show = insertText(fwi_res, position, text_str, 'FontSize',50);

tic;
his_res = histogram_hdr(img);
toc;
text_str = cell(4,1);
text_str{1} = ['IM2011'];
text_str{2} = ['Time: ' num2str(toc) ' seconds'];
text_str{3} = ['GMSD: ' num2str(GMSD(rgb2gray(hdr_res), rgb2gray(his_res)))];
text_str{4} = ['NIQE: ' num2str(niqe_metric(his_res))];
position = [20 20;20 120;20 220;20 320]; 
his_show = insertText(his_res, position, text_str, 'FontSize',50);

tic;
slh_res = SLH(img);
toc;
text_str = cell(4,1);
text_str{1} = ['PARK2019'];
text_str{2} = ['Time: ' num2str(toc) ' seconds'];
text_str{3} = ['GMSD: ' num2str(GMSD(rgb2gray(hdr_res), rgb2gray(slh_res)))];
text_str{4} = ['NIQE: ' num2str(niqe_metric(slh_res))];
position = [20 20;20 120;20 220;20 320]; 
slh_show = insertText(slh_res, position, text_str, 'FontSize',50);

tic;
lra_res = lra_hdr(img);
toc;
text_str = cell(4,1);
text_str{1} = ['WANG2015'];
text_str{2} = ['Time: ' num2str(toc) ' seconds'];
text_str{3} = ['GMSD: ' num2str(GMSD(rgb2gray(hdr_res), rgb2gray(lra_res)))];
text_str{4} = ['NIQE: ' num2str(niqe_metric(lra_res))];
position = [20 20;20 120;20 220;20 320]; 
lra_show = insertText(lra_res, position, text_str, 'FontSize',50);


tic;
ying_res = Ying_2017_CAIP(img);
toc;
text_str = cell(4,1);
text_str{1} = ['YING2017'];
text_str{2} = ['Time: ' num2str(toc) ' seconds'];
text_str{3} = ['GMSD: ' num2str(GMSD(rgb2gray(hdr_res), rgb2gray(ying_res)))];
text_str{4} = ['NIQE: ' num2str(niqe_metric(ying_res))];
position = [20 20;20 120;20 220;20 320]; 
ying_show = insertText(ying_res, position, text_str, 'FontSize',50);


tic;
cla_res = CLA(img);
toc;
text_str = cell(4,1);
text_str{1} = ['LEE2012'];
text_str{2} = ['Time: ' num2str(toc) ' seconds'];
text_str{3} = ['GMSD: ' num2str(GMSD(rgb2gray(hdr_res), rgb2gray(cla_res)))];
text_str{4} = ['NIQE: ' num2str(niqe_metric(cla_res))];
position = [20 20;20 120;20 220;20 320]; 
cla_show = insertText(cla_res, position, text_str, 'FontSize',50);

%% show result
padding = zeros(size(img), 'uint8') + 255;
imshow([img, fwi_show, his_show, slh_show; lra_show, ying_show, cla_show, hdr_res ]);
%imshow([img, fwi_res, his_res; slh_res, lra_res, ying_res]);