function res = lra_hdr(img)
%% pre-processing
gray = rgb2gray(img);
L = double(gray)/255;

%% Local Region Segementation

function R_map = local_region_segmentation(L, PBmin, PBmed, PBmax)
    [height, width] = size(L);
    map_size = height * width;
    H = imhist(L);
    PB = cumsum(H)/map_size;
    for i = 1:256
        if PB(i) <= PBmin
            Llow = i;
        elseif PB(i) <= PBmed
            Lmid = i;
        elseif PB(i) <= PBmax
            Lhigh = i;
        end
    end
    R_map = zeros(height,width);
    for i = 1:height
        for j = 1:width
            if L(i,j) <= Llow
                R_map(i,j) = 1;
            elseif L(i,j) <= Lmid
                R_map(i,j) = 2;
            elseif L(i,j) <= Lhigh
                R_map(i,j) = 3;
            else
                R_map(i,j) = 4;
            end
        end
    end

end

region_map = local_region_segmentation(gray,0.25,0.5,0.75);

%% Pseudo exposure generation

function Lw = pseudo_exposure(L,u,Ev,P,Lsmax)

    [height,width] = size(L);
    Lad = 1.0 + exp(u*Ev);
    Lw = zeros(height,width);
    for i = 1:height
        for j = 1:width
            if L(i,j) < 1.0
                Lw(i,j) = ((10^(-P)) * Lad * L(i,j)) / (1 - L(i,j));
            else
                Lw(i,j) = 10^(-P) * Lad * Lsmax;
            end
        end
    end
end

Lw1 = pseudo_exposure(L, 0.85, -1.0, 1.6, 382.5);
Lw2 = pseudo_exposure(L, 0.85, -0.5, 1.3, 382.5);
Lw3 = pseudo_exposure(L, 0.85, 0.0, 1.0, 382.5);
Lw4 = pseudo_exposure(L, 0.85, 0.5, 0.8, 382.5);
Lw5 = pseudo_exposure(L, 0.85, 1.0, 0.75, 382.5);
%imshow([Lw1,Lw2,Lw3;Lw4,Lw5,L]);

%% Calculate med for each luminance map
[height,width] = size(L);
c1 = 0;
c2 = 0;
c3 = 0;
c4 = 0;
Lmed1 = 0;
Lmed2 = 0;
%Lmed3 = 0.5 * (max(max(Lw3)) - min(min(Lw3)));
Lmed3 = mean(mean(Lw3));
Lmed4 = 0;
Lmed5 = 0;
region_img = zeros(height, width, 3);
for i = 1:height
    for j = 1:width
        if region_map(i,j) == 1
            Lmed5 = Lmed5 + Lw5(i,j);
            c1 = c1 + 1;
            region_img(i,j,1) = 255;
        elseif region_map(i,j) == 2
            Lmed4 = Lmed4 + Lw4(i,j);
            c2 = c2 + 1;
            region_img(i,j,2) = 255;
        elseif region_map(i,j) == 3
            Lmed2 = Lmed2 + Lw2(i,j);
            c3 = c3 + 1;
            region_img(i,j,3) = 255;
        elseif region_map(i,j) == 4
            Lmed1 = Lmed1 + Lw1(i,j);
            c4 = c4 + 1;
            region_img(i,j,1) = 255;
            region_img(i,j,2) = 255;
            region_img(i,j,3) = 255;
        end
    end
end
Lmed1 = Lmed1 / c4;
Lmed2 = Lmed2 / c3;
Lmed4 = Lmed4 / c2;
Lmed5 = Lmed5 / c1;

%% Weight function

function w = weight_function(L, alpha, Lmed)

    w = exp(-alpha * ((L - Lmed).^2 / (Lmed^2)));
end

w1 = weight_function(Lw1, 1, Lmed1);
w2 = weight_function(Lw2, 0.1, Lmed2);
w3 = weight_function(Lw3, 0.1, Lmed3);
w4 = weight_function(Lw4, 0.1, Lmed4);
w5 = weight_function(Lw5, 1, Lmed5);
%imshow([w1,w2,w3;w4,w5,L])

%% Fusion
Lw_comb = (w1 .* Lw1 + w2 .* Lw2 + w3 .* Lw3 + w4 .* Lw4 + w5 .* Lw5) ...
    ./ (w1+w2+w3+w4+w5);

%% Tone mapping

function img_out = local_tone_mapping(L_in,L_out,img_in,gamma)

    L_out = log(1 + L_out);
    img_out = img_in;
    for i = 1:3
        img_out(:,:,i) = (double(img_in(:,:,i)) ./ L_in) .^ gamma .* L_out;
    end
end

function img_out = color_reconstruction(img_in, L_out, L_in)

    [height,width] = size(L_in);
    img_out = zeros(height, width, 3);
    for i = 1:3
        img_out(:,:,i) = 0.5 * ((L_out ./ L_in) .* (double(img_in(:,:,i)) + ...
            255 * L_in) + double(img_in(:,:,i)) - 255 * L_in);
    end

    img_out = uint8(img_out);
end

img_tone = local_tone_mapping(L, Lw_comb, img, 1.3);
img_tone2 = color_reconstruction(img, L, Lw_comb);
%imshow([L Lw1 Lw2 Lw3 Lw4 Lw5 Lw_comb]);
res = imbilatfilt(img_tone);
end

