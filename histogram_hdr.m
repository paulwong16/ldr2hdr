function HDR = histogram_hdr(img, w)
%% io
img_hsv = rgb2hsv(img);
img_v = img_hsv(:,:,3);
img_histogram = imhist(img_v);

%% w: weighted factor in Eq. (1)
if  ~exist( 'w', 'var' )
    w = 0.6;
end

%% WHS function
function [H0,H1] = WHS(H,w)

    M = sum(H);
    tau_map = zeros(1, length(H));

    for t = 1:length(H)
        tau_map(t) = abs(w - (1/M * sum(H(1:t))));
    end

    [~, tau] = min(tau_map);
    H0 = H;
    H1 = H;

    for p = 1:length(H)
        if p <= tau
            H0(p) = H(p);
            H1(p) = 0;
        else
            H1(p) = H(p);
            H0(p) = 0;
        end
    end

end

%% Stretches the histogram of the V channel
[H0, H1] = WHS(img_histogram, w);

function [min_,max_] = his_stat_helper(H)
    
    for i=1:length(H)
        if H(i) ~= 0
            min_ = i;
            break
        end
    end

    for j=length(H):-1:1
        if H(j) ~= 0
            max_ = j;
            break
        end
    end
    
end

[min_H0, max_H0] = his_stat_helper(H0);
[min_H1, max_H1] = his_stat_helper(H1);
img0_v = (img_v * 255 - min_H0) / (max_H0 - min_H0);
img1_v = (img_v * 255 - min_H1) / (max_H1 - min_H1);
img0_ = img_hsv;
img1_ = img_hsv;
img0_(:,:,3) = img0_v;
img1_(:,:,3) = img1_v;
img = double(img);
img0_ = double(uint8(255 * hsv2rgb(img0_)));
img1_ = double(uint8(255 * hsv2rgb(img1_)));

%% LDR image fusion
function [weight_O, weight_N, weight_U] = gaussian_shaped(O, N, U, p)

    weight_O = 10/2 * exp(((255 - O + mean(mean(O))).^2) / (127.5^2) * 0.5);
    weight_N = 10/5 * exp(((255 - N + mean(mean(N))).^2) / (127.5^2) * 0.5);
    weight_U = 10/8 * exp(((255 - U + mean(mean(U))).^2) / (127.5^2) * 0.5);
end

[w_o,w_n,w_u] = gaussian_shaped(img0_v, img_v, img1_v, 1);
w_sum = w_o + w_n + w_u;
w_o = w_o ./ w_sum;
w_n = w_n ./ w_sum;
w_u = w_u ./ w_sum;
HDR = (w_o .* img0_ + w_n .* img + w_u .* img1_) ./ (w_o + w_n + w_u);
HDR = uint8(HDR);
end

