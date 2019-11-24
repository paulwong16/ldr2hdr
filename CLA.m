function res = CLA(img)

gray = rgb2gray(img);

%% Generating image sets
    function img_sets = img_set_generator(img)
        N = 10;
        M = 5;
        [height, width] = size(img);
        Y_k = zeros(height, width, 2*N+1, 'uint8');
        for k = -N:N
            Y_k(:, :, k+N+1) = double(img)*(sqrt(2)^(k));
        end
        mini = inf;
        anc = 0;
        
        for i = 1:2*N+1
            q = abs(mean(mean(Y_k(:,:, i)))-128);
            if q < mini
                mini = q;
                anc = i;
            end
        end
        img_sets = zeros(height, width, 2*M+1, 'uint8');
        for i = 1:2*M+1
            idx = anc - M - 1;
            img_sets(:,:,i) = Y_k(:,:,idx+i);
        end
    end

img_sets = img_set_generator(gray);


%% Threshold search
thresholds = multithresh(gray,2);
thd0 = thresholds(1);
thd1 = thresholds(2);
%% JND-based Contrast Measure

    function C_t_k = jnd_contrast(img_sets)
        [height, width, set_nums] = size(img_sets);
        C_t_k = zeros(height, width, set_nums, 'double');
        for k = 1:set_nums
            img_ = img_sets(:,:,k);
            C_t_k(:,:,k) = JND_CON(img_);
        end
    end

C_t_k = jnd_contrast(img_sets);

%% Well-exposedness measure

    function E_t_k = exposedness(img, img_set, thd0, thd1)
        img = double(img);
        sigma_L=32;
        sigma_M=64;
        sigma_H=32;
        [height,width] = size(img);
        mu_l = 0;
        mu_l_num = 0;
        mu_h = 0;
        mu_h_num = 0;
        region_map = zeros(height, width);
        for i = 1:height
            for j = 1:width
                if img(i,j) < thd0
                    region_map(i,j) = 1;
                    mu_l = mu_l + img(i,j);
                    mu_l_num = mu_l_num + 1;
                elseif img(i,j) > thd1
                    region_map(i,j) = 3;
                    mu_h = mu_h + img(i,j);
                    mu_h_num = mu_h_num + 1;
                else
                    region_map(i,j) = 2;
                end
            end
        end
        
        mu_l = mu_l / mu_l_num;
        mu_h = mu_h / mu_h_num;
        r_l = mu_l_num / (height*width);
        target = img;
        
        for i = 1:height
            for j = 1:width
                switch region_map(i,j)
                    case 1
                        if mu_l > 64
                            target(i,j) = 64;
                        elseif mu_l >= 32 && mu_l <= 64
                            target(i,j) = mu_l;
                        elseif mu_l < 32 && r_l > 0.5
                            target(i,j) = 64;
                        elseif mu_l < 32 && r_l >= 0.25 && r_l <= 0.5
                            target(i,j) = r_l * 128;
                        elseif mu_l < 32 && r_l < 0.25
                            target(i,j) = 32;
                        end
                    case 2
                        target(i,j) = 128;
                    case 3
                        if mu_h > 224
                            target(i,j) = mu_h;
                        elseif mu_h >= 192 && mu_h <= 224
                            target(i,j) = 224;
                        else
                            target(i,j) = 192;
                        end
                end      
            end
        end
        
        [height,width,img_set_size] = size(img_set);
        E_t_k = zeros(height,width,img_set_size);
        for k = 1:img_set_size
            tmp = (double(img_set(:,:,k)) - target) .^2;
            tmp_l = exp(-(tmp / (2 * sigma_L^2)));
            tmp_m = exp(-(tmp / (2 * sigma_M^2)));
            tmp_h = exp(-(tmp / (2 * sigma_H^2)));
            for i = 1:height
                for j = 1:width
                    switch region_map(i,j)
                        case 1
                            E_t_k(i,j,k) = tmp_l(i,j);
                        case 2
                            E_t_k(i,j,k) = tmp_m(i,j);
                        case 3
                            E_t_k(i,j,k) = tmp_h(i,j);
                    end
                end
            end
        end
        
    end

E_t_k = exposedness(gray, img_sets, thd0, thd1);

%% Weight generation

    function W_k = weight_generator(C_t_k, E_t_k)
        W_k = double(C_t_k) .* E_t_k;
        weight_sum = sum(W_k, 3);
        W_k = W_k ./ weight_sum;
    end

W_k = weight_generator(C_t_k, E_t_k);



%% Fuse luminance map

    function F = fuse_sets(img_set, weight_set)
        
        W_1 = impyramid(weight_set, 'reduce');
        W_2 = impyramid(W_1, 'reduce');
        W_3 = impyramid(W_2, 'reduce');
        
        coeff_1 = zeros([size(W_1), 4]);
        coeff_2 = zeros([size(W_2), 4]);
        coeff_3 = zeros([size(W_3), 4]);
        
        sets = size(img_set, 3);
        
        for k = 1:sets
            [coeff_1(:,:,k,1), coeff_1(:,:,k,2), coeff_1(:,:,k,3), coeff_1(:,:,k,4)] = dwt2(img_set(:,:,k),'haar');
        end
        for k = 1:sets
            [coeff_2(:,:,k,1), coeff_2(:,:,k,2), coeff_2(:,:,k,3), coeff_2(:,:,k,4)] = dwt2(coeff_1(:,:,k,1),'haar');
        end
        for k = 1:sets
            [coeff_3(:,:,k,1), coeff_3(:,:,k,2), coeff_3(:,:,k,3), coeff_3(:,:,k,4)] = dwt2(coeff_2(:,:,k,1),'haar');
        end
        
        tmp3 = coeff_3 .* W_3;
        tmp3 = sum(tmp3, 3);
        tmp2_1 = idwt2(tmp3(:,:,1,1), tmp3(:,:,1,2), tmp3(:,:,1,3), tmp3(:,:,1,4), 'haar');
        tmp2_2_4 = coeff_2 .* W_2;
        tmp2 = sum(tmp2_2_4, 3);
        [height, width] = size(W_2(:,:,1));
        tmp2_1 = tmp2_1(1:height, 1:width);
        tmp1_1 = idwt2(tmp2_1, tmp2(:,:,1,2), tmp2(:,:,1,3), tmp2(:,:,1,4), 'haar');
        tmp1_2_4 = coeff_1 .* W_1;
        tmp1 = sum(tmp1_2_4, 3);
        [height, width] = size(W_1(:,:,1));
        tmp1_1 = tmp1_1(1:height, 1:width);
        F = idwt2(tmp1_1, tmp1(:,:,1,2), tmp1(:,:,1,3), tmp1(:,:,1,4), 'haar');
        
    end

F = fuse_sets(img_sets, W_k);
%F = sum(double(img_sets) .* W_k, 3);
%imshow(uint8(F));

%% Color reconstruction

    function img_out = color_reconstruction(img_in, L_out, L_in)
        
        L_in = double(L_in);
        [height,width] = size(L_in);
        img_out = zeros(height, width, 3);
        for i = 1:3
            img_out(:,:,i) = 0.5 * ((L_out ./ L_in) .* (double(img_in(:,:,i)) + ...
                L_in) + double(img_in(:,:,i)) - L_in);
        end

        img_out = uint8(img_out);
    end

[height_, width_] = size(gray);
F = F(1:height_, 1:width_);
res = color_reconstruction(img, F, gray);
res = imbilatfilt(res);

end

