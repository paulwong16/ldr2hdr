function res = CLA(img)

gray = rgb2gray(img);

%% Generating image sets
    function img_sets = img_set_generator(img)
        N = 10;
        M = 5;
        [height, width] = size(img);
        Y_k = zeros(height, width, 2*N+1);
        for k = -N:N
            Y_k(:,:,k + N + 1) = double(uint8(double(img) * (sqrt(2)^k)));
        end
        mini = inf;
        anc = 0;
        
        for i = 1:2*N+1
            if abs(mean(mean(Y_k(:,:,i))) - 128) < mini
                mini = abs(mean(mean(Y_k(:,:,i))) - 128);
                anc = i;
            end
        end
        img_sets = zeros(height, width, 2*M+1);
        for i = 1:2*M+1
            idx = anc - M - 1;
            img_sets(:,:,i) = Y_k(:,:,idx+i);
        end
    end

img_sets = img_set_generator(gray);

%% Threshold search

[thd0, thd1] = multithresh(gray,2);

%% JND-based Contrast Measure

    function C_t_k = jnd_contrast(img_sets)
        C_t_k = img_sets;
        [~, ~, set_nums] = size(img_sets);
        for k = 1:set_nums
            img_ = img_sets(:,:,k);
            tmp = JND_pixel(img_, 'Chou');
            C_t_k(:,:,k) = tmp;
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
            tmp = (img_set(:,:,k) - target) .^2;
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
        W_k = C_t_k .* E_t_k;
        weight_sum = sum(sum(W_k));
        W_k = W_k ./ weight_sum;
    end

W_k = weight_generator(C_t_k, E_t_k);

imshow(uint8(E_t_k(:,:,1)));

%% Fuse luminance map

    function F = fuse_sets(img_set, weight_set)
        L = 3;
        [height, width, sets] = size(img_set);
        coeff = zeros(height/2, width/2, sets, 4);
        for k = 1:sets
            [coeff(:,:,k,1), coeff(:,:,k,2), coeff(:,:,k,3), coeff(:,:,k,4)] = dwt2(img_set(:,:,k),'haar');
        end
        F_ = coeff .* impyramid(weight_set, 'reduce');
        F_ = sum(F_, 3);
        F = idwt2(F_(:,:,1), F_(:,:,2), F_(:,:,3), F_(:,:,4), 'haar');
    end

F = fuse_sets(img_sets, W_k);

res = E_t_k(:,:,1);

end

