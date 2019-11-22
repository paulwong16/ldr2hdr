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
        sigma_L=32;
        sigma_M=64;
        sigma_H=32;
        
    end

E_t_k = exposedness(gray, img_sets, thd0, thd1);

%% Weight generation

    function W_k = weight_generator(C_t_k, E_t_k)
        W_k = C_t_k .* E_t_k;
        weight_sum = sum(sum(W_k));
        W_k = W_k ./ weight_sum;
    end

W_k = weight_generator(C_t_k, E_t_k);
        
res = C_t_k;

end

