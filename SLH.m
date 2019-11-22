function hdr = SLH(img)
%% Pre-processing

L = double(rgb2gray(img)) / 255;

%% Weighted Least Square Filter
function OUT = wlsFilter(IN, lambda, alpha, L)

    if(~exist('L', 'var')),
        L = log(IN+eps);
    end

    if(~exist('alpha', 'var')),
        alpha = 1.2;
    end

    if(~exist('lambda', 'var')),
        lambda = 1;
    end

    smallNum = 0.0001;

    [r,c] = size(IN);
    k = r*c;

    % Compute affinities between adjacent pixels based on gradients of L
    dy = diff(L, 1, 1);
    dy = -lambda./(abs(dy).^alpha + smallNum);
    dy = padarray(dy, [1 0], 'post');
    dy = dy(:);

    dx = diff(L, 1, 2); 
    dx = -lambda./(abs(dx).^alpha + smallNum);
    dx = padarray(dx, [0 1], 'post');
    dx = dx(:);


    % Construct a five-point spatially inhomogeneous Laplacian matrix
    B(:,1) = dx;
    B(:,2) = dy;
    d = [-r,-1];
    A = spdiags(B,d,k,k);

    e = dx;
    w = padarray(dx, r, 'pre'); w = w(1:end-r);
    s = dy;
    n = padarray(dy, 1, 'pre'); n = n(1:end-1);

    D = 1-(e+w+s+n);
    A = A + A' + spdiags(D, 0, k, k);

    % Solve
    OUT = A\IN(:);
    OUT = reshape(OUT, r, c);
    
end
I = wlsFilter(L);

%% Calculate Reflectance

R = log(L) - log(I);

%% Selective Reflectance Scaling

function res = srs(reflectance, illuminace)

    res = reflectance;
    r_R = 0.5;
    mean_I = mean(mean(illuminace));
    [height, width] = size(illuminace);
    for i = 1:height
        for j = 1:width
            if illuminace(i,j) > mean_I
                res(i,j) = reflectance(i,j) * (illuminace(i,j) / mean_I) ^ r_R;
            end
        end
    end

end

R_ = srs(R, L);

%% Virtual Illumination Generation

function I_k = vig(illuminace,inv_illuminace)

    inv_illuminace = inv_illuminace / (max(max(inv_illuminace)));
    mi = mean(mean(illuminace));
    maxi = max(max(illuminace));
    v1 = 0.2;
    v3 = mi;
    v2 = 0.5 * (v1 + v3);
    v5 = 0.8;
    v4 = 0.5 * (v3 + v5);
    
    function fv = scale_fun(v_, mean_i_, max_i_)

        r = 1.0 - mean_i_/max_i_;
        fv = r*(1/(1+exp(-1.0*(v_ - mean_i_))) - 0.5 );
    end

    fv_1 = scale_fun(v1, mi, maxi);
    fv_2 = scale_fun(v2, mi, maxi);
    fv_3 = scale_fun(v3, mi, maxi);
    fv_4 = scale_fun(v4, mi, maxi);
    fv_5 = scale_fun(v5, mi, maxi);
    I_1 = (1 + fv_1) * (illuminace + fv_1 * inv_illuminace);
    I_2 = (1 + fv_2) * (illuminace + fv_2 * inv_illuminace);
    I_3 = (1 + fv_3) * (illuminace + fv_3 * inv_illuminace);
    I_4 = (1 + fv_4) * (illuminace + fv_4 * inv_illuminace);
    I_5 = (1 + fv_5) * (illuminace + fv_5 * inv_illuminace);
    I_k = cat(3, I_1, I_2, I_3, I_4, I_5);
end

I_k = vig(L, 1.0 - L);

%% Tone Reproduction

function result_ = tonemap_slh(image, L, R_, I_k)

    L_k = exp(R_) .* I_k;
    W_k = I_k;
    for i = 1:5
        if i < 3
            W_k(:,:,i) = I_k(:,:,i) / max(max(I_k(:,:,i)));
        else
            temp = 0.5 * (1 - I_k(:,:,i));
            W_k(:,:,i) = temp / max(max(temp));
        end
    end
    L_ = sum(W_k .* L_k, 3) ./ sum(W_k,3);
    ratio = L_ ./ L;
    out = double(image) .* ratio;
    result_ = uint8(out);
end

hdr = tonemap_slh(img, L, R_, I_k);

end

