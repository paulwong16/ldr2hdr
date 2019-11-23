function C = JND_CON(img)
%JND_CON Summary of this function goes here
%   Detailed explanation goes here
img = double(img);
[row, col] = size(img);
C = img;
for i = 2:row-1
    for j = 2:col-1
        tmp = img(i-1:i+1, j-1:j+1);
        maxi = max(max(tmp));
        mini = min(min(tmp));
        avg = mean(mean(tmp));
        diff = maxi - mini;
        JND = computeJND(avg);
        if diff < JND
            C(i, j) = 1/256;
        else
            C(i, j) = (diff + 1)/256;
        end
    end
    C(1, :) = C(2, :);
    C(:, 1) = C(:, 2);
    C(row, :) = C(row-1, :);
    C(:, col) = C(:, col-1);
end