function IM = mlog(IM, SIZE, SIGMA)
%MLOG Laplacian of Gaussian filtering of 3D image in the frequency domain
%
%   Parameters
%   ==========
%   IM      - double (input 3D image to be filtered)
%   SIGMA   - double (standard deviation of the Gaussian filter)
%   SIZE    - double (filter size)
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

%Define the LoG filter
SIZE = floor((SIZE + 1) / 2 * 2);
FILTER = zeros(SIZE, SIZE, SIZE);

for i = 1 : SIZE
    for j = 1 : SIZE
        for k = 1 : SIZE
            RADIUS_SQUARED = (i - SIZE / 2) ^ 2 + (j - SIZE / 2) ^ 2 + ...
                (k - SIZE / 2) ^ 2;
            FILTER(i, j, k) = (2 * RADIUS_SQUARED - 3) .* ...
                exp( - RADIUS_SQUARED / (2 * SIGMA ^ 2));
        end
    end
end

FILTER = FILTER - mean(FILTER(:));
FILTER = FILTER / std(FILTER(:));

%Apply the filter to the 3D image in frequency domain
IM = ifftn(fftn(IM) .* fftn(FILTER, size(IM)));

IM(round(1 : end - SIZE / 2), round(1 : end - SIZE / 2), ...
    round(1 : end - SIZE / 2)) = ...
    IM(round(SIZE / 2 + 1 : end), round(SIZE / 2 + 1 : end), ...
    round(SIZE / 2 + 1 : end));

IM = - IM .* (IM < 0);
IM = IM / max(IM(:));

end
