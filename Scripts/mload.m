function IM = mload(FNAMEFMT, DECREASE, RATIO)
%MLOAD Loads a 3D image from a series of 2D images with filenames matching
%the format FNAMEFMT. The 3D image has uniform resolution equal to 1/DECREASE
%of the original xy resolution of 2D images
%
%   Parameters
%   ==========
%   FNAMEFMT        - string (file name format like 's_C001Z*.tif')
%   DECREASE        - double (decrease in resolution, e.g. 4-fold)
%   RATIO           - double (ratio between xy and z resolution, e.g. 1/4)
%   IM              - array (loaded 3D image)
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

FILES_LIST = dir(FNAMEFMT);
NUM = length(FILES_LIST);
STEP = DECREASE / RATIO;

%Allocate an array for the 3D image
TMP = imread(FILES_LIST(1).name);
TMP = TMP(:, :, 1);
if DECREASE ~= 1
    TMP = imresize(TMP, 1 / DECREASE);
end
[M, N] = size(TMP);
IM = zeros(M, N, ceil(NUM / DECREASE));

%Load the images
for i = 1 : STEP : NUM
    fprintf('%d ', i)
    
    if i == round(i)
        TMP = imread(FILES_LIST(round(i)).name);
        TMP = TMP(:, :, 1);
    else %Interpolate
        TMP = imread(FILES_LIST(floor(i)).name) .* (ceil(i) - i);
        TMP = TMP + imread(FILES_LIST(ceil(i)).name) .* (i - floor(i));
        TMP = TMP(:, :, 1);
    end
    if DECREASE ~= 1 %Decrease resolution
        TMP = imresize(TMP, 1 / DECREASE);
    end
    IM(:, :, 1 + round(i / STEP)) = double(TMP);
end

fprintf('\n')

imagesc(max(IM, [], 3)); axis image; colormap hot; colorbar;

end
