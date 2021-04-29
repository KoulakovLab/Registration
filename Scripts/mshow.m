function mshow(IM1, IM2)
%MSHOW shows differences between a pair of images
%
%   Parameters
%   ==========
%   IM1, IM2    - arrays (containing 3D images to be displayed)
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

%Normalize the intensities of the images
IM1 = (IM1 - min(IM1(:))) / (max(IM1(:)) - min(IM1(:)));
IM2 = (IM2 - min(IM2(:))) / (max(IM2(:)) - min(IM2(:)));

%Crop the images to the same size
X_SIZE = min(size(IM1,1), size(IM2,1));
Y_SIZE = min(size(IM1,2), size(IM2,2));
Z_SIZE = min(size(IM1,3), size(IM2,3));

IM_DISPLAY = zeros(X_SIZE + Z_SIZE, Y_SIZE + Z_SIZE, 3);
IM1 = IM1(1 : X_SIZE, 1 : Y_SIZE, 1 : Z_SIZE);
IM2 = IM2(1 : X_SIZE, 1 : Y_SIZE, 1 : Z_SIZE);

%Plot the image projections
IM1_DISPLAY = [max(IM1, [], 3), squeeze(max(IM1, [], 2)); ...
    squeeze(max(IM1, [], 1))', zeros(Z_SIZE, Z_SIZE)];
IM2_DISPLAY = [max(IM2, [], 3), squeeze(max(IM2, [], 2)); ...
    squeeze(max(IM2, [], 1))', zeros(Z_SIZE, Z_SIZE)];

IM_DISPLAY(:,:,1) = max(IM1_DISPLAY, [], 3);
IM_DISPLAY(:,:,2) = max(IM2_DISPLAY, [], 3);

image(IM_DISPLAY);
axis image

end
