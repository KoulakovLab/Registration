function [IMs, AGES_ADJ] = mwarp(IMs, AGES, AGE_FLAG, SIZE_FLAG)
%MWARP aligns the 3D images in sequence defined by their ages. The first
%image of each age group serves as a reference image. The reference images
%are aligned to each other in sequence (e.g. P0#1 <- P1#1 <- P2#1 ...),
%then the images within each age group are aligned to the transformed
%reference images (e.g. P0#2 <- P0#1, P0#3 <- P0#1, P1#2 <- P1#1 ...)
%   USAGE: [IMsOUT, AGES_ADJ] = mwarp(...
%               {IM1P, IM2P, IM3P}, [1 1 2], 'young', 'halves')
%
%   Parameters
%   ==========
%   IMs             - list (of the 3D images, e.g. {IM1P, IM2P, ...})
%   AGES            - array (of experimentally recorded ages of samples)
%   AGES_ADJ        - array (of the adjusted ages of samples)
%   AGE_FLAG        - string('young' | 'adult')
%   SIZE_FLAG       - string('full' | 'halves')
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

addpath(genpath('Scripts'));

if nargin < 3
    AGE_FLAG = 'young';
end
if nargin < 4
    SIZE_FLAG = 'full';
end

%Find reference images, i.e. the first images in each age group
ref = ones(1, length(AGES));
for j = 2 : length(AGES)
    if AGES(j) == AGES(j - 1)
        ref(j) = ref(j - 1);
    else
        ref(j) = j;
    end
end

%If working with 'halves', split images into separate images of hemispheres
if strcmp(SIZE_FLAG, 'halves')
    ref = [ref, zeros(size(ref))];
    AGES = [AGES, zeros(size(AGES))];
    
    for i = length(IMs) : -1 : 1 %inverse order to prevent overwriting
        IMs{2 * i} = flip(IMs{i}(ceil(end / 2) + 1 : end, :, :));
        IMs{2 * i - 1} = IMs{i}(1 : round(end / 2), :, :);
        ref(2 * i) = ref(i); ref(2 * i - 1) = ref(i);
        AGES(2 * i) = AGES(i); AGES(2 * i - 1) = AGES(i);
    end
end

%Align the images to reference images
tic
for i = (2 + strcmp(SIZE_FLAG, 'halves')) : length(IMs)
    
    if strcmp(SIZE_FLAG, 'full') || (mod(i, 2) == 1)
        TF = msimanneal(IMs{i}, IMs{ref(i - 1)}, AGE_FLAG, SIZE_FLAG);
    end
    
    %In case of 'halves' same TF applies to both hemispheres
    IMs{i} = mtransform(IMs{i}, TF);
    figure, mshow(IMs{i} .^ 0.5, IMs{ref(i - 1)} .^ 0.5); drawnow;
    toc; 
end

%Compute the adjusted ages based on similarities between the images
AGES_ADJ = mdate(IMs, AGES);

end
