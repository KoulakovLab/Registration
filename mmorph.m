function mmorph(IMs, AGES_ADJ, SIZE_FLAG)
%MMORPH displays the image development dynamics
%   USAGE: mmorph({IM1T, IM2T, IM3T}, [1.3 0.7 2.4], 'halves')
%
%   Parameters
%   ==========
%   IMs             - list (of the 3D images, e.g. {IM1P, IM2P, ...})
%   AGES_ADJ        - array (of the adjusted ages of samples)
%   SIZE_FLAG       - string ('full' | 'halves')
%   FILE_PREFIX     - string (file name for the output video)
%   SCALE           - double (relative growth per day, e.g. 0.02 = 2%/day)
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

FILE_PREFIX = 'growth';
SCALE = 0.02;

if nargin < 3
    SIZE_FLAG = 'full';
end

%Prepare the video files
writerObjMean = VideoWriter(strcat(FILE_PREFIX, '_mean.avi'), ...
    'Uncompressed AVI');
open(writerObjMean);
writerObjDiff = VideoWriter(strcat(FILE_PREFIX, '_diff.avi'), ...
    'Uncompressed AVI');
open(writerObjDiff);

frame_mean.colormap = hot;
frame_diff.colormap = [];

%Prepare the 3D images
for i = 1 : length(IMs)
    IMs{i} = (IMs{i} - min(IMs{i}(:))) / (max(IMs{i}(:)) - min(IMs{i}(:)));
    IMs{i} = IMs{i} .^ 1.5;
end

%Write the video
[a0, b0, c0] = size(IMs{1});
if strcmp(SIZE_FLAG, 'halves')
    a0 = a0 * 2;
end
a = ceil(a0 * (1 + real(max(AGES_ADJ(:))) * SCALE));
b = ceil(b0 * (1 + real(max(AGES_ADJ(:))) * SCALE));
c = ceil(c0 * (1 + real(max(AGES_ADJ(:))) * SCALE));

scaledXY = zeros(a, b, 3);
scaledXZ = zeros(a, c, 3);
scaledYZ = zeros(c, b, 3);
combXYZ = zeros(a + c, b + c, 3);

IM_mean = zeros(size(IMs{1}));
IM_diff = zeros(size(IMs{1}));

for j = 1 : 0.01 : max(AGES_ADJ(:))
    IM_mean = IM_mean * 0;
    IM_diff = IM_diff * 0;
    for i = 1 : length(IMs)
        IM_mean = IM_mean + IMs{i} * (...
            (exp(-(j - AGES_ADJ(i)) ^ 2 * 2)) / ...
                sum(exp(-(j - AGES_ADJ) .^ 2 * 2)));
            
        IM_diff = IM_diff + IMs{i} * (...
            (exp(-(j - AGES_ADJ(i)) ^ 2 * 2)) / ...
                sum(exp(-(j - AGES_ADJ) .^ 2 * 2)) - ...
            (exp(-(j - 1 - AGES_ADJ(i)) ^ 2 * 2)) / ...
                sum(exp(-(j - 1 - AGES_ADJ) .^ 2 * 2)));
    end

    %Compute the image projections: 1) maximum intensity of the difference;
    %2) mean image; 3) minimum intensity of the difference
    projXY = zeros(a0, b0, 3);
    projXZ = zeros(a0, c0, 3);
    projYZ = zeros(c0, b0, 3);

    if strcmp(SIZE_FLAG, 'halves')
        projXY(:, :, 2) = [max(IM_mean, [], 3); flip(max(IM_mean, [], 3))];
        projXZ(:, :, 2) = [squeeze(max(IM_mean, [], 2)); ...
                            flip(squeeze(max(IM_mean, [], 2)))];
        projXY(:, :, 1) = [max(IM_diff, [], 3); flip(max(IM_diff, [], 3))];
        projXZ(:, :, 1) = [squeeze(max(IM_diff, [], 2)); ...
                            flip(squeeze(max(IM_diff, [], 2)))];
        projXY(:, :, 3) = [min(IM_diff, [], 3); flip(min(IM_diff, [], 3))];
        projXZ(:, :, 3) = [squeeze(min(IM_diff, [], 2)); ...
                            flip(squeeze(min(IM_diff, [], 2)))];
    else
        projXY(:, :, 2) = max(IM_mean, [], 3);
        projXZ(:, :, 2) = squeeze(max(IM_mean, [], 2));
        projXY(:, :, 1) = max(IM_diff, [], 3);
        projXZ(:, :, 1) = squeeze(max(IM_diff, [], 2));
        projXY(:, :, 3) = min(IM_diff, [], 3);
        projXZ(:, :, 3) = squeeze(min(IM_diff, [], 2));
    end
    projYZ(:, :, 2) = flip(squeeze(max(IM_mean, [], 1))');
    projYZ(:, :, 1) = flip(squeeze(max(IM_diff, [], 1))');
    projYZ(:, :, 3) = flip(squeeze(min(IM_diff, [], 1))');

    %Scale the projection images to reflect growth
    projXY = imresize(projXY, 1 + j * SCALE); 
    [A, B, ~] = size(projXY);
    scaledXY(1 + floor((a - A) / 2) : floor((a + A) / 2), ...
             floor(1 + (b - B) / 2) : floor((b + B) / 2), :) = projXY;
    projXZ = imresize(projXZ, 1 + j * SCALE); 
    [A, C, ~] = size(projXZ);
    scaledXZ(1 + floor((a - A) / 2) : floor((a + A) / 2), ...
             floor(1 + (c - C) / 2) : floor((c + C) / 2), :) = projXZ;
    projYZ = imresize(projYZ, 1 + j * SCALE); 
    [C, B, ~] = size(projYZ);
    scaledYZ(floor(1 + (c - C) / 2) : floor((c + C) / 2), ...
             1 + floor((b - B) / 2) : floor((b + B) / 2), :) = projYZ;

    %Combine the scaled projection images into the video frame
    combXYZ(:, :, 1) = [scaledXY(:, :, 1), scaledXZ(:, :, 1); ...
                    scaledYZ(:, :, 1), zeros(c)];
    combXYZ(:, :, 2) = abs([scaledXY(:, :, 2), scaledXZ(:, :, 2); ...
                    scaledYZ(:, :, 2), zeros(c)]);
    combXYZ(:, :, 3) = -[scaledXY(:, :, 3), scaledXZ(:, :, 3); ...
                    scaledYZ(:, :, 3), zeros(c)];

    %Apply colormaps to the video frames
    displ_diff = repmat(abs(combXYZ(:, :, 2)) .^ 0.5 / 30, 1, 1, 3);
    
    displ_diff(:, :, 1) = displ_diff(:, :, 1) + combXYZ(:, :, 3);
    displ_diff(:, :, 2) = displ_diff(:, :, 2) + combXYZ(:, :, 1) + ...
                          combXYZ(:, :, 3);
    displ_diff(:, :, 3) = displ_diff(:, :, 3) + combXYZ(:, :, 1);
    
    displ_diff = (min(1, max(0, (1 - 3 * displ_diff)))) .^ 2;
    displ_mean = (min(1, max(0, (combXYZ(:, :, 2) .^ 0.5))));

    %Write the video frames
    figure(1), imagesc(displ_mean); colormap hot; axis image; drawnow
    if (j == round(j))
        hgsave(strcat(FILE_PREFIX, '_mean_P', num2str(j, '%d'), '.fig'));
    end
    
    figure(2), image(displ_diff); axis image; drawnow
    if (j == round(j))
        hgsave(strcat(FILE_PREFIX, '_diff_P', num2str(j, '%d'), '.fig'));
    end
    
    frame_mean.cdata = uint8(63 * displ_mean);
    writeVideo(writerObjMean, frame_mean);
    frame_diff.cdata = displ_diff;
    writeVideo(writerObjDiff, frame_diff);
end

close(writerObjMean);
close(writerObjDiff);

end