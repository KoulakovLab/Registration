function IMs = mprepare(IMs)
%MPREPARE prepares the brain images for registration by 1) removing
%background; 2) equalizing intensity within and across samples; 3)
%pre-aligning samples along the image axis; 4) cropping the images
%
%   Parameters
%   ==========
%   IMs             - list (of the 3D images, e.g. {IM1P, IM2P, ...})
%   MIDTONE_THLD    - double (fraction of max intensity separating tissue
%                             autofluorescence and the signal)
%   SHAPE_THLD      - double (fraction of non-zero voxels in the binarized
%                             image separating brain tissue from background)
%   SHAPE_MARGIN    - double (fraction of the sample size to be added on
%                             each side ot the sample as the empty margins)
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

MIDTONE_THLD = 0.15;
SHAPE_THLD = 0.05;
SHAPE_MARGIN = 0.10;

close all

sizes = zeros(length(IMs), 3, 2);
threshold_colormap = hot;
threshold_colormap(1, :) = [0, 0, 1];

for i = 1 : length(IMs)
    %Threshold to remove background
    figure
    max_fluor = max(IMs{i}(:));
    im_display = [IMs{i}(:, :, round(end * 1 / 8)), ...
                  IMs{i}(:, :, round(end * 3 / 8)); ...
                  IMs{i}(:, :, round(end * 5 / 8)), ...
                  IMs{i}(:, :, round(end * 7 / 8))];
    width_im = size(im_display, 2);
    height_im = size(im_display, 1);
    im_selection = repmat(10 .^ ((1 : width_im) / width_im * 3 - 3) * ...
                                max_fluor, 100, 1);
    threshold = 0;
    ysel = 0;
    while true
        im_threshold = log10([im_selection; ...
            im_display .* (im_display > threshold)]);
        imagesc(im_threshold), axis image, colormap(threshold_colormap)
        title('Select threshold on colormap; click outside colormap to exit');
        hold on
        line([0 width_im], [100 100], 'color', 'w', 'linewidth', 3);
        line([0 width_im], (100 + height_im / 2) * [1 1], 'color', 'w');
        line(width_im * [0.5 0.5], [100 height_im + 100], 'color', 'w');
        line(ysel * [1 1], [0 100], 'color', 'w', 'linewidth', 5);
        hold off
        
        [ysel, xsel] = ginput(1);
        if(xsel < 1 || ysel < 1 || xsel > 100 || ysel > width_im)
            break
        else
            threshold = im_selection(1, round(ysel));
        end
    end
    IMs{i} = IMs{i} .* (IMs{i} > threshold);

    %Equalization of brightness in optical layers
    mean_backgr = mean(IMs{i}(find((IMs{i} < MIDTONE_THLD * max_fluor) .* ...
        (IMs{i} > 0))));

    for k = 1 : size(IMs{i}, 3)
        tmp = IMs{i}(:, :, k);
        curr_backgr = mean(tmp(find((tmp < MIDTONE_THLD * max_fluor) .* ...
            (tmp > 0))));
        IMs{i}(:, :, k) = tmp / curr_backgr * mean_backgr;
    end
    
    %Histogram equalization across samples
    tmp = histeq(IMs{i}(IMs{i} > 0) / max(IMs{i}(:)), ...
        imhist(IMs{1}(IMs{1} > 0) / max(IMs{1}(:)), 1000)) * max(IMs{1}(:));
    IMs{i}(IMs{i} > 0) = tmp;
    
    %PCA-based rotation and scaling
    ix = find(IMs{i} > 0);
    [x, y, z] = ind2sub(size(IMs{i}), ix);
    [a, ~, sz] = pca([y, x, z]);
    if i == 1
        sz_ref = sz;
    else
        a = a * diag((sz_ref(:) ./ sz(:)) .^ 0.5);
    end
    IMs{i} = imwarp(IMs{i}, affine3d([[a, [0; 0; 0]]; 0 0 0 1]));

    %Automatic crop
    im_bw = double(IMs{i} > 0);
    all_dims = 1 : 3;
    
    for dim = all_dims
        hidden_dims = all_dims(~ismember(all_dims, dim));
        im_bw_proj = squeeze(sum(sum(im_bw, hidden_dims(1)), hidden_dims(2)));
        im_bw_proj = (im_bw_proj > SHAPE_THLD * max(im_bw_proj(:)));
        sizes(i, dim, 1) = find(im_bw_proj, 1, 'first');
        sizes(i, dim, 2) = find(im_bw_proj, 1, 'last');
        margin = round((sizes(i, dim, 2) - sizes(i, dim, 1)) * SHAPE_MARGIN);
        sizes(i, dim, 1) = sizes(i, dim, 1) - margin;
        sizes(i, dim, 2) = sizes(i, dim, 2) + margin;
    end

    im_tmp = zeros(sizes(i, 1, 2) - sizes(i, 1, 1), ...
            sizes(i, 2, 2) - sizes(i, 2, 1), ...
            sizes(i, 3, 2) - sizes(i, 3, 1));
        
    im_cut = IMs{i}(max(1, sizes(i, 1, 1)) : min(end, sizes(i, 1, 2)), ...
            max(1, sizes(i, 2, 1)) : min(end, sizes(i, 2, 2)), ...
            max(1, sizes(i, 3, 1)) : min(end, sizes(i, 3, 2)));
        
    im_tmp(max(0, -sizes(i, 1, 1)) + (1 : size(im_cut, 1)), ...
           max(0, -sizes(i, 2, 1)) + (1 : size(im_cut, 2)), ...
           max(0, -sizes(i, 3, 1)) + (1 : size(im_cut, 3))) = im_cut;
    
    IMs{i} = im_tmp;

    %Manual flip
    figure
    [xsize, ysize, zsize] = size(IMs{i});
    while true
        mshow(IMs{max(1, i - 1)} .^ 0.5, IMs{i} .^ 0.5), axis image
        title('Click on any image to flip; click outside the box to exit')
        [ysel, xsel] = ginput(1);

        if(xsel < 1 || ysel < 1 || xsel > xsize + zsize || ...
                ysel > ysize + zsize)
            break
        elseif (xsel < xsize) && (ysel < ysize)
            IMs{i} = flip(IMs{i}, 1); IMs{i} = flip(IMs{i}, 2);
        elseif (xsel < xsize) && (ysel > ysize)
            IMs{i} = flip(IMs{i}, 1); IMs{i} = flip(IMs{i}, 3);
        elseif (xsel > xsize) && (ysel < ysize)
            IMs{i} = flip(IMs{i}, 2); IMs{i} = flip(IMs{i}, 3);
        end
    end
end

%Format all brains to a same size
max_size = max(sizes(:, :, 2) - sizes(:, :, 1), [], 1) + 1;
for i = 1 : length(IMs)
    im_tmp = zeros(max_size);
    curr_size = size(IMs{i});
    margin = max(0, floor((max_size - curr_size) / 2));
    im_tmp(margin(1) + (1 : curr_size(1)), ...
        margin(2) + (1 : curr_size(2)), ...
        margin(3) + (1 : curr_size(3))) = IMs{i};
    IMs{i} = im_tmp;
end

end
