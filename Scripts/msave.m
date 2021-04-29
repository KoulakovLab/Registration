function msave(IM, CHANNEL)
%MSAVE Saves a 3D image as a series of 2D images with filenames defined by
%the format FNAMEFMT.
%
%   Parameters
%   ==========
%   IM              - array (loaded 3D image)
%   CHANNEL         - number (of a channel)
%   FNAMEFMT        - string (file name format like 's_C001Z*.tif')
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

if nargin < 2
    CHANNEL = 1;
end

FNAMEFMT='s_C00%dZ%03d.tif';

IM = IM / max(IM(:));
for k = 1 : size(IM, 3)
    imwrite(IM(:, :, k), sprintf(FNAMEFMT, CHANNEL, k));
end

end
