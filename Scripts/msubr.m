function IMOUT = msubr(IM, SCALE)
%MSUBR decreases resolution of the image
%
%   Parameters
%   ==========
%   IM        - array (input 3D image)
%   SCALE     - double (decrease in resolution, e.g. 4-fold)
%   IMOUT     - array (low-resolution 3D image)
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

IM = imresize(IM, 1 / SCALE, 'bilinear');
IM = permute(IM, [3 1 2]);
IM = imresize(IM, [round(size(IM, 1) / SCALE), size(IM, 2)], 'bilinear');
IMOUT = permute(IM, [2 3 1]);

IMOUT = (IMOUT - min(IMOUT(:))) / (max(IMOUT(:)) - min(IMOUT(:)));
IMOUT = IMOUT / sum(IMOUT(:));

end
