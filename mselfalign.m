function IMs = mselfalign(IMs, AGE_FLAG)
%MSELFALIGN aligns the 3D image to its mirror reflection to make it
%symmetrical w.r.t. the middle plane of the 3D image
%   USAGE: IMsOUT = mselfalign({IM1P, IM2P, IM3P}, 'young')
%
%   Parameters
%   ==========
%   IMs             - list (of the 3D images, e.g. {IM1P, IM2P, ...})
%   AGE_FLAG        - string('young' | 'adult')
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

addpath(genpath('Scripts'));

if nargin < 2
    AGE_FLAG = 'young';
end

tic
for i = 1 : length(IMs)
    %Aligning the image to its mirror reflection
    TF = msimanneal(IMs{i}, flip(IMs{i}), AGE_FLAG, 'full');
    
    %Applying 1/2 transformation to make the image symmetrical
    TF.Vx = (TF.Vx + TF.Xs) / 2;
    TF.Vy = (TF.Vy + TF.Ys) / 2;
    TF.Vz = (TF.Vz + TF.Zs) / 2;
    
    %Rendering and saving the result
    IMs{i} = mtransform(IMs{i}, TF);
    figure, mshow(IMs{i}, flip(IMs{i}));
    toc; drawnow;
end

end
