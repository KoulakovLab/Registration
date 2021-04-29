function IMOUT = mtransform(IM, TF, MODE)
%MTRANSFORM transforms the 3D image
%
%   Parameters
%   ==========
%   IM          - array (input 3D image)
%   TF          - structure (transformation)
%   MODE        - string ('relative' | 'absolute')
%   IMOUT       - array (transformed 3D image)
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

[M, N, K] = size(IM);
IMOUT = zeros(M, N, K);

%Scaling the transformation to the size of the 3D image
if (nargin < 3) || strcmp(MODE, 'relative')
    TF.Vx = (TF.Vx - min(TF.Xs(:))) / (max(TF.Xs(:)) - min(TF.Xs(:))) * N + 0.5;
    TF.Vy = (TF.Vy - min(TF.Ys(:))) / (max(TF.Ys(:)) - min(TF.Ys(:))) * M + 0.5;
    TF.Vz = (TF.Vz - min(TF.Zs(:))) / (max(TF.Zs(:)) - min(TF.Zs(:))) * K + 0.5;
    TF.Xs = (TF.Xs - min(TF.Xs(:))) / (max(TF.Xs(:)) - min(TF.Xs(:))) * N + 0.5;
    TF.Ys = (TF.Ys - min(TF.Ys(:))) / (max(TF.Ys(:)) - min(TF.Ys(:))) * M + 0.5;
    TF.Zs = (TF.Zs - min(TF.Zs(:))) / (max(TF.Zs(:)) - min(TF.Zs(:))) * K + 0.5;
elseif ~strcmp(MODE, 'absolute')
    warning('The mode should be either relative(default), of absolute')
end

[A, B, C] = size(TF.Xs);

for a = 1 : A - 1
    for b = 1 : B - 1
        for c = 1 : C - 1
            
            %Selecting a subregion
            TF_CELL.Vx = TF.Vx(a : a + 1, b : b + 1, c : c + 1);
            TF_CELL.Vy = TF.Vy(a : a + 1, b : b + 1, c : c + 1);
            TF_CELL.Vz = TF.Vz(a : a + 1, b : b + 1, c : c + 1);
            TF_CELL.Xs = TF.Xs(a : a + 1, b : b + 1, c : c + 1);
            TF_CELL.Ys = TF.Ys(a : a + 1, b : b + 1, c : c + 1);
            TF_CELL.Zs = TF.Zs(a : a + 1, b : b + 1, c : c + 1);
            
            [X, Y, Z] = meshgrid(...
                round(max(1, TF_CELL.Xs(1))) : round(min(N, TF_CELL.Xs(end))), ...
                round(max(1, TF_CELL.Ys(1))) : round(min(M, TF_CELL.Ys(end))), ...
                round(max(1, TF_CELL.Zs(1))) : round(min(K, TF_CELL.Zs(end))));
            
            %Interpolating voxel origin coordinates
            coord = [X(:), Y(:), Z(:), X(:) .* Y(:), ...
                Y(:) .* Z(:), Z(:) .* X(:), ...
                X(:) .* Y(:) .* Z(:), ones(length(X(:)), 1)];
            
            tform = [TF_CELL.Xs(:), TF_CELL.Ys(:), TF_CELL.Zs(:), ...
                TF_CELL.Xs(:) .* TF_CELL.Ys(:), ...
                TF_CELL.Ys(:) .* TF_CELL.Zs(:), ...
                TF_CELL.Zs(:) .* TF_CELL.Xs(:), ...
                TF_CELL.Xs(:) .* TF_CELL.Ys(:) .* TF_CELL.Zs(:), ...
                ones(8, 1)] \ [TF_CELL.Vx(:), TF_CELL.Vy(:), TF_CELL.Vz(:)];
            
            coord_new = coord * tform;
            
            Xn = coord_new(:, 1); Xn = Xn(:);
            Yn = coord_new(:, 2); Yn = Yn(:);
            Zn = coord_new(:, 3); Zn = Zn(:);
            
            %Selecting voxels with existing origins
            mask = (Xn >= 1) .* (Xn <= N) .* (Yn >= 1) .* ...
                (Yn <= M) .* (Zn >= 1) .* (Zn <= K);
            ind = find(mask);
            
            Xn = Xn(ind);
            Yn = Yn(ind);
            Zn = Zn(ind);
            
            %Interpolating voxel intensities
            flXn = floor(Xn);
            flYn = floor(Yn);
            flZn = floor(Zn);
            
            IMOUT_VECTOR = zeros(size(X));
            
            IMOUT_VECTOR(ind) = ...
                IM(sub2ind([M, N, K], flYn, flXn, flZn)) .* ...
                (flYn + 1 - Yn) .* (flXn + 1 - Xn) .* (flZn + 1 - Zn) + ...
                IM(sub2ind([M, N, K], min(flYn + 1, M), flXn, flZn)) .* ...
                (Yn - flYn) .* (flXn + 1 - Xn) .* (flZn + 1 - Zn) + ...
                IM(sub2ind([M, N, K], flYn, min(flXn + 1, N), flZn)) .* ...
                (flYn + 1 - Yn) .* (Xn - flXn) .* (flZn + 1 - Zn) + ...
                IM(sub2ind([M, N, K], flYn, flXn, min(flZn + 1, K))) .* ...
                (flYn + 1 - Yn) .* (flXn + 1 - Xn) .* (Zn - flZn) + ...
                IM(sub2ind([M, N, K], min(flYn + 1, M), min(flXn + 1, N), flZn)) .* ...
                (Yn - flYn) .* (Xn - flXn) .* (flZn + 1 - Zn) + ...
                IM(sub2ind([M, N, K], flYn, min(flXn + 1, N), min(flZn + 1, K))) .* ...
                (flYn + 1 - Yn) .* (Xn - flXn) .* (Zn - flZn) + ...
                IM(sub2ind([M, N, K], min(flYn + 1, M), flXn, min(flZn + 1, K))) .* ...
                (Yn - flYn) .* (flXn + 1 - Xn) .* (Zn - flZn) + ...
                IM(sub2ind([M, N, K], min(flYn + 1, M), min(flXn + 1, N), ...
                min(flZn + 1, K))) .* (Yn - flYn) .* (Xn - flXn) .* (Zn - flZn);
            
            IMOUT(max(1, round(TF_CELL.Ys(1))) : min(M, round(TF_CELL.Ys(end))), ...
                max(1, round(TF_CELL.Xs(1))) : min(N, round(TF_CELL.Xs(end))), ...
                max(1, round(TF_CELL.Zs(1))) : min(K, round(TF_CELL.Zs(end)))) = ...
                IMOUT_VECTOR;
        end
    end
end

end
