function TF_BEST = msimanneal(IM1, IM2, AGE_FLAG, SIZE_FLAG)
%MSIMANNEAL Computes the trasformation TF_BEST to align 3D image IM1 to IM2
%using attention-gated simulated annealing procedure 
%
%   Parameters
%   ==========
%   IM1             - array (3D image to be transformed)
%   IM2             - array (target 3D image)
%   AGE_FLAG        - string ('young' | 'adult')
%   SIZE_FLAG       - string ('full' | 'halves')
%   TF_BEST         - structure (best transformation)
%   DIVS            - array (number of cells in transformation grid for
%                            each stage of alignment)
%   TEMP_START      - array (initial temperatures for simulated annealing)
%   TEMP_END        - array (final temperatures for simulated annealing)
%   COEFF_START     - array (initial attemped displacement magnitude,
%                            relative to transformation grid cell size)
%   COEFF_DROP      - double (coefficient drop on an unseccessful attempt)
%   N_STEPS         - array (numbers of simulated annealing steps per stage)
%   SCALE           - array (image downsampling, e.g. 8 -> 1/8 of the
%                            original resolution)
%   LAMBDA          - double (regularization coefficient for deformation
%                             energy)
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

if nargin < 3
    AGE_FLAG = 'young';
end
if nargin < 4
    SIZE_FLAG = 'full';
end

switch SIZE_FLAG
    case 'full'
        DIVS = [1 1 1; 4 2 1; 8 4 2; 8 8 4];
    case 'halves'
        DIVS = [1 1 1; 4 1 1; 8 2 2; 8 4 4];
end
TEMP_START = [1e-3, 1e-3, 1e-4, 1e-5];
TEMP_END = TEMP_START / 30;
COEFF_START = [0.2, 0.2, 0.2, 0.2];
N_STEPS = [2e3, 2e3, 2e3, 2e3];
SCALE = [8, 8, 8, 4];
COEFF_DROP = 0.9966;
switch AGE_FLAG
    case 'young'
        LAMBDA = 1e-3;
    case 'adult'
        LAMBDA = 1e-2;
end

%Preprocess the images
switch AGE_FLAG
    case 'young'
        im1contours = mlog(IM1, 10, 5);
        im2contours = mlog(IM2, 10, 5);
    case 'adult'
        im1contours = double(mlog(IM1, 10, 5) > 0.01);
        im2contours = double(mlog(IM2, 10, 5) > 0.01);
end
im1binarized = (IM1 > 0);
im2binarized = (IM2 > 0);

switch AGE_FLAG
    case 'young'
        im1filtered = 0.5 * im1binarized / max(im1binarized(:)) + ...
            0.5 * im1contours / mean(im1contours(im1contours > 0.01));
        im2filtered = 0.5 * im2binarized / max(im2binarized(:)) + ...
            0.5 * im2contours / mean(im2contours(im1contours > 0.01));
    case 'adult'
        im1filtered = 0.2 * im1binarized / max(im1binarized(:)) + ...
            0.8 * im1contours / mean(im1contours(im1contours > 0.01));
        im2filtered = 0.2 * im2binarized / max(im2binarized(:)) + ...
            0.8 * im2contours / mean(im2contours(im1contours > 0.01));
end

%%%%%%%%%%%%%%%%%%%%%%%
% Simulated annealing %
%%%%%%%%%%%%%%%%%%%%%%%

PHASE_BORDERS = cumsum(N_STEPS);
PHASE_BORDERS = PHASE_BORDERS(1 : end - 1);
PROB = 1 / 8 : 1 / 8 : 1;
COST_LOG = zeros(1, sum(N_STEPS));
TENSION_LOG = zeros(1, sum(N_STEPS));
figure

for st = 1 : length(N_STEPS)
    for nst = 1 : N_STEPS(st)       
        if nst == 1
            
            %Resize images
            IM1 = msubr(im1filtered, SCALE(st));
            IM2 = msubr(im2filtered, SCALE(st));

            [M, N, K] = size(IM1);
            coeff = COEFF_START(st);

            %Define the coordinates of the transformation grid nodes
            Xs = (0.5 : N / DIVS(st, 1) : N + 0.5);
            Ys = (0.5 : M / DIVS(st, 2) : M + 0.5);
            Zs = (0.5 : K / DIVS(st, 3) : K + 0.5);

            [Xs, Ys, Zs] = meshgrid(Xs, Ys, Zs);

            TF_TMP.Xs = Xs;
            TF_TMP.Ys = Ys;
            TF_TMP.Zs = Zs;

            if st == 1 %On the 1st stage, no displacements
                TF_TMP.Vx = Xs;
                TF_TMP.Vy = Ys;
                TF_TMP.Vz = Zs;

                TF_BEST = TF_TMP;
            else %On the later stages, interpolate previous displacements
                SCALE_X = N / max(TF.Xs(:) - 0.5);
                SCALE_Y = M / max(TF.Ys(:) - 0.5);
                SCALE_Z = K / max(TF.Zs(:) - 0.5);
                TF_TMP.Vx = interp3((TF.Xs - 0.5) * SCALE_X + 0.5, ...
                    (TF.Ys - 0.5) * SCALE_Y + 0.5, (TF.Zs - 0.5) * SCALE_Z + 0.5, ...
                    (TF.Vx - 0.5) * SCALE_X + 0.5, TF_TMP.Xs, TF_TMP.Ys, TF_TMP.Zs);
                TF_TMP.Vy = interp3((TF.Xs - 0.5) * SCALE_X + 0.5, ...
                    (TF.Ys - 0.5) * SCALE_Y + 0.5, (TF.Zs - 0.5) * SCALE_Z + 0.5, ...
                    (TF.Vy - 0.5) * SCALE_Y + 0.5, TF_TMP.Xs, TF_TMP.Ys, TF_TMP.Zs);
                TF_TMP.Vz = interp3((TF.Xs - 0.5) * SCALE_X + 0.5, ...
                    (TF.Ys - 0.5) * SCALE_Y + 0.5, (TF.Zs - 0.5) * SCALE_Z + 0.5, ...
                    (TF.Vz - 0.5) * SCALE_Z + 0.5, TF_TMP.Xs, TF_TMP.Ys, TF_TMP.Zs);
            end
            
            TF = TF_TMP;
            IM1_TFORMED = mtransform(IM1, TF);
            TENSION = mregularize(TF);
            COST = corr(IM1_TFORMED(:), IM2(:)) - LAMBDA * TENSION;
            COST_BEST = COST;
        end

        T = TEMP_START(st) * (TEMP_END(st) / (TEMP_START(st) + 1e-30)) ^ ...
            (nst / N_STEPS(st));
        TF_NEW = TF;

        %Set probabilities to choose a grid node for attempted displacement
        %as an L1 norm of diffences between two images in the grid cells 
        %adjacent to that node
        if ((~mod(nst - 1, 100)) && (sum(DIVS(st, :)) > 3))
            IM_DIFF = abs(IM1_TFORMED - IM2);
            XS_INT = min(N, max(1, round(unique(TF.Xs(:)))));
            YS_INT = min(M, max(1, round(unique(TF.Ys(:)))));
            ZS_INT = min(K, max(1, round(unique(TF.Zs(:)))));
            PROB = zeros(DIVS(st, [2 1 3]));
            for ii = 1 : DIVS(st, 1)
                for jj = 1 : DIVS(st, 2)
                    for kk = 1 : DIVS(st, 3)
                        IM_DIFF_CELL = IM_DIFF(YS_INT(jj) : YS_INT(jj + 1), ...
                            XS_INT(ii) : XS_INT(ii + 1), ...
                            ZS_INT(kk) : ZS_INT(kk + 1));
                        PROB(jj, ii, kk) = sum(IM_DIFF_CELL(:));
                    end
                end
            end
            PROB = convn(PROB, ones(2, 2, 2));
            PROB = cumsum(PROB(:) / sum(PROB(:)));
        end

        %Pick a node and attempt a displacement of that node
        NODE_IX = find(rand <= PROB, 1, 'first');
        TF_NEW.Vx(NODE_IX) = TF.Vx(NODE_IX) + randn(1) * coeff * N / DIVS(st, 1);
        switch SIZE_FLAG
            case 'full'
                TF_NEW.Vy(NODE_IX) = ...
                    TF.Vy(NODE_IX) + randn(1) * coeff * M / DIVS(st, 2);
            case 'halves'
                TF_NEW.Vy(NODE_IX) = ...
                    TF.Vy(NODE_IX) + randn(1) * coeff * M / DIVS(st, 2) * ...
                    (TF.Ys(NODE_IX) ~= M);
        end
        TF_NEW.Vz(NODE_IX) = TF.Vz(NODE_IX) + randn(1) * coeff * K / DIVS(st, 3);

        %Extract the region of interest (ROI) consisting of the grid cells
        %adjacent to the displaced node (the other parts of the 3D image
        %are not affected by a displacement of a single node)
        ROI_IX = find(...
            (abs(Xs - Xs(NODE_IX)) < 1.5 * (N - 1) / DIVS(st, 1)) .* ...
            (abs(Ys - Ys(NODE_IX)) < 1.5 * (M - 1) / DIVS(st, 2)) .* ...
            (abs(Zs - Zs(NODE_IX)) < 1.5 * (K - 1) / DIVS(st, 3)));

        %Extract the parts of the transformations covering this ROI
        [sz1x, sz1y, sz1z] = ind2sub(size(TF_NEW.Vx), min(ROI_IX));
        [sz2x, sz2y, sz2z] = ind2sub(size(TF_NEW.Vx), max(ROI_IX));

        TF_NEW_ROI.Xs = zeros(1 + abs(sz2x - sz1x), 1 + abs(sz2y - sz1y), ...
                       1 + abs(sz2z - sz1z));
        TF_NEW_ROI.Ys = TF_NEW_ROI.Xs;
        TF_NEW_ROI.Zs = TF_NEW_ROI.Xs;
        TF_NEW_ROI.Vx = TF_NEW_ROI.Xs;
        TF_NEW_ROI.Vy = TF_NEW_ROI.Xs;
        TF_NEW_ROI.Vz = TF_NEW_ROI.Xs;
        
        TF_ROI = TF_NEW_ROI; %the old transformation, limited to the ROI

        TF_NEW_ROI.Xs(:) = TF_NEW.Xs(ROI_IX);
        TF_NEW_ROI.Ys(:) = TF_NEW.Ys(ROI_IX);
        TF_NEW_ROI.Zs(:) = TF_NEW.Zs(ROI_IX);
        TF_NEW_ROI.Vx(:) = TF_NEW.Vx(ROI_IX);
        TF_NEW_ROI.Vy(:) = TF_NEW.Vy(ROI_IX);
        TF_NEW_ROI.Vz(:) = TF_NEW.Vz(ROI_IX);
        
        TF_ROI.Xs(:) = TF.Xs(ROI_IX);
        TF_ROI.Ys(:) = TF.Ys(ROI_IX);
        TF_ROI.Zs(:) = TF.Zs(ROI_IX);
        TF_ROI.Vx(:) = TF.Vx(ROI_IX);
        TF_ROI.Vy(:) = TF.Vy(ROI_IX);
        TF_ROI.Vz(:) = TF.Vz(ROI_IX);

        %Apply transformation to the ROI and put the transformed ROI on the
        %full image (same effect as transforming the entire 3D image using
        %the entire transformation, but faster)
        IM1_TF_NEW_ROI = mtransform(IM1, TF_NEW_ROI, 'absolute');                    
        IM1_TFORMED_NEW = IM1_TFORMED;
        IM1_TFORMED_NEW(...
            max(1, round(TF_NEW_ROI.Ys(1))) : min(M, round(TF_NEW_ROI.Ys(end))), ...
            max(1, round(TF_NEW_ROI.Xs(1))) : min(N, round(TF_NEW_ROI.Xs(end))), ...
            max(1, round(TF_NEW_ROI.Zs(1))) : min(K, round(TF_NEW_ROI.Zs(end)))) = ...
            IM1_TF_NEW_ROI(...
                max(1, round(TF_NEW_ROI.Ys(1))) : min(M, round(TF_NEW_ROI.Ys(end))), ...
                max(1, round(TF_NEW_ROI.Xs(1))) : min(N, round(TF_NEW_ROI.Xs(end))), ...
                max(1, round(TF_NEW_ROI.Zs(1))) : min(K, round(TF_NEW_ROI.Zs(end))));

        %Compute the change in the deformation energy ('tension') as a
        %result of applying transformation to the ROI
        TENSION_NEW = TENSION - ...
            mregularize(TF_ROI, TF) + mregularize(TF_NEW_ROI, TF);
        COST_NEW = corr(IM1_TFORMED_NEW(:), IM2(:)) - LAMBDA * TENSION_NEW;

        %Accepting of declining the attempted transformation
        if exp((COST_NEW - COST) / T) > rand(1)
            COST = COST_NEW;
            TENSION = TENSION_NEW;
            TF = TF_NEW;
            IM1_TFORMED = IM1_TFORMED_NEW;
            coeff = coeff / COEFF_DROP;

            %Looking whether new transformation is the best
            if COST > COST_BEST
                COST_BEST = COST;
                TF_BEST = TF;
            end
        else
            coeff = coeff * COEFF_DROP;
        end

        %Saving current parameters of simulated annealing
        COST_LOG(sum(N_STEPS(1 : st - 1)) + nst) = COST;
        TENSION_LOG(sum(N_STEPS(1 : st - 1)) + nst) = TENSION;

        %Plotting the parameters
        if ~mod(nst, 200)
            subplot(2, 2, 1), mshow(IM1, IM2); axis image; axis off
            title('Original'); 
            subplot(2, 2, 3), mshow(IM1_TFORMED,IM2); axis image; axis off
            title('Transformed')
            subplot(2, 2, [2; 4]), plot(COST_LOG, 'k'), hold on
            plot(TENSION_LOG, 'k--')
            line([PHASE_BORDERS; PHASE_BORDERS], [0 1], ...
                'Color', 'black', 'LineStyle', '--');
            axis([0, length(COST_LOG), 0, 1])
            legend('Objective function', 'Deformation energy', ...
                'Location', 'West')
            xlabel('Iterations')
            hold off; drawnow
        end
    end
end

end
