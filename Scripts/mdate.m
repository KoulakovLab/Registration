function [AGES_ADJ, DIST] = mdate(IMs, AGES)
%MDATE establishes adjusted ages of the brain samples using correlations
%
%   Parameters
%   ==========
%   IMs             - list (of the 3D images, e.g. {IMs{i}P, IM2P, ...})
%   AGES            - array (of experimentally recorded ages of samples)
%   AGES_ADJ        - array (of the adjusted ages of samples)
%   DIST            - array (of correlation distances between samples)
%
%   Author
%   ======
%   Sergey Shuvaev, 2014-2021. sshuvaev@cshl.edu

%Filter images, same as throughout the alignment
for i = 1 : length(IMs)
    IMcontour = mlog(IMs{i}, 10, 5);
    IMbinarized = (IMs{i} > 0);
    IMs{i} = 0.5 * IMbinarized / max(IMbinarized(:)) + ...
        0.5 * IMcontour / mean(IMcontour(IMcontour > 0.01));
    IMs{i} = msubr(IMs{i}, 8);
end

%Correlation-based distance matrix
DIST = 0.5 * eye(length(IMs));
for i = 1 : length(IMs)
    for j = i + 1 : length(IMs)
        DIST(i, j) = corr(IMs{i}(:), IMs{j}(:));
    end
end
DIST = DIST + DIST';
DIST = 1 - DIST;

%Computing the adjusted ages
score = cmdscale(DIST);
ages_embedding = score(:, 1)';
p = polyfit(ages_embedding, AGES, 1);
AGES_ADJ = real(p(1) * ages_embedding + p(2));

%Plot
figure, imagesc(DIST); axis square; colormap parula; colorbar
figure, plot(ages_embedding, AGES, 'r+'), hold on;
plot(ages_embedding, AGES_ADJ, 'go');

end
