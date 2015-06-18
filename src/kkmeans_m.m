function [label, medoid] = kkmeans_m(points, kernelf, k, n)

% get kernel matrix
K = zeros(n, n);
for i = 1:n
    for j = 1:n
        K(i, j) = kernelf(points(:, i), points(:, j));
    end
end

cluster = kkmedoid(points', k, K);
label = cluster.label;
medoid = cluster.medoids;

end