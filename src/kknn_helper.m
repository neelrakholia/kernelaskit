% helper function for kknn
% p:      a row of the cell
% test:   a row of previously discovered neighbors
function [nn, deval] = kknn_helper(p, test, train, query, sigma, K)
p = cell2mat(p);
test = cell2mat(test);
search_inds = unique([p, test]);
nn = kknn(train, search_inds, query, sigma, K, ...
    numel(search_inds));
deval = numel(search_inds);
end