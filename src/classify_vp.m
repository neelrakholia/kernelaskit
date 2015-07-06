% computes the left, the right data, and the radius
% data:    data values
% nsize:   number of entries
% sigma:   bandwidth
function [datal, datar, rad, cent] = classify_vp(data, nsize, sigma)
% select a random point
perm = randperm(nsize);
index = perm(1);
cent = data(:, index);

% linear search through all the points
dist = distk(data, cent, sigma);

% compute midpoint
mid = floor(nsize/2);

% sort distance and find the elements in each node
[m,ind] = sort(dist);
datal = data(:,ind(1:mid));
datar = data(:,ind(mid+1:end));

% compute radius
rad = distk(cent, data(:, mid), sigma);
    
end