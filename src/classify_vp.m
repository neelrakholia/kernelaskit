% computes the left, the right data, and the radius
% data:    data values
% nsize:   number of entries
% sigma:   bandwidth
% diseval: number of distance evaluations
function [datal, datar, indl, indr, rad, cent, diseval] = classify_vp(data, ...
    indi, nsize, sigma, diseval)
% the number of random points to select
rand = 1; 

% update distance evaluations
diseval = diseval + rand*nsize;

% select random points
perm = randperm(nsize);
index = perm(1:rand);

% linear search through all the points
cent = data(:, index(rand));
bestp = cent;
dist = distk(data, bestp, sigma);

% compute midpoint
mid = ceil(nsize/2);

% sort distance and find the elements in each node
[~,ind] = sort(dist);
datal = data(:,ind(1:mid));
indl = indi(ind(1:mid));
datar = data(:,ind(mid+1:end));
indr = indi(ind(mid+1:end));
cent = bestp;

% compute radius
rad = distk(bestp, datal(:, end), sigma);
    
end