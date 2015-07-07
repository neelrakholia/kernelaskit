% computes the left, the right data, and the radius
% data:    data values
% nsize:   number of entries
% sigma:   bandwidth
function [datal, datar, rad, cent] = classify_vp(data, nsize, sigma)
% the number of random points to select
rand = 20;

% select random points
perm = randperm(nsize);
index = perm(1:rand);

varm = 0;
bestp = 0;

% select best point
for i = 1:rand
    cent = data(:, index(i));
    dist = distk(data, cent, sigma);
    
    % calculate variance
    varl = var(dist);
    
    if(varl > varm)
        bestp = cent;
        varm = varl;
    end
end

% linear search through all the points
dist = distk(data, bestp, sigma);

% compute midpoint
mid = floor(nsize/2);

% sort distance and find the elements in each node
[m,ind] = sort(dist);
datal = data(:,ind(1:mid));
datar = data(:,ind(mid+1:end));

% compute radius
rad = distk(bestp, data(:, mid), sigma);
    
end