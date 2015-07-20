% evaluates the rbf kernel distance
function d = rbf(X, x, sigma)
% evaluate euclidean part
[~,N] = size(X);
[~,n] = size(x);
X2 = sum(X.^2,1);
x2 = sum(x.^2,1);
dotProd = X'*x;
d = repmat(x2,N,1)+repmat(X2',1,n)-2*dotProd;

% evaluate kernel
d = exp(-d/(2*sigma^2));
d = sqrt(-2*d + 2);
d = d./(d + 1);
end