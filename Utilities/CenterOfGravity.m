function x0 = CenterOfGravity(x)

% x is a vector or matrix and this function returns the mean weighted value
% of the values in the vector, or equivalently the center of gravity of the
% histogram of values in the vector.

[n,xout] = hist(double(x(:)),round(numel(x)/100));

x0 = sum(n.*xout) / sum(n);