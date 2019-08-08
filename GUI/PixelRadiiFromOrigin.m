function [r theta] = PixelRadiiFromOrigin(PixelMatrix,x0,y0)

[y x] = ind2sub(size(PixelMatrix),1:numel(PixelMatrix));

r = sqrt((y-y0).^2+(x-x0).^2);        % radius of each pixel to origin (x0, y0)
r = reshape(r, size(PixelMatrix));

theta = 180/pi*atan2(y-y0,x-x0);
theta = reshape(theta, size(PixelMatrix));