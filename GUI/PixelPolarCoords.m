function [r theta] = PixelPolarCoords(PixelMatrix,Center)

% This function converts Cartesian coordinates of a video frame (the rows
% and columns, given by the indices of each cell via 'ind2sub') into polar
% coordinates.
%
% As a result, this subroutine is responsible for identifying the
% rotational origin of the polar coordinate system, i.e., the place where
% theta = 0.  
%
% To match with normal convention, the zero angle should be at the 3
% o'clock position, with increasing angles moving counterclockwise.  Note
% that older versions of this code did not carefully account for this, and
% as a result used the 9 o'clock position (still counterclockwise though).  
%
% Since the element (1,1) of a matrix 

if nargin==0
    % Load a sample mean image file and use the image center as the circle
    % center
       load('High Speed Image Analysis\Sample Data Files\SampleMeanImage.mat');
       Center = [round(size(MeanImage,1)/2) round(size(MeanImage,2)/2)];
       
       PixelMatrix = MeanImage;
end
x0 = Center(1);
y0 = Center(2);

[y x] = ind2sub(size(PixelMatrix),1:numel(PixelMatrix));

r = sqrt((y-y0).^2+(x-x0).^2);        % radius of each pixel to origin (x0, y0)
r = reshape(r, size(PixelMatrix));

% The y-values are given a negative sign when computing theta because
% MATLAB uses a convention where the top left corner of the image is the
% origin -- thus, the y-axis is inverted.  To compute sensible azimuthal
% angles following the convention that 0 is at 3 o'clock and positive
% angles go counterclockwise, one must have the extra negative sign.  For
% a demonstration, run 'ChannelThetaBin.m' with no input arguments
% (nargin=0)
theta = 180/pi*atan2(-(y-y0),x-x0);
theta = reshape(theta, size(PixelMatrix));


if nargin==0
    % Now plot the values of theta to illustrate the origin and
    % directionality (clockwise vs. counterclockwise)
    
    pcolor(theta); shading flat; daspect([1 1 1]); colorbar; colormap jet;
    clear all
end