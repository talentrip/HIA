function [a,b,R] = TaubinSVDThresholdFit(MeanImage)

%% Comments
%--------------------------------------------------------------------------
% [A B R] = TaubinSVDThresholdFit(MEANIMG)
%
%   This function takes in a rectangular 2D matrix MEANIMG corresponding to
%   an image, with one element per pixel and the value of each element
%   corresponding to the brightness of the pixel.  Given this MEANIMG
%   matrix, the program uses a thresholded Taubin circle-fitting technique
%   to calculate the center (A,B) and radius R of the best-fit circle to
%   the image.
%
%   For more information, see:
%   Chernov, ch. 5
%   G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
%                  Space Curves Defined By Implicit Equations, With
%                  Applications To Edge And Range Image Segmentation",
%      IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
%
%   The actual Taubin fit calculation by singular value decomposition is
%   from code provided by Chernov on his website at
%   (http://people.cas.uab.edu/~mosya/cl/).
%
%   The program is intended for use with images of Hall-effect thrusters to
%   convert rectangular x-y (row-column) pixel coordinates to polar
%   coordinates in terms of r and theta for analysis of rotating
%   instabilities.
%
%   This program is valid for cases even where less than the full 2*pi of
%   the thruster channel is captured in the image, unlike the Kasa fit
%   which tends to bias toward underestimation of the radius.
%
%   Michael McDonald
%   University of Michigan
%   9/10/12
%   v1.0

%% Initialize, hit switches

% The 'iterate' switch triggers a second iteration of circle fitting that
% uses the first iteration to refine the fit by zeroing extraneous pixels
% (especially the cathode) from <0.5*R1 and >1.5*R1, where R1 is the first
% iteration's computed radius.  Its use is recommended.
iterate = 1;
% The 'visualize' switch generates a plot to illustrate the iterative
% refinement process and show the final computed circle fit
visualize = 0;

% For debugging the algorithm, run with no arguments passed in and use the
% sample mean image file included with the software
if nargin==0
    clc; close all;
    load('High Speed Image Analysis\Sample Data Files\SampleMeanImage.mat');
    visualize = 1;    
end

%% Determine Threshold Pixel Intensity Ti and Set Pixels < Ti to Zero

% Compute a "center of mass" pixel intensity using a histogram of all the
% pixel values in the mean image

Ti = CenterOfGravity(MeanImage);

% Compute a Taubin fit to the circle
[a,b,R] = TaubinSVD(MeanImage,Ti);


%% Iteratively refine

if iterate
    % Compute polar coordinates in image given calculated center
    [r theta] = PixelRadiiFromOrigin(MeanImage,a,b);
    % Set pixels far inside and far outside this circle to zero
    MeanImage(r<.5*R) = 0;
    MeanImage(r>1.5*R) = 0;
    % Recalculate circle without these extraneous points affecting things
    [a,b,R] = TaubinSVD(MeanImage,Ti);
end

%% Plot results to illustrate iterative process and check results
if visualize
    VisualizeCircleFit
end

function [a,b,R] = TaubinSVD(MeanImage,Ti)

% Find indices of image pixels above threshold
BrightSpots = find(MeanImage>Ti);

% Use only these bright pixels as inputs to compute the Taubin fit
[y x] = ind2sub(size(MeanImage),BrightSpots) ;

% Center data
centroid = [mean(x) mean(y)];   % the centroid of the data set
X = x - centroid(1);  %  centering data
Y = y - centroid(2);  %  centering data

% Do Taubin SVD method, following Chernov code obtained from
% http://people.cas.uab.edu/~mosya/cl/.  Chernov notes that the code is
% optimized for stability, not speed.  However, it still takes fractions of
% a second to compute, and stability is a good thing for an automated
% procedure.
Z = X.*X + Y.*Y;
Zmean = mean(Z);
Z0 = (Z-Zmean)/(2*sqrt(Zmean));
ZXY = [Z0 X Y];
[U,S,V]=svd(ZXY,0);
A = V(:,3);
A(1) = A(1)/(2*sqrt(Zmean));
A = [A ; -Zmean*A(1)];
Par = [-(A(2:3))'/A(1)/2+centroid , sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];
a = Par(1);
b = Par(2);
R = Par(3);