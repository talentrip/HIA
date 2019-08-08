function CircleFitParams = AutoCircleFit(MeanImage,Thruster,visualize)

if nargin<3
    visualize = 0;
end
if nargin<2
   Thruster = 'unknown'; 
end
if nargin==0
    clc;
    visualize = 1;
  load('High Speed Image Analysis\Sample Data Files\SampleMeanImage.mat');
end
MeanImage = double(MeanImage);

if strcmp(Thruster,'H6')
    [a b R] = TaubinSVDThresholdFit(MeanImage);
    center = [a b];
    IR = 0.75*R;
    OR = 1.25*R;
elseif strcmp(Thruster,'X2SingleInner')
    [a b R] = TaubinSVDThresholdFit(MeanImage);
    center = [a b];
    IR = 0.8*R;
    OR = 1.18*R;
elseif strcmp(Thruster,'X2SingleOuter')
    [a b R] = TaubinSVDThresholdFit(MeanImage);
    center = [a b];
    IR = 0.9*R;
    OR = 1.1*R;
elseif strcmp(Thruster,'X2DualInner')
    [a b R] = TaubinSVDThresholdFit(MeanImage);
    center = [a b];
    IR = 0.35*R;
    OR = 0.55*R;
elseif strcmp(Thruster,'X2DualOuter')
    [a b R] = TaubinSVDThresholdFit(MeanImage);
    center = [a b];
    IR = 0.97*R;
    OR = 1.18*R;
elseif strcmp(Thruster,'Default75-125')
    [a b R] = TaubinSVDThresholdFit(MeanImage);
    center = [a b];
    IR = 0.75*R;
    OR = 1.25*R;
    disp('')
    disp(' ---------- WARNING! ---------- ')
    disp('Default channel ID applied; region between 75% and 125% of calculated radius will be analyzed')
    disp(' ------------------------------ ')
    disp('')
    visualize = 1;
else % if strcmp(Thruster,'ManualFit')  
    [center R IR OR] = CircleFitGui(MeanImage);
    a = center(1);
    b = center(2);
end

CircleFitParams.CircleCenter = center;
CircleFitParams.R = R;
CircleFitParams.IR = IR;
CircleFitParams.OR = OR;

%% Plot results to check if results are reasonable
if visualize
    h = pcolor(MeanImage); shading interp;  daspect([1 1 1]);
%     get(h);
    hold on

    t = 0:0.01:2*pi;%       IR = .8*R;          OR = 1.25*R;
    xi = IR*cos(t) + a;    yi = IR*sin(t) + b;
    xo = OR*cos(t) + a;    yo = OR*sin(t) + b;
    xr = R*cos(t) + a;    yr = R*sin(t) + b;
    plot(xi,yi,'r')
    plot(xo,yo,'r')
    plot(xr,yr,'k')
        plot3(a,b,1e6,'y.','MarkerSize',16)
    hold off
end
