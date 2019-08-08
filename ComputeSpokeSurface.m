function [spokeSurface,spokeSurfFilename] = ComputeSpokeSurface(fn,ImageMatrix,ImageInfo,Thruster)
%% Comments
%
% function [CMSUM,spokeSurface,THETAS,NumFrames,PixelsPerAzBin] =
% IVSTHETAVSTIMESURF(FILEDIR,FORMAT,INTCAT,NumAzBins]
%
%   This function takes images of Hall thruster firings at high speed and
%   collapses them radially to create a set of aggreage brightnesses as a
%   function of azimuthal position, which it then displays in a surface
%   plot showing the time evolution.  This is designed to track spoke
%   instabilities, which should show up as parallel ridges in the final
%   surface plot.
%
%   Specifically, the program accepts files of type FORMAT located in
%   FILEDIR.  It subtracts out the average brightness from each frame, then
%   plots the fluctuating (i.e., AC-coupled) component of the pixel
%   intensities as a function of theta.  Theta resolution is dictated by
%   NumAzBins.  Manual input of the thruster center pixel (hard-coded)
%   is required, as is manual location of the inner and outer channel
%   radius (RMIN and RMAX, in the code).
%
%   by Michael McDonald
%
%   Version history:
%   v1p0: 11/17/09
%   v3p5: 01/25/10
%   v4p0: 4/30/10
%       - updated syntax to match GUI, revamped processing to use bsxfun
%       singleton expansion, removed save-as-you-go feature since
%       processing is now much faster

%% Initialize

cores = feature('Numcores');

% Switches
showplots = 0;
% Have to AC Couple video (subtract off mean) or else zero-frequency
% component swamps everything else.
ACCoupleVideo =0;
IsolateChannel = 1;
AAFilter = 0;
NumAzBins = 180;
NormalizeVideo =0;
saveresults = 1;
ParallelizeSpokeSurfaceComputation = 0;

% Memory Management: judge max sizes for arrays based on memory available
% on computer.
% Check memory available
mem = memory;
% Set max number of array elements to permit binary singleston expansion
% ('bsxfun') based on the max possible array size returned by the 'memory'
% command
BSXsize = 1/10 * mem.MaxPossibleArrayBytes;

%% Load data into an array

tic
t0 = toc;

fn = RemoveFileExtension(fn);
if isempty(whos('ImageMatrix'))
    disp('Importing .mat data file...')
    [ImageInfo,ImageMatrix] = ImportMatVideo(fn);
    NumFrames = size(ImageMatrix,3);
    disp(['Done importing ' num2str(NumFrames) ' images with ' num2str(numel(ImageMatrix)) ' elements.'])
    toc
    disp(' ')
end

% ImageInfo
fps = ImageInfo.Fps;

NumFrames = size(ImageMatrix,3);
MeanImage = mean(ImageMatrix,3);

%% Calculate center

% Automatic Circle Fitting to bright pixels in image (assumed to be
% discharge channel)
CircleFitParams = AutoCircleFit(MeanImage,Thruster,1);

x0 = CircleFitParams.CircleCenter(1);
y0 = CircleFitParams.CircleCenter(2);
R = CircleFitParams.R;
r = PixelPolarCoords(MeanImage,CircleFitParams.CircleCenter);
% 此函数判断通道内外径
RFilter = r>CircleFitParams.IR & r<CircleFitParams.OR;

%% Assign azimuthal bins.
% ChannelThetaBin函数规定了周向角度分配个数，将周向隔离出的放电部分分为bins个模块，结构体AzBinPixels为bin中像素值大于0的indice值
% PixelsPerAzBin为每个bin中像素值个数，即矩阵元素个数，
disp('Assigning azimuthal bins...')
[r,XXX,AzBinPixels,AzBinAngles,PixelsPerAzBin] = ChannelThetaBin(MeanImage,[x0 y0],RFilter,NumAzBins);
disp('Done assigning azimuthal bins.')
toc
pause(.1)
disp(' ')

%% Isolate Channel

if IsolateChannel
    disp('Isolating Channel...')
    x = int16(RFilter);
    if numel(ImageMatrix)>BSXsize
        for i=1:NumFrames; ImageMatrix(:,:,i) = ImageMatrix(:,:,i).*x; end
    else
        ImageMatrix = bsxfun(@times,ImageMatrix,x);
    end   
    % This quantity is saved later -- disregard red MLint underlining
    IsolatedImage = mean(ImageMatrix,3);   
    disp('Done isolating channel.')
    toc
    pause(1)
    disp(' ')
end

%% Anti-aliasing filter, if desired

if AAFilter
    numel(MeanImage)
    sum(RFilter(:))
    disp('Anti-aliasing filter in progress...')
    [B,A] = butter(9,.95);
    for j = 1:size(ImageMatrix,1)
        j
        parfor k=1:size(ImageMatrix,2)
            if RFilter(sub2ind(size(MeanImage),j,k))
                
                ImageMatrix(j,k,:) = filtfilt(B,A,double(ImageMatrix(j,k,:)));
            else
            end
        end
    end
    toc
end

%% Normalize summed intensity for breathing mode

if NormalizeVideo
    disp('Normalizing out breathing mode...')
    
    FrameIntensity = sum(sum(ImageMatrix));
    if numel(ImageMatrix)>BSXsize
        Normalization = single(mean(FrameIntensity(:))/FrameIntensity);
        ImageMatrix = single(ImageMatrix);
        for i = 1:NumFrames; ImageMatrix(:,:,i) = ImageMatrix(:,:,i)*Normalization(i); end
    else
        Normalization = mean(FrameIntensity(:))/FrameIntensity;
        ImageMatrix = bsxfun(@times,double(ImageMatrix),Normalization);
    end
    MeanImage = mean(ImageMatrix,3);
    disp('Done normalizing.')
    toc
    disp(' ')
end

%% Subtract off mean image to look at AC-coupled images

if ACCoupleVideo
    disp('Subtracting off mean image...')
    if NormalizeVideo
        x = single(MeanImage);
    else
        x = int16(MeanImage);
    end
    if numel(ImageMatrix)>BSXsize
        for i=1:NumFrames; ImageMatrix(:,:,i) = ImageMatrix(:,:,i) - x; end
    else
        ImageMatrix = bsxfun(@minus,double(ImageMatrix),x);
        %ImageMatrix = bsxfun(@minus,ImageMatrix,x);
    end
    
    disp('Done subtracting mean image.')
    toc
    disp(' ')
    
end

%% Compute Spoke Surface
% Add up (or just average) the intensities in each bin and use these averaged- or summed-over-bin
% intensities to collapse a 2D matrix of intensity vs. r and theta into a 1-D vector of
% average or summed intensities vs. theta.

% count # of pixels in each bin to make the averages more accurate.

disp('Computing spoke surface...')
spokeSurface = zeros(NumFrames,NumAzBins);
% 多核运算
if ParallelizeSpokeSurfaceComputation
    cores = feature('NumCores');
    isOpen = matlabpool('size') > 0;
    if ~isOpen
        matlabpool open
    end
    disp(['Using ' num2str(cores) ' cores...'])
    % parfor 并行运算for循环
    parfor(j = 1:NumFrames,cores)        
        image = ImageMatrix(:,:,j); % 某一帧图像矩阵
        [spokeSurface(j,:)] = thetafnX(AzBinPixels,NumAzBins,PixelsPerAzBin,image);
    end
else
    for j=1:NumFrames
        if mod(j,round(NumFrames/4)) == 0; disp(num2str(j)); end
        image = ImageMatrix(:,:,j);
        [spokeSurface(j,:)] = thetafnX(AzBinPixels,NumAzBins,PixelsPerAzBin,image);
    end
end
% 在这设断点，得到中间变量spokeSurface ，为周向分配后计算每一bin模块的平均像素值，一个角度对应一个值
% 是一个帧数*bins数的二维矩阵
disp('Done computing spoke surface.')
toc
disp(' ')

%% Save the results of the spoke surface matrix

% Need to add the capability to save normalized spoke surfaces here, and
% make the naming a little less ridiculous with if loops.  To do this, make
% a string for each norm/AC/isolate case, and just append them in one
% naming line
if saveresults

    if IsolateChannel
        if ACCoupleVideo
            spokeSurfFilename = [fn '-spokeSurface-IC1-N0-AC1-' Thruster];
            save(spokeSurfFilename,'spokeSurface','MeanImage','IsolatedImage','AzBinPixels','ImageInfo','CircleFitParams');
        else
            spokeSurfFilename = [fn '-spokeSurface-IC1-N0-AC0-' Thruster];
            save(spokeSurfFilename,'spokeSurface','MeanImage','IsolatedImage','AzBinPixels','ImageInfo','CircleFitParams');
        end
    else
        if ACCoupleVideo
            spokeSurfFilename = [fn '-spokeSurface-IC0-N0-AC1-' Thruster];
            save(spokeSurfFilename,'spokeSurface','MeanImage','AzBinPixels','ImageInfo','CircleFitParams');
        else
            spokeSurfFilename = [fn '-spokeSurface-IC0-N0-AC0-' Thruster];
            save(spokeSurfFilename,'spokeSurface','MeanImage','AzBinPixels','ImageInfo','CircleFitParams');
        end
    end
end

%% Plot the results

if(1)
    figure('Name',['Spoke Surface - ' Thruster],'NumberTitle','off')

    t = (0:(NumFrames-1))/fps;
        npts = round(fps/1e3);
    
    surf(AzBinAngles,t(1:npts),spokeSurface(1:npts,:))
    shading interp
    title('Spoke Surface')

    xlabel('Azimuthal Angle, deg')
xlim([-180 180])
ylabel('Time, seconds')
ylim([0 1e-3])
zlabel('Average Pixel Value')
    
end
tf = toc;
disp(['Total time elapsed: ' num2str(tf-t0) ' seconds.'])
% beep

