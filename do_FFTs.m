function [mrange,freqs,magZ] = do_FFTs(spokeSurface,Nchop,framerate,plotfft)

%% Initialize

% close all
if nargin ==0
    fn = 'C:\Users\Michael\Documents\PEPL\Spoke Analysis\SampleSpokeSurface_H6-600V.mat';
    load(fn)
    framerate = ImageInfo.Fps;
    NumAzBins = length(AzBinPixels);

   Nchop = 1;
   plotfft = 0;
   clc; close all
end
numthetabins = size(spokeSurface,2);
Nfiles = size(spokeSurface,1);

dt = Nfiles/framerate;
mrange = 1:16;
% modeamps = zeros(1,mrange);

% number of files per FFT analysis
N = floor(Nfiles/Nchop);

% 1/2 that number because only half the FFT matrix is non-redundant
Mfreq = floor(N/2);

% 1/2 the theta bins for the same reason
Mtheta = floor(numthetabins/2);
Z=zeros(N,numthetabins);


% return
%% Calculate the 2D FFTs

% ----------------- Windowing -----------------------------
% choose windowing type
wintime = window(@hamming,N); % column vector
% using non-rectangular windows in theta gives a spurious m=1 mode, since
% it forces a central peak going to zero on each extreme.  
winspace = window(@rectwin,numthetabins)'; % row vector
% create 2D window
win = wintime*winspace; % matrix with N rows, numthetabins columns

% ---------------------- FFT Calculation ---------------------------------
% for i =1:Nchop; Z(:,:,i) = fft2(spokeSurface((i-1)*N+1:i*N,:)); end
for i = 0 : (Nchop-1)
    Zi = fft2( spokeSurface( i*N + 1 : (i+1)*N , :) .*win ) ;
    Z = Z + Zi; 
end



    Z = fft2( spokeSurface );
    N = size(Z,1);
    M = size(Z,2);
    magZ = Z.*conj(Z);
   
    FullRMS = sqrt(sum(magZ(:)))/(M*N)
    DC = sqrt(sum(magZ(1,1)))/(M*N)
    Space = sqrt(sum(magZ(1,2:end)))/(M*N)
    Time = sqrt(sum(magZ(2:end,1)))/(M*N)
    Spacetime = sqrt(sum(sum(magZ(2:end,2:end))))/(M*N)
    
    % Convert to power-spectral density real values from complex FFT
% coefficients. Division by both Nfiles and numthetabins sets the units to
% amplitude^2 / (bin-Hertz) 
magZ = magZ/(dt*Nfiles*numthetabins);

if plotfft
figure()
pcolor(log(magZ)); shading interp;
xlabel('Azimuthal Mode Number')
ylabel('Frequency, Hz')
title('Full 2D DFT')
end

% eliminate redundant second half of FFT matrix (theta) and above Nyquist
% values
magZ = magZ(1:Mfreq,1:Mtheta);

% Make sure all values returned are positive.  This line of code is
% probably redundant...
magZ = abs(magZ(:,mrange));

% return
%% Plot FFT

freqs = (1:Mfreq)*framerate*Nchop/Nfiles;
if plotfft
    % Make surf plot of 2D FFT
    figure('WindowStyle','docked')
    pcolor(mrange-1,freqs,log(magZ));
    shading interp    
end
% return

%% Find peaks

thisFile = mfilename;
isaSubroutine = CheckIfSubroutine(thisFile);

if ~isaSubroutine
Zsearch = magZ;
Npeaks = 5;
freqPeaks = zeros(1,Npeaks);
modePeaks = zeros(1,Npeaks);
ampPeaks  = zeros(1,Npeaks);
 for i=1:Npeaks
     [freqPeaks(i) modePeaks(i) ampPeaks(i) Zsearch] = peakSearchAndCrop(Zsearch);
 end

 Modes = modePeaks-1
ModeFreqs = freqs(freqPeaks)
Ratio = ampPeaks / max(ampPeaks(:))

figure()
surf(mrange-1,freqs,log(magZ));shading interp; view(55,45)
hold on
plot3(Modes,ModeFreqs,log(ampPeaks),'bx','MarkerSize',20)
hold off
end



function [freq mode peakamp Z] = peakSearchAndCrop(Z)

peakamp = max(Z(:));
[freq mode] = ind2sub(size(Z),find(Z==peakamp))

DeltaF = 500;
Z(max(freq-DeltaF,1):min(freq+DeltaF,size(Z,1)),:) = 0;
Z(:,mode) = 0;