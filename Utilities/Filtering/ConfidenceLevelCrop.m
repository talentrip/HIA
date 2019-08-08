function [F,P,freqbins] = ConfidenceLevelCrop(P,F,freqbins,xlimskHz)
%% Comments
%
% ConfidenceLevelCrop takes the frequency vector F and the power spectral
% density P and crops both down to a desired range of frequencies specified
% by xlimskHz.  It also chops the freqbins vector to reflect the newly
% shortened P and F.

%% Code
% F = frequency vector output of an fft call
% P = power spectral density output of an fft call

% numptsBeforeCrop = numel(F);

throwout = F>xlimskHz(2);
P(throwout)=[];
F(throwout)=[];
% numptsAfterCrop = numel(F);

%     % Now downsample the points that are to be plotted to a maximum fo
%     % 50000 points
% xtrim = ceil(numptsAfterCrop/10000);


freqbins(freqbins>xlimskHz(2)) = [];
freqbins(freqbins<xlimskHz(1)) = [];