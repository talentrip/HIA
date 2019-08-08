function AzimProfile = thetafnX(AzBinPixels,NumAzBins,PixelsPerAzBin,data)

% Preallocate matrix
AzimProfile = zeros(1,NumAzBins);

% Sum up all the pixels in each azimuthal bin
for i = 1:NumAzBins
    region = data(AzBinPixels(i).indices);
    AzimProfile(i) = sum(region);
end

% Now find the average pixel value in each bin by dividing the sum by the
% number of pixels in each bin
AzimProfile = AzimProfile(:)./PixelsPerAzBin(:);
% AzimProfile是一个1*180的向量存储不同周向位置的平均像素值