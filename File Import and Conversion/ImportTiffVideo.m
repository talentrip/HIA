% function do_importfiles_multipage(filename)
function [ImageInfo, ImageMatrix] = ImportTiffVideo(fn,saveall,savefirstN,smallNfrm,CompressData)
%% Comments
%


%% Initialization


textout = 1;

if nargin==0
    [filename, pathname] = uigetfile( ...
        {'*.tif','Multipage TIFF files (*.tif)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file');
    fn = [pathname filename]
    savefirstN = 1;
    smallNfrm = 100;
    close all
    saveall = 1;
    CompressData = 0;
end

if CompressData
    FileInfo.CompressionVersion = '1.1';
    FileInfo.DeleteMRAW = 0;
    CompressFilename = [RemoveFileExtension(fn) '_Compress'];
else
    FileInfo.CompressionVersion = 'none';
    FileInfo.DeleteMRAW = 0;
end

if textout; tic; end

%% Load the files

% If there is no *.mat file or if the user requests a forced re-load of
% the original images (if, for instance, the *.mat file has been
% corrupted), then load the multipage TIFF file and save as a *.mat
% file for future use.

if textout;     disp(' '); end
if textout;     disp('Reading multi-page TIFF file...')    ; end

% Get info from multi-page TIFF files
info = imfinfo(fn);
numFrames = numel(info); %1000;
% initialize image matrix
ImageMatrix = zeros(info(1).Height, info(1).Width, numFrames,'int16');

% save each TIFF image file into a page of a 3D matrix
if textout;     disp(' ')   ; end
if textout;     disp([num2str(numFrames) '-page multipage TIFF file found, importing all pages...']); end  % status output to user
for i = 1:numFrames
    if mod(i,round(numFrames/4)) == 0; disp(num2str(i)); end
    ImageMatrix(:,:,i) = ImportTiffSinglePage(fn, 'Index', i, 'Info', info);
end
if textout;     toc; end

bitdepth = round(log(double(max(ImageMatrix(:))))/log(2));

% This takes about 150 seconds for a 10,000 page 128x112 pixel image on laptop

% save 'ImageMatrix' to MATLAB binary file for quick future loading


if textout;    toc; end



%% Save First N Frames

    ImageInfo.NumFrames = size(ImageMatrix,3);
    ImageInfo.ImageWidth = size(ImageMatrix,2);
    ImageInfo.ImageHeight = size(ImageMatrix,1);
    ImageInfo.BitDepth = bitdepth;
    
    answer = inputdlg('Please enter video framerate (fps):',...
        'Video Framerate Required',1,{'1000'});
    ImageInfo.Fps = str2num(answer{1});
    
    ImageInfo.Duration = ImageInfo.NumFrames/ImageInfo.Fps;


% Save a subset of the video in a smaller file, if the input length of the
% shorter file (smallNfrm) is less than the length of the full file
if savefirstN && smallNfrm < size(ImageMatrix,3)
    % store compression version for full image matrix
    dum = FileInfo.CompressionVersion;
    FileInfo.CompressionVersion = 'none';
    % Create file name to save first N frames separately
    smallfn = [RemoveFileExtension(fn) '_small.mat'];
    % Store first N images in new matrix
    ImageMatrixSmall = ImageMatrix(:,:,1:smallNfrm);
    % Compute mean image (double floating precision)
    MeanImage = mean(ImageMatrix,3);
    % Save to specified file
    save(smallfn,'ImageInfo','FileInfo','MeanImage','ImageMatrixSmall','-v7.3')
    % Restore compression version for saving large file
    FileInfo.CompressionVersion = dum;
    clear dum
end

if saveall
    
    disp(['All files imported.  Saving... '])
    if CompressData
        disp(['Compressing... '])
        MeanImage = int16(mean(ImageMatrix,3));
        FrameIntensity = sum(sum(ImageMatrix));
        Normalization = mean(FrameIntensity(:))/FrameIntensity;
        ImageMatrixAC = zeros(size(ImageMatrix),'int16');
        
        % Subtract off normalized mean image to create AC-coupled image matrix
        %   -- this step subtracts off a brighter mean image from brighter
        %   images (cases where normalization(i) < 1) and vice versa
        for i=1:size(ImageMatrix,3)
            % Doing it this way preserves integer mathematics reversibility
            ImageMatrixAC(:,:,i) = ImageMatrix(:,:,i) - MeanImage/Normalization(i);
            
            % Doing it this way introduces reconstruction errors:
            % ImageMatrixAC(:,:,i) = ImageMatrix(:,:,i)*Normalization(i) - MeanImage;
        end
        
        % Do the data compression using Haar transforms over image rows and columns
        % (but not in the third dimension, i.e., time)
        [Final8BitSparsity,ImageMatrixAC] = VideoCompress_v1p1(ImageMatrixAC);
        
        % Store indices and excess values of each pixel outside the int8 range
        if Final8BitSparsity<0.01
            LoIndices = sparse(ImageMatrixAC(:)<intmin('int8'));
            HiIndices = sparse(ImageMatrixAC(:)>intmax('int8'));
        else
            LoIndices = single(find(ImageMatrixAC(:)<intmin('int8')));
            HiIndices = single(find(ImageMatrixAC(:)>intmax('int8')));
        end
        LoValues = ImageMatrixAC(LoIndices) - int16(intmin('int8'));
        HiValues = ImageMatrixAC(HiIndices) - int16(intmax('int8'));
        
        % Now store image matrix into 8-bit 'int8' format since the larger values
        % have been removed and stored elsewhere
        ImageMatrixAC = int8(ImageMatrixAC);
        
        % Save compressed image data
        save(CompressFilename,'ImageInfo','FileInfo','ImageMatrixAC','MeanImage','Normalization',...
            'LoIndices','HiIndices','LoValues','HiValues','-v7.3')
        disp(['Saved compressed image matrix as ' CompressFilename])
    else
        % Save uncompressed image data
        save(RemoveFileExtension(fn),'ImageInfo','FileInfo','ImageMatrix','-v7.3')
        disp(['Saved image matrix as ' RemoveFileExtension(fn)])
    end       
end
