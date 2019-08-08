function [ImageInfo, ImageMatrix] = ImportMatVideo(filename,CreateMRAW)


%% Import MATLAB data file and CIH header file
% Loading in the .mat file provides:
%   HiIndices       double
%   HiValues        int16
%   ImageInfo       struct
%   ImageMatrixAC   int8
%   LoIndices       double
%   LoValues        int16
%   MeanImage       int16
%   Normalization   double
%   FileInfo        struct

%% Initialization
if nargin==0
    [file, path] = uigetfile('.mat');
    filename = fullfile(path,file);
    CreateMRAW = 0;
elseif nargin==1
        CreateMRAW = 0;

end

%% Decompression
% MAT file will contain compression algorithm information in variable
% 'FileInfo.CompressionVersion'
load(filename)

% If this file is parent function, show contents of just-loaded file
thisFile = mfilename;
isaSubroutine = CheckIfSubroutine(thisFile);
if ~isaSubroutine; whos; end

if strcmp(FileInfo.CompressionVersion,'1.1')
    % Reconstruct Image Matrix
    ImageMatrix = zeros(size(ImageMatrixAC),'int16');
    ImageMatrixAC = int16(ImageMatrixAC);
    ImageMatrixAC(LoIndices) = ImageMatrixAC(LoIndices) + LoValues;
    ImageMatrixAC(HiIndices) = ImageMatrixAC(HiIndices) + HiValues;
    
    % Do the data de-compression using Haar transforms over image rows and columns
    % (but not in the third dimension, i.e., time)
    ImageMatrixAC = VideoDecompress_v1p1(ImageMatrixAC);
    
    P = size(ImageMatrix,3);
    for i=1:P
        if mod(i,round(P/10)) == 0; disp(num2str(i)); end
        % This technique preserves integer mathematics reversibility
        ImageMatrix(:,:,i) = ImageMatrixAC(:,:,i) + MeanImage/Normalization(i);
        
        % Doing it this way introduces reconstruction errors:
        % ImageMatrix(:,:,i) = (ImageMatrixAC(:,:,i) + MeanImage) / Normalization(i);
    end
    
elseif strcmp(FileInfo.CompressionVersion,'none')
    % Usually, do nothing; ImageMatrix is uncompressed and you're good to go.
    
    % In some cases, if loading a small snippet of a large video, the image
    % matrix is saved as 'ImageMatrixSmall'.  Cover these cases:
    if isempty(whos('ImageMatrix')) && ~isempty(whos('ImageMatrixSmall'))
       ImageMatrix = ImageMatrixSmall; 
    end
else
    disp('Error, unrecognized compression version!')
    disp('')
    return
end

%% Create MRAW File if desired
if CreateMRAW
    ExportMrawVideo(filename,ImageMatrix,ImageInfo);
    
end