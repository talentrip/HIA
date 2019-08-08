function IntegratedFASTCAMAnalysis(fn,Thruster)
% This program integrates several steps of FASTCAM video analysis: File
% import (from either MRAW or MAT formats) into an image matrix, computation of a spoke surface
% from the image matrix, calculation of the DFT of the resulting spoke
% surface, and plotting.
%
% Inputs:
% fn = filename string of the video file.  Must end in either '.mraw' or
% '.mat'
%
% Thruster = string giving thruster type to be analyzed.  Must be one of
% the following options: 'H6', 'X2SingleInner', 'X2SingleOuter',
% 'X2DualInner', 'X2DualOuter', or 'Other'.  This input guides the
% automatic circle fitting algorithm to center the annular thruster image
% and analyze only the discharge channel.
%
% If the program is run without inputs, the user will be prompted to supply
% the file and thruster type.

%% Initialization

% Switches and options
saveall = 1;
savefirstN = 1;
smallNfrm = 500;
CompressData = 0;

% Prompt user for a video file to analyze and thruster type if one is not provided
if nargin==0
    close all
    % Video file prompt
    [filename, pathname] = uigetfile( ...
        {'*.mat','MATLAB MAT-files (*.mat)'; ...
        '*.mraw','Photron MRAW files (*.mraw)'; ...
        '*.tif','Multipage TIFF files (*.tif)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file');
    fn = [pathname filename]
    
    % Thruster type prompt
    set(0,'DefaultFigureWindowStyle','normal')
    Thrusters = {'H6';'X2SingleInner';'X2SingleOuter';'X2DualInner';'X2DualOuter';'X2Dual';'Default75-125';'ManualFit'};
    ThrusterLabels = {'H6';'X2 Single Channel - Inner';'X2 Single Channel - Outer';'X2 Dual Channel - Inner';'X2 Dual Channel - Outer';'X2 Dual - Analyze Both';'Default 75%/125%';'Manual Fit'};
    choice = menu('Discharge Channel Cropping: Choose Pre-set thruster, default 75%-125% or manual fit',ThrusterLabels);
    Thruster = Thrusters{choice}
    % A brief pause allows the user menu to clear
    pause(0.1)
    
end

%% Import Video
tic
% Video files may be either MRAW or MAT files; treat accordingly
if strcmp(fn(end-4:end),'.mraw')
    [ImageInfo, ImageMatrix] = ImportMrawVideo(fn,saveall,savefirstN,smallNfrm,CompressData);
elseif strcmp(fn(end-3:end),'.tif')
        [ImageInfo, ImageMatrix] = ImportTiffVideo(fn,saveall,savefirstN,smallNfrm,CompressData);
elseif strcmp(fn(end-3:end),'.mat')
    CreateMRAW = 0;
    % Check that MAT file was created recently enough to have 'FileInfo'
    % saved, because that specifies a compression version.  If it doesn't
    % have it, might as well re-import from the MRAW to create it.
    load(fn,'FileInfo');
    
    % Check File Compression Version
    FileInfoExists = ~isempty(whos('FileInfo'));
    if FileInfoExists
        if strcmp(FileInfo.CompressionVersion,'1.1') || strcmp(FileInfo.CompressionVersion,'none')
            GoodCompressionVersion = 1;
        else
            GoodCompressionVersion = 0;
        end
    end
    
    % If FileInfo exists AND CompressionVersion is v1.1 or uncompressed,
    % import MAT file.  Otherwise, delete MAT file and try to import from
    % MRAW directly.
    if FileInfoExists && GoodCompressionVersion
        load(fn,'ImageInfo')
        
        if strcmp(FileInfo.CompressionVersion,'1.1')
            
            disp(['Importing compressed MAT file with ' num2str(ImageInfo.NumFrames) ' images...'])
            
            [ImageInfo, ImageMatrix] = ImportMatVideo(fn,CreateMRAW);
        else
            % If the file is uncompressed, it can be loaded directly into memory with the
            % 'load' command
            disp(['Importing uncompressed MAT file with ' num2str(ImageInfo.NumFrames) ' images...'])
            load(fn)
        end
        % In case of loading a small image matrix subset, the ImageMatrix
        % may be labeled as ImageMatrixSmall.  No prob, just re-name.
        if ~isempty(whos('ImageMatrixSmall'));  ImageMatrix = ImageMatrixSmall;  end
    else
        disp('Unrecognized or obsolete file compression version.  Checking for MRAW file...')
        % Create probable MRAW file name
        fn_MRAW = RemoveFileExtension(fn);
        fn_MRAW = [fn_MRAW '.mraw'];
        % Check that file exists; if so import MRAW then delete obsolete
        % MAT file
        MrawExists = ~isempty(dir(fn_MRAW));
        if MrawExists
            disp('MRAW file found, importing...')
            try
                [ImageInfo, ImageMatrix] = ImportMrawVideo(fn_MRAW,saveall,savefirstN,smallNfrm,CompressData);
                disp('Import successful, obsolete MAT file overwritten.')
            catch
                disp('Import not successful.  Manually check MRAW file.')
                return
            end
        end
    end
    
else
    disp('Unknown file type, check that it is MRAW or MAT format')
    return
end
toc

%% Compute Spoke Surface from Video

% In the case of the X2, allow analysis of both channels without havin to
% re-load the image matrix from the hard drive
if strcmp(Thruster,'X2Dual')
    Thruster = 'X2DualInner';
    [spokeSurface spokeSurfFilename] = ComputeSpokeSurface(fn,ImageMatrix,ImageInfo,Thruster);
    spokeSurfFilename
    Thruster = 'X2DualOuter';
    [spokeSurface spokeSurfFilename] = ComputeSpokeSurface(fn,ImageMatrix,ImageInfo,Thruster);
    spokeSurfFilename
else
    Thruster
    [spokeSurface spokeSurfFilename] = ComputeSpokeSurface(fn,ImageMatrix,ImageInfo,Thruster);
    spokeSurfFilename
end
% return
%% Perform 2D DFTs on Spoke Surface and Plot

disp('Rotating spokes propagate in the ExB direction; look for spokes in both cw and ccw directions:')
disp('Computing 2D DFT of spoke surface in forward (cw) direction.')
Video2DDFT(spokeSurfFilename,0)
pause(0.1)

disp('Computing 2D DFT of spoke surface in reverse (ccw) direction')
Video2DDFT(spokeSurfFilename,1)
pause(0.1)
