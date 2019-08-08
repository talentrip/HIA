function [ output_args ] = ComputeSpokeSurfaceBatch( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



%% Initialization and Constants
close all; clc

if nargin==0
    fdir = 'G:\FASTCAM files\10-08-09 H6\600V 10A';
    Thruster = 'H6';
end


% Add a backslash if necessary
if ~strcmp(fdir(end),'\');  fdir = [fdir '\'];  end

% Get filenames in given directory corresponding to PFV header files (.cih)
% [files,err] = get_filenames(fdir,'cih');
files = subdir(fullfile(fdir,'*.mat'));
disp([num2str(length(files)) ' files found...'])

tic
for i=1:numel(files)
    fn = RemoveFileExtension(files(i).name);
    
    disp(['File ' num2str(i) ': ' fn])
    
    % Check that this is an image matrix file (and not, for example, a
    % spoke surface file, which also has a .mat extension)
    load(files(i).name,'FileInfo');
    FileInfoExists = ~isempty(whos('FileInfo'));

    
    % Check that this is not one of the small subset files saved off a
    % larger image matrix.  Better to compute the full spoke surface from
    % the full image matrix.
    IsSmallFile = strcmp(fn(end-5:end),'_small');
    
    if FileInfoExists && ~IsSmallFile
        if strcmp(FileInfo.CompressionVersion,'1.1') || strcmp(FileInfo.CompressionVersion,'none')
            GoodCompressionVersion = 1;
        else
            GoodCompressionVersion = 0;
        end
               
        % If CompressionVersion is v1.1 or uncompressed, import MAT file.
        % Otherwise, delete MAT file and try to import from MRAW directly.
        if GoodCompressionVersion
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
        else
            disp('Unrecognized or obsolete file compression version.  Checking for MRAW file...')
            % Create probable MRAW file name
            fn_MRAW = [fn '.mraw'];
            % Check that file exists; if so import MRAW then delete obsolete
            % MAT file
            MrawExists = ~isempty(dir(fn_MRAW));
            if MrawExists
                disp('MRAW file found, importing...')
                try
                    [ImageInfo, ImageMatrix] = ImportMrawVideo(fn_MRAW,saveall,savefirstN,smallNfrm,CompressData);
                catch
                    disp('Import not successful.  Manually check MRAW file.')
                    return
                end
            end
        end        
        [spokeSurface spokeSurfFilename] = ComputeSpokeSurface(fn,ImageMatrix,ImageInfo,Thruster);
        spokeSurfFilename
            clear FileInfo
    else
        % skip file
        disp('Small file; skipping...')
    end
    
end

