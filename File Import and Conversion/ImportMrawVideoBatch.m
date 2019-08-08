function ImportMrawVideoBatch(fdir)
%% Comments
% This code will recursively find all MRAW files in the given directory and
% any subdirectories, check for any existing .MAT file conversions, and if
% none are found will import and convert the MRAW to MAT format.  If a MAT
% file is found, the code will check if it is up to date and if not, will
% delete the MRAW file and create a new up to date one.

%% Initialization and Constants
close all; clc

if nargin==0
    fdir = 'G:\10-08-09 H6\300V 10A';    
end

saveall = 1;
savefirstN = 1;
smallNfrm = 500;
CompressData = 0;

%% Loop through filenames

% Add a backslash if necessary
if ~strcmp(fdir(end),'\');  fdir = [fdir '\'];  end

% Get filenames in given directory corresponding to PFV header files (.cih)
% [files,err] = get_filenames(fdir,'cih');
files = subdir(fullfile(fdir,'*.mraw'));
disp([num2str(length(files)) ' files found...'])

tic
for i=1:numel(files)
        
    fn = RemoveFileExtension(files(i).name);   
    mraw_filename = [fn '.mraw'];
    mat_filename = [fn '.mat'];
    
    disp(['File ' num2str(i) ': ' fn])
    if isempty(dir(mat_filename))
        disp('No .mat file found.  Importing MRAW file...')
        ImportMrawVideo(mraw_filename,saveall,savefirstN,smallNfrm,CompressData);
    else
        disp('.mat file for this file found, checking...')
        
        % Check that MAT file was created recently enough to have 'FileInfo'
        % saved, because that specifies a compression version.  If it doesn't
        % have it, might as well re-import from the MRAW to create it.
        load(mat_filename,'FileInfo');
        FileInfoExists = ~isempty(whos('FileInfo'));
        
        if FileInfoExists
            disp('.mat file checks out; moving on to next file')
        else
            disp('Older .mat file, lacks FileInfo from .cih header file.')
            disp('Now importing MRAW file...')
            try
                ImportMrawVideo(mraw_filename,saveall,savefirstN,smallNfrm,CompressData);
                disp('Import successful, moving obsolete MAT file to Recycle Bin...')
            catch
                disp('Import not successful.  Manually check MRAW file.')
                return
            end
        end
    end
    
    % To allow text outputs to update
    pause(1)
end


toc




