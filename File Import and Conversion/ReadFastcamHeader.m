function ImageInfo = ReadFastcamHeader(filename)
%% ImageInfo = ReadFastcamHeader(filename)
%       This file takes in a *.cih header file created by Photron's PFV
%       software in conjunction with a saved video file and pulls out
%       relevant video information in the ImageInfo struct array.  For
%       example, an output of this function might look like:
%
% ImageInfo = 
%             Date: '2008/3/9'                      [ Date video was taken ]
%       CameraType: 'FASTCAM-1024PCI model 100K'    [ Photron camera model ]
%              Fps: 109500          [ Framerate in frames per second (fps) ]
%     ShutterSpeed: 9.1324e-006     [ Shutter speed, usually 1 / fps ]
%        NumFrames: 43695           [ Number of frames or images in video ]           
%       ImageWidth: 128             [ Width of image in pixels ]
%      ImageHeight: 16              [ Height of image in pixels ]         
%         BitDepth: 10              [ Bit depth, i.e., # of bits per pixel ]
%         Duration: 0.3990          [ Length of video in seconds ]
%
%       This information is typically used to convert the Photron .mraw
%       file into a MATLAB *.mat file for faster file import and analysis.

%% Read in the header file.
% The *.cih file is just a plain ASCII text file (the .cih extension stands
% for Camera Information Header file). The header data is saved into a cell
% format by the 'textscan' function, with each line of the text file
% turning into a row in the cell array.
fid = fopen(filename);
HeaderData = textscan(fid, '%s %s',...
                        'delimiter', ':', ...
                        'commentStyle', '#',...
                        'CollectOutput', 1);
fclose(fid);

%% Pull out the desired information
% Declare a struct array to fill in the desired information.  The loop
% goes through every row of the cell array looking for an exact match to
% the hard-coded string observed to label that quantity in the *.cih file.
% When it finds a match, it saves the corresponding qantity into a new
% field in the struct array.
%
% Note: this is not very robust, and if Photron changes their output
% routines at all this script may no longer work.
ImageInfo = struct();
for i=1:length(HeaderData{1,1})
    % Frames per second
    if strcmp(HeaderData{1,1}(i,1),'Record Rate(fps) '); ImageInfo.Fps = str2double(char(HeaderData{1,1}(i,2)));           end
    % Image Width
    if strcmp(HeaderData{1,1}(i,1),'Image Width ');      ImageInfo.ImageWidth = str2double(char(HeaderData{1,1}(i,2)));    end
    % Image Height
    if strcmp(HeaderData{1,1}(i,1),'Image Height ');     ImageInfo.ImageHeight = str2double(char(HeaderData{1,1}(i,2)));   end
    % Date
    if strcmp(HeaderData{1,1}(i,1),'Date ');      ImageInfo.Date = char(HeaderData{1,1}(i,2));    end
    % Camera Type
    if strcmp(HeaderData{1,1}(i,1),'Camera Type ');      ImageInfo.CameraType = char(HeaderData{1,1}(i,2));             end    
    % Total Frames
    if strcmp(HeaderData{1,1}(i,1), 'Total Frame ');     ImageInfo.NumFrames = str2double(char(HeaderData{1,1}(i,2)));     end            
    % BitDepth
    if strcmp(HeaderData{1,1}(i,1),'EffectiveBit Depth ');      ImageInfo.BitDepth = str2double(char(HeaderData{1,1}(i,2)));    end    
    % ShutterSpeed
    % -since this is written in the header file literally as '1/100000' or
    % similar, the 'eval' function is necessary to convert this to a number
    if strcmp(HeaderData{1,1}(i,1), 'Shutter Speed(s) ');      ImageInfo.ShutterSpeed = eval(char(HeaderData{1,1}(i,2)));    end                
end
% Duration
ImageInfo.Duration = ImageInfo.NumFrames/ImageInfo.Fps;

%% Error checking
% As a precaution against missing some of this information if the string
% label for a field is updated in future versions of Photron's PFV
% software, check that all of the fields above have been created and filled
% in, and if not return a warning message to the user.

DataFound = length(fieldnames(ImageInfo));

if DataFound<9
    disp(['Possible Error: Information not found in following header file:' sprintf('\r\t') filename])
end


