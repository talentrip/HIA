function SpokeToolboxGUIv6
clc
% Bug Reports:
%
% - Trying to replay simulation after playing once fails (no playback)
%
% Feature To-Do List:
%
% - Show superimposed m=0,1,2,3 1D FFTs in lower right corner
% - Make analysis steps happen on the fly instead of for whole video at
% once, as much as possible

%% Load GUI Window and Assign Callback Functions

% Note: the '@' indicates a so-called "nested function" that is explicitly
% located later in this m-file (see "Callbacks" code section)

% clc
homedir = userpath%'C:\SVN\MATLAB';
homedir = homedir(1:end-1)
GUIdir = fullfile(homedir, 'High Speed Image Analysis\GUI');
GUIFigFile = fullfile(GUIdir,'Spoke_GUI_Front_Panel_v6.fig');
SettingsFile = fullfile(GUIdir,'GUIsettings.mat');
fig = openfig(GUIFigFile);
% To reveal the tags of all the objects in the GUI figure window,
% un-suppress the output of the handles = guihandles(fig) command.
handles = guihandles(fig);

set(handles.AbortButton, 'CallBack', @AbortCallback);
set(handles.StopButton, 'CallBack', @StopCallback,'UserData',1);
set(handles.PlayButton, 'CallBack', @PlayCallback);
set(handles.File, 'CallBack', @File_Callback);
set(handles.PlaybackSlider, 'CallBack', @SliderCallback);
set(handles.Experimental, 'CallBack', @ExperimentalCallback);
set(handles.Simulation, 'CallBack', @SimulationCallback);
set(handles.RateFreqDisplay, 'CallBack', @RateFreqDisplayCallback);
set(handles.RateFreqSlider, 'CallBack', @RateFreqSliderCallback);
set(handles.DurationDisplay, 'CallBack', @DurationDisplayCallback);
set(handles.DurationSlider, 'CallBack', @DurationSliderCallback);
set(handles.AnalysisDisplay1, 'CallBack', @AnalysisDisplay1Callback); % Analysis Frequency Display Callback for surf plot
set(handles.AnalysisSlider1, 'CallBack', @AnalysisSlider1Callback); % Analysis Frequency Callback for surf plot
% set(handles.AnalysisDisplay2, 'CallBack', @AnalysisDisplay2Callback); % Analysis Frequency Display Callback for FFT
% set(handles.AnalysisSlider2, 'CallBack', @AnalysisSlider2Callback); % Analysis Frequency Callback for FFT
set(handles.Close, 'CallBack', @CloseCallback);
set(handles.OpenVideo, 'CallBack', @OpenVideoCallback);
set(handles.Help, 'CallBack', @HelpCallback);
set(handles.RecentFileMenu, 'Callback', @RecentFileMenuCallback);
set(handles.OpenFileButton, 'Callback', @OpenVideoCallback);
set(handles.NormalizeVideoCheckbox, 'Callback', @NormalizeVideoCallback)
set(handles.ACCoupleCheckbox, 'Callback', @ACCoupleCallback)
set(handles.IsolateChannelCheckbox, 'Callback', @IsolateChannelCallback)
set(handles.CalibrateImageCheckbox, 'Callback', @CalibrateImageCallback) % calibrate image checkbox
set(handles.NormFluctCheckbox, 'Callback', @NormFluctCallback) % Log color scaling checkbox
set(handles.CircleFitButton, 'Callback', @CircleFitToolCallback)
set(handles.ResetButton, 'Callback', @ResetButtonCallback);
set(handles.AzimuthSlider, 'Callback', @ViewAngleSliderCallback);
set(handles.ElevationSlider, 'Callback', @ViewAngleSliderCallback);
% set(handles.AzimuthSlider1,'Callback', @ViewAngleSliderCallback1);  % Callback for Azimuth Slider of FFT plot
% set(handles.ElevationSlider1,'Callback',@ViewAngleSliderCallback1);  % Callback for Elevation Slider of FFT plot
set(handles.ColormapButtonGroup, 'SelectionChangeFcn', @ColormapChangeCallback);

%% Assign Necessary Default Values to Objects

set(fig,'WindowStyle','docked',...
    'Renderer', 'OpenGL',...
    'RendererMode','manual');
set(0,'CurrentFigure',handles.BaseGUIFigure)


FrameRate = 1;
FASTCAMFrameRate = 65100;
Duration = 10;
analysisfreq1 = 25;
FFT2Dfreq = 100;
Nchop = 1;      % Nchop necessary for FFT surf plot

set(handles.RateFreqDisplay,'String',num2str(FrameRate))
set(handles.RateFreqSlider,'Value',FrameRate)
set(handles.DurationDisplay,'String',num2str(Duration))
set(handles.DurationSlider,'Value',Duration)
set(handles.AnalysisDisplay1, 'String', num2str(analysisfreq1)) % Set value of analysis frequency display for surf plot
set(handles.AnalysisSlider1, 'Value', analysisfreq1) % Set value of analysis frequency slider for surf plot
clc% set(handles.AnalysisDisplay2, 'String', num2str(FFT2Dfreq)) % Set value of analysis frequency display for FFT
% set(handles.AnalysisSlider2, 'Value', FFT2Dfreq) % Set value of analysis frequency slider for FFT
AzViewAngle = 330;
ElViewAngle = 60;
AzViewAngle1 = 330; % Initial azimuth angle for FFT plot
ElViewAngle1 = 60;  % Initial elevation angle for FFT plot
set(handles.AzimuthSlider,'Value',AzViewAngle)
set(handles.ElevationSlider,'Value',ElViewAngle)
% set(handles.AzimuthSlider1,'Value',AzViewAngle1)     % Azimuth slider for FFT plot
% set(handles.ElevationSlider1,'Value',ElViewAngle1)   % Elevation slider for FFT plot
ColormapSelection = handles.ColorButton;
set(handles.CircleFitButton,'Enable','off')
set(handles.NormalizeVideoCheckbox,'Enable','off')
set(handles.ACCoupleCheckbox,'Enable','off')
set(handles.IsolateChannelCheckbox,'Enable','off')
set(handles.CalibrateImageCheckbox,'Enable','off')  % Calibrate Image Checkbox
set(handles.NormFluctCheckbox,'Enable','off')
set(handles.Experimental, 'Value', 1);
set(handles.PlayButton,'Enable','off'); %disable Play control
set(handles.StopButton,'enable','off');
set(handles.Simulation,'Value',0);
set(handles.RBKButton,'enable','off')
set(handles.RBWButton,'enable','off')
set(handles.ColormapButtonGroup,'SelectedObject',handles.ColorButton);
DisableSimControls;

%% Initialize global variables

NumAzBins = 120;
AzBinAngles = [];
SSPlotSwitch = 0;       % Switch for updating the spoke surface (the contour plot)
MFPlotSwitch = 0;       % Switch for updating movie frames
LFPlotSwitch = 0;       % Switch for updating line frames
FFTPlotSwitch = 0;      % Switch for updating the FFT surface plot
ColorBitDepthSwitch = 0;
PlaybackMode = 'Experimental';

% Switches are created for the Spoke Surface and FastCam Viewer in the hope
% that, while the first call to display the image will necessarily be slow,
% since the figure, axes, view, etc. must all be updated, that a low-level
% function may be found that will allow the image to be updated rapidly
% after the first time.  For now, all attempts at this have failed.

CircleCenter = [];
CircleRadius = [];
RFilter = [];
Mask2D = [];
Mask3D = [];
handles.FileName = [];
ImageMatrix = [];
MeanImage = [];
MaxImage = [];
MinImage = [];
PixelRange = [];
% maxval = 0;
% minval = 0;
FrameIntensity = [];
MaxPixel = 0;
MinPixel = 0;
NewMaxPixel = 0;
NewMinPixel = 0;
NumFrames = 0;
Normalization = [];
cmavgline = [];
SliderFraction = 0;
IsolateChannel = get(handles.IsolateChannelCheckbox,'Value');
ACCouple = get(handles.ACCoupleCheckbox,'Value');
NormalizeVideo= get(handles.NormalizeVideoCheckbox,'Value');
CalibrateImage = get(handles.CalibrateImageCheckbox,'Value'); % Calibrate Image Checkbox
NormFluct = get(handles.NormFluctCheckbox,'Value');

SpokeSurfaceHandle=[];
FFTSurfaceHandle=[];
FFTZLims = [];
MovieFrameHandle=[];
LineFrameHandle=[];
cores = feature('Numcores');

CurrentFrame = 0;
DisplayColorBitDepth = 6;
MaxRecentFiles = 5;

%% Main Function
% If any operations were specified to be carried out by default upon
% opening the program, they would be called here.  In the debug, in lieu of
% selecting the file to be opened with the File --> Open menu, it is just
% loaded automatically (useful for timing purposes).

% Load any saved settings
forceload = 0;
if(isempty(dir(SettingsFile)) || forceload)
    defaultmsg = {'No Recently Viewed Files'};
    save(SettingsFile, 'defaultmsg')
    set(handles.RecentFileMenu,'String',defaultmsg)
else
    load(SettingsFile,'recentFiles');
    set(handles.RecentFileMenu, 'String', recentFiles)
end

debug = 0;
if debug
    set(handles.BaseGUIFigure,'HandleVisibility','on')
    get(0,'CurrentFigure')
    %     handles.FileName = 'C:\Users\Mike\Documents\PEPL\Gerbil Videos and Analysis\rawimg_small.mat';
    %     CircleCenter = [59.2366   58.8150];
    %     CircleRadius = 44.0117;
    %     CircleIR = 35.9816;
    %     CircleOR = 54.2673;
    
    handles.FileName = 'C:\Users\Mike\Desktop\Video Cases\600V 10A s.mat';
    CircleCenter = [66.0845 60.1776];
    CircleRadius = 48.0763;
    CircleIR = 35.8767;
    CircleOR = 58.9167;
    
    OpenVideoCallback(0,0,handles.FileName);
    
    % Copied from CircleFit code
    
    
    r = PixelPolarCoords(MeanImage,CircleCenter);
    RFilter = r>CircleIR & r<CircleOR;
    Mask2D = ~RFilter;
    MeanImage(Mask2D)=0;
    pcolor(handles.Axis5,MeanImage);
    shading(handles.Axis5,'flat')
    daspect(handles.Axis5,[1 1 1]);
    set(handles.PlayButton,'Enable','on');
    get(0,'CurrentFigure')
    
    % NumFrames = 7;
    PlayCallback(0,0);
    AbortCallback(0,0);
    
    
end

%% Callbacks
% This section should contain the callback functions of all buttons on the
% front panel.

% Playback Controls
    function PlayCallback(h, dummy)
        ResetSwitches;
        set(handles.PlayButton,'enable','off');
        set(handles.StopButton,'enable','on');
        
        if strcmp(PlaybackMode,'Simulation')
            spokeFreq = FrameRate;
            defaultFPS = 10;
            NumFrames = Duration*defaultFPS;
            imagePixels = 120;
            CircleCenter = imagePixels/2*[1 1];
            channelPixelRad = 45;
            channelPixelWidth = 15;
            spokeWidthDeg = 30;
            ImageMatrix = SpokeSimulator(imagePixels,channelPixelRad,channelPixelWidth,spokeFreq,spokeWidthDeg,NumFrames,0.2);
            size(ImageMatrix)
            MeanImage = mean(ImageMatrix,3);
            RFilter = true(size(MeanImage));
            % need rest of circle fitting for the simulated
            % spokes so they can be plotted and analyzed too
        end
        
        % Bin in theta
        [r,theta,azBinPixels,AzBinAngles,PixelsPerAzBin] = do_thetabinning(MeanImage,CircleCenter,RFilter,NumAzBins);
        
        %         switch PlaybackMode
        %             case handles.ExperimentalButton
        %             case handles.SimulationButton
        %         end
        
        
        %% Playback Loop
        t2=0; t1=0; tic
        
        % Initialize cmavg matrices
        cmavg1 = zeros(analysisfreq1,NumAzBins);
        cmavg2 = zeros(FFT2Dfreq,NumAzBins);
        
        while (CurrentFrame < NumFrames && get(handles.StopButton,'UserData'))
            tic

            MinPixel = NewMinPixel;
            MaxPixel = NewMaxPixel;
            % Update Frame Counter
            CurrentFrame = CurrentFrame + 1;
            set(handles.PlaybackSlider,'Value',CurrentFrame/NumFrames); %updates sliders position --> slider will change accordingly
            set(handles.FrameCounter,'String',[num2str(CurrentFrame) ' / ' num2str(NumFrames)]);
            
            % Update Current Image and filter according to checkboxes
            currentImage = ImageMatrix(:,:,CurrentFrame);
            cmavgline = thetafn(azBinPixels,NumAzBins,PixelsPerAzBin,currentImage);
            
            if IsolateChannel == 1;
                currentImage(Mask2D) = 0;
            end
            
            % Update Video Frame
            if 1  %strcmp(PlaybackMode,'Experimental')
                MovieFramePlot(currentImage);
                LineFramePlot(AzBinAngles,cmavgline)
            else
                pcolor(handles.Axis5,currentImage);%displays the gerbil at location thetapicker(:,:,CurrentFrame)
            end
            
            %% Periodic Analysis modules (toolkits) happen here
            %--------------------------------------------------------------------------
            % Generate and Update Spoke Surf Plot and FFT Surf Plot
            
            if CurrentFrame <= analysisfreq1
                [cmavg1(CurrentFrame,:)] = thetafn(azBinPixels,NumAzBins,PixelsPerAzBin,currentImage);
            elseif CurrentFrame > analysisfreq1
                [cmavg1((analysisfreq1 + 1),:)] = thetafn(azBinPixels,NumAzBins,PixelsPerAzBin,currentImage);
                cmavg1 = cmavg1((2:analysisfreq1 + 1),:);
            end
            
            SpokeSurfPlot(cmavg1); % Spoke Surf Plot Function
            
            if CurrentFrame <= FFT2Dfreq
                [cmavg2(CurrentFrame,:)] = thetafn(azBinPixels,NumAzBins,PixelsPerAzBin,currentImage);
            elseif CurrentFrame > FFT2Dfreq
                [cmavg2((FFT2Dfreq + 1),:)] = thetafn(azBinPixels,NumAzBins,PixelsPerAzBin,currentImage);
                cmavg2 = cmavg2((2:FFT2Dfreq + 1),:);
            end
            
%             if CurrentFrame >= FFT2Dfreq
%                 FFTSurfPlot(cmavg2,Nchop); % FFT Surf Plot Function
%             end
            
            % Pause briefly to maintain proper framerate
            if toc<1/FrameRate
              pause( 1/FrameRate - toc )  
            end
        end
        ResetSwitches;
        set(handles.StopButton,'UserData',1);
        set(handles.PlayButton,'enable','on');
        set(handles.StopButton,'enable','off');
        if CurrentFrame==NumFrames; set(handles.PlaybackSlider,'Value',0); end
    end

    function StopCallback(h, dummy)
        set(handles.StopButton,'UserData',0);
        set(handles.PlayButton,'enable','on');
        set(handles.StopButton,'enable','off');
    end

    function SliderCallback(h, dummy)
        %Video position fraction is defined by the slider position
        SliderFraction = get(gcbo,'Value')/((get(gcbo,'Max'))-(get(gcbo,'Min')));
        CurrentFrame = SliderFraction*NumFrames;
        CurrentFrame = round(CurrentFrame);
    end

% Spoke Simulation Controls
    function SimulationCallback(h, dummy)
        PlaybackMode = 'Simulation';  %Simulation radio button has been pressed update handles.playback_mode accordingly
        EnableSimControls;
        set(handles.RateFreqSliderLabel,'String','Simulated Spoke Angular Frequency');
        set(handles.PlayButton,'Enable','on');
    end

    function RateFreqDisplayCallback(h, dummy)
        val = str2double(get(handles.RateFreqDisplay,'string'));     %gets the number displayed on the screen
        upLim = 20;
        if (~isnumeric(val) || val < 1 || val > upLim)%the following statement is called if the user types in a number greater than 18, less than 1 or doesn't type in a number
            if ~isnumeric(val); val=upLim/2; end
            if val<1; val=1; end
            if val>upLim; val=upLim; end
        end
        FrameRate = val;
        set(handles.RateFreqSlider,'Value',FrameRate);
    end

    function RateFreqSliderCallback(h, dummy)
        val = get(handles.RateFreqSlider,'Value');    %returns position of slider as a Value
        val = ceil(val);     %rounds value of slider to integer toward infinity
        set(handles.RateFreqDisplay,'String',val)
        FrameRate = val;
    end

    function DurationDisplayCallback(h, dummy)
        Duration = str2double(get(handles.DurationDisplay,'string'));    %gets the number displayed on the edit text screen
        upLim = 60;
        if (~isnumeric(Duration) || Duration < 1 || Duration > 50) %called if the user types in a number greater than 50, less than 1 or doesn't type in a number
            if ~isnumeric(val); val=upLim/2; end
            if val<1; val=1; end
            if val>upLim; val=upLim; end
        end
        Duration = val;
        set(handles.DurationSlider,'Value',Duration);
    end

    function DurationSliderCallback(h, dummy)
        slider_value = get(handles.DurationSlider,'Value');    %returns position of slider as a Value
        slider_value = round(slider_value);     %rounds value of slider to integer
        set(handles.DurationDisplay,'String',slider_value);  %updates the value obtained by the slider to the corresponding edit text box
        handles.duration = slider_value;    %updates the handles.duration variable to be used in other subroutines
    end

    function ColormapChangeCallback(h,button)
        ColormapSelection = button.NewValue;
    end

% Spoke Experimental Controls

    function ExperimentalCallback(h, dummy)
        PlaybackMode = 'Experimental';  %Experimental radio button has been pressed update handles.playback_mode accordingly
        DisableSimControls;
        set(handles.RateFreqSliderLabel,'String','Playback Framerate (fps)');
        if ~isempty(handles.FileName)
            handles.FileName
            OpenVideoCallback(0,0,handles.FileName);
        end
    end

    function OpenVideoCallback(h, dummy, givenfile)
        forceExit = 0;
        if nargin==3
            handles.FileName = givenfile;
        else
            % Prompt user to select image matrix (*.mat) file to open
            [FileName,PathName] = uigetfile('*.mat','Select a .mat file');
            handles.FileName = fullfile(PathName,FileName);
            if isempty(dir(handles.FileName))
                forceExit = 1;
            end
        end
        
        if forceExit; return; end
        % Load file into memory
        clear ImageMatrix
        ImageMatrix = load(handles.FileName);

        if isfield(ImageMatrix,'MeanImage');  MeanImage = ImageMatrix.MeanImage;  end
        if isfield(ImageMatrix,'ImageMatrixSmall');  ImageMatrix = ImageMatrix.ImageMatrixSmall;  end
        if isfield(ImageMatrix,'img');  ImageMatrix = ImageMatrix.img;  end
        ImageMatrix = double(ImageMatrix);
        NumFrames = size(ImageMatrix,3);
        CurrentFrame = 0;
        
        % CheckFigureProperties;
        ResetMinMaxPixels;
% whos
        PixelRange = MaxImage - MinImage;
        
        
        
        RFilter = logical(MeanImage);
        
        FrameIntensity = sum(sum(ImageMatrix));
        %         maxval = max(max(ImageMatrix(:,:,1)));
        %         minval = min(min(ImageMatrix(:,:,1)));
        Normalization = mean(FrameIntensity)/FrameIntensity;
        
        if ColorBitDepthSwitch
            ImageMatrix = (ImageMatrix-MinPixel)/(MaxPixel-MinPixel) * 2^DisplayColorBitDepth;
        end
        
        % User interface adjustment
        set(handles.CircleFitButton,'Enable','on')
        set(handles.FrameCounter,'String',['0 / ' num2str(NumFrames)]);
        CurrentFrame = 0;
        set(handles.PlaybackSlider,'Value',0)
        set(handles.NormalizeVideoCheckbox,'Enable','on','Value',0)
        set(handles.ACCoupleCheckbox,'Enable','on','Value',0)
        set(handles.IsolateChannelCheckbox,'Enable','on','Value',0)
        set(handles.CalibrateImageCheckbox,'Enable','on','Value',0) % Calibrate Image Checkbox
                set(handles.NormFluctCheckbox,'Enable','on','Value',0)

        % Set Limits Of Analysis Frequency Slider and Display
        set(handles.AnalysisDisplay1,'Min',25)
        set(handles.AnalysisDisplay1,'Max',NumFrames)
        set(handles.AnalysisSlider1,'Min',25)
        set(handles.AnalysisSlider1,'Max',NumFrames)
%         set(handles.AnalysisDisplay2,'Min',10)
%         set(handles.AnalysisDisplay2,'Max',NumFrames)
%         set(handles.AnalysisSlider2,'Min',10)
%         set(handles.AnalysisSlider2,'Max',NumFrames)
        
        
        % Update recent file list
        recentFiles = get(handles.RecentFileMenu, 'String');
        if strcmp(recentFiles,'No Recently Viewed Files')
            clear recentFiles
            recentFiles{1} = 'Recently Viewed Files (Select to Open)...';
        end
        alreadyOnList = cellfun(@(x) strcmp(x,handles.FileName),recentFiles);
        recentFiles(alreadyOnList) = [];
        recentFiles(3:end+1) = recentFiles(2:end);
        currentFile = 2;
        recentFiles(currentFile) = {handles.FileName};
        recentFiles(MaxRecentFiles+1:end)=[];
        set(handles.RecentFileMenu, 'String', recentFiles);
        set(handles.RecentFileMenu, 'Value', currentFile);
        
        % Show rough fit in GUI window
        [CircleCenter(1) CircleCenter(2) CircleRadius] = TaubinSVDThresholdFit(MeanImage);
        pcolor(handles.Axis5,MeanImage);
        shading(handles.Axis5,'flat')
        daspect(handles.Axis5,[1 1 1])
    end

% Function for Surf Plot Analysis Frequency Display
    function AnalysisDisplay1Callback(h, dummy)
        
        val = str2double(get(handles.AnalysisDisplay1,'string'));     %gets the number displayed on the screen
        upLim = NumFrames;
        if (~isnumeric(val) || val < 25 || val > upLim)
            if ~isnumeric(val); val = floor(upLim/2); end
            if val < 25; val = 25; end
            if val>upLim; val = upLim; end
        end
        analysisfreq1 = val;
        set(handles.AnalysisSlider1,'Value',analysisfreq1);
    end

% Function for Surf Plot Analysis Frequency Slider
    function AnalysisSlider1Callback(h, dummy)
        
        val = get(handles.AnalysisSlider1,'Value');    %returns position of slider as a Value
        val = ceil(val);     %rounds value of slider to integer toward infinity
        set(handles.AnalysisDisplay1,'String',val)
        analysisfreq1 = val;
    end

% Function for FFT Analysis Frequency Display
    function AnalysisDisplay2Callback(h, dummy)
        
        val = str2double(get(handles.AnalysisDisplay2,'string'));     %gets the number displayed on the screen
        upLim = NumFrames;
        if (~isnumeric(val) || val < 100 || val > upLim)
            if ~isnumeric(val); val = floor(upLim/2); end
            if val < 100; val = 100; end
            if val>upLim; val = upLim; end
        end
        FFT2Dfreq = val;
        set(handles.AnalysisSlider2,'Value',FFT2Dfreq);
    end

% Function for FFT Analysis Frequency Slider
    function AnalysisSlider2Callback(h, dummy)
        
        val = get(handles.AnalysisSlider2,'Value');    %returns position of slider as a Value
        val = ceil(val);     %rounds value of slider to integer toward infinity
        set(handles.AnalysisDisplay2,'String',val)
        FFT2Dfreq = val;
    end

    function RecentFileMenuCallback(h, dummy)
        % Determine the selected data set.
        str = get(handles.RecentFileMenu, 'String');
        val = get(handles.RecentFileMenu, 'Value');
        % If the user has selected a different file, one that is not the
        % initial 'Recent Files' message, open teh new video
        if val~=1
            OpenVideoCallback(0,0,str{val})
        end
    end

    function ViewAngleSliderCallback(h, dummy)
        AzViewAngle = get(handles.AzimuthSlider, 'Value');
        ElViewAngle = get(handles.ElevationSlider, 'Value');
        view(handles.Axis2,[AzViewAngle,ElViewAngle]);
    end

    function ViewAngleSliderCallback1(h,dummy)
        AzViewAngle1 = get(handles.AzimuthSlider1, 'Value');
        ElViewAngle1 = get(handles.ElevationSlider1, 'Value');
        view(handles.Axis8,[AzViewAngle1,ElViewAngle1]);
    end

    function CircleFitToolCallback(h, dummy)
        if isempty(handles.FileName); OpenVideoCallback; end
        [CircleCenter CircleRadius CircleIR CircleOR] = CircleFitGui(MeanImage,GUIdir);
        r = PixelPolarCoords(MeanImage,CircleCenter);
        RFilter = r>CircleIR & r<CircleOR;
        Mask2D = ~RFilter;
        sampleImage = MeanImage;
        sampleImage(Mask2D) = 0;
        set(handles.BaseGUIFigure,'CurrentAxes',handles.Axis5);
        image(sampleImage,'CDataMapping','scaled');
        shading(handles.Axis5,'flat')
        daspect(handles.Axis5,[1 1 1]);
        set(handles.PlayButton,'Enable','on');
    end

% Calibrate Image
    function CalibrateImageCallback(h,dummy)
        CalibrateImage = get(handles.CalibrateImageCheckbox,'Value');
        
        
        if CalibrateImage == 1  % Calibrate Image Matrix
            
            % Convert Original Pixel Value to Brightness
            ImageMatrix = -2.11 * log((1/-118895.47)*(ImageMatrix-65473));
            % Convert Brightness to Calibrated Pixel Value
            ImageMatrix = 22178*ImageMatrix - 23640;
            CurrentImage = ImageMatrix(:,:,CurrentFrame);
            set(MovieFrameHandle,'CData',CurrentImage)
            
        elseif CalibrateImage == 0  % Uncalibrate Image Matrix
            % Convert Calibrated Pixel Value to Brightness
            ImageMatrix = (1/22178)*(ImageMatrix+23640);
            % Convert Brightness to Original Pixel Value
            ImageMatrix = (65473)-(118895.47*exp((-1/2.11)*ImageMatrix));
            CurrentImage = ImageMatrix(:,:,CurrentFrame);
            set(MovieFrameHandle,'CData',CurrentImage)
            
        end
    end

% Normalize Image
    function NormalizeVideoCallback(h, dummy)
        NormalizeVideo = get(handles.NormalizeVideoCheckbox,'Value');
        
        if NormalizeVideo
            % ImageMatrix = ImageMatrix * Normalization
            ImageMatrix = bsxfun(@times,ImageMatrix,Normalization);
            set(handles.CalibrateImageCheckbox,'enable','off')
        else
            % ImageMatrix = ImageMatrix / Normalization
            ImageMatrix = bsxfun(@times,ImageMatrix,1/Normalization);
            set(handles.CalibrateImageCheckbox,'enable','on')
        end
        CurrentImage = ImageMatrix(:,:,CurrentFrame);
        set(MovieFrameHandle,'CData',CurrentImage)
        
        RefreshGlobalImageParameters;
        
        
    end

% Filter Mean Image
    function ACCoupleCallback(h, dummy)
        ACCouple = get(handles.ACCoupleCheckbox,'Value');
        
        if ACCouple
            % Before AC-Coupling, video must be normalized, otherwise it's just a pain to figure out what you're seeing (overall oscillations dominate azimuthal oscillations)
            if ~NormalizeVideo
                set(handles.NormalizeVideoCheckbox,'Value',1);
                NormalizeVideoCallback;
            end
            % Subtracts MeanImage from all images in ImageMatrix
            ImageMatrix = bsxfun(@minus,ImageMatrix,MeanImage);
            ResetMinMaxPixels;
            % Take care of color choice radio buttons
            set(handles.ColormapButtonGroup,'SelectedObject',handles.RBWButton)
            set(handles.BWButton,'enable','off');   set(handles.RBWButton,'enable','on');   %set(handles.RBKButton,'enable','on');   
            ColormapSelection = handles.RBWButton;
        else
            % Restore to non-normalized values before de-AC Coupling if NormFluct option is set
            if NormFluct
                set(handles.NormFluctCheckbox,'Value',0);
                NormFluctCallback;
            end
            % Restores (adds) MeanImage to all images in ImageMatrix
            ImageMatrix = bsxfun(@plus,ImageMatrix,MeanImage);
            ResetMinMaxPixels;
            set(handles.ColormapButtonGroup,'SelectedObject',handles.ColorButton)
            set(handles.BWButton,'enable','on');    set(handles.RBKButton,'enable','off');  set(handles.RBWButton,'enable','off')
            ColormapSelection = handles.ColorButton;
            
        end
        CurrentImage = ImageMatrix(:,:,CurrentFrame);
        CurrentImage = MakeMeanMap(CurrentImage);
        try
        set(MovieFrameHandle,'CData',CurrentImage)
        catch ME
           min(CurrentImage(:))
           max(CurrentImage(:))
                   CurrentImage = ImageMatrix(:,:,CurrentFrame);
           min(CurrentImage(:))
           max(CurrentImage(:))
                   MinPixel
                   MaxPixel
           rethrow(ME)
        end
    end

    function IsolateChannelCallback(h, dummy)
        % This function applied frame-by-frame within playback loop under
        % an if-statement
        IsolateChannel = get(handles.IsolateChannelCheckbox,'Value');
        ResetMinMaxPixels;
    end

    function NormFluctCallback(h, dummy)
        % This function applies the natural log to the visible part of the
        % image.  It is only intended to be used with the AC-Coupling
        % feature.  The initial idea of it was to better visualize small
        % oscillations in front of the poles without getting washed out by
        % channel oscillations (just as channel oscillations are washed out
        % by the bright cathode).
        NormFluct = get(handles.NormFluctCheckbox,'Value');
        PixelRange = range(ImageMatrix,3);
        
        if NormFluct
            % Only makes sense to normalize fluctuations in AC-coupled mode
            if ~ACCouple
                set(handles.ACCoupleCheckbox,'Value',1);
                ACCoupleCallback;
            end
            
            % Compute normalized fluctuation for each pixel            
            % AC-coupled ImageMatrix = AC-coupled ImageMatrix ./ PixelRange;
            ImageMatrix = bsxfun(@elementdivide,ImageMatrix,PixelRange);
        else
            % AC-coupled ImageMatrix = AC-coupled ImageMatrix .* PixelRange;
            ImageMatrix = bsxfun(@elementmultiply,ImageMatrix,PixelRange);
        end
        ResetMinMaxPixels;
    end
% GUI Window Controls

    function HelpCallback(h, dummy)
        msgbox( ...
            ['To watch a video, select either the simulation or experimental radio button under movie preferences.' ...
            '  Select video preference parameters (angular velocity or duration) by using the slider or text boxes.' ...
            '  Click the Play button to start the video.  The contour plot and other tools will update in real time.'...
            '  The duration parameter is defined for both experimental and simulation mode, where the angular velocity parameter is strictley defined for simulation mode.'...
            '  Stop Button and Playback Slider'...
            '  To rewind the video press the stop button and use the playback slider to go back to the desired part of the video, similarly for fast forward.'...
            '  Loading an experimental video'...
            '  If the experimental video radio button is pressed, you will be asked to locate the desired experimental video you wish to play.'...
            '  The file you select must be a .mat file.  This tool only recognizes videos that have been converted to .mat video files.'], ...
            '  Instructions','help')
    end

    function AbortCallback(h, dummy)
        StopCallback;
        recentFiles = get(handles.RecentFileMenu, 'String');
        save(SettingsFile, 'recentFiles')
        delete(handles.BaseGUIFigure);
        return
    end

    function File_Callback(h, dummy)
    end

    function ResetButtonCallback(h,dummy)
        
        handles.angfreq = 10; %sets default for angular frequency
        handles.duration = 1; %sets default for duration
        handles.exp = 1;
        handles.sim = 2;
        handles.framenumbercounter = 0;
        handles.totalnumberframes = 0;
        
        set(handles.CircleFitButton,'Enable','off')
        set(handles.NormalizeVideoCheckbox,'Enable','off')
        set(handles.ACCoupleCheckbox,'Enable','off')
        set(handles.IsolateChannelCheckbox,'Enable','off')
        set(handles.CalibrateImageCheckbox,'Enable','off') % calibrate image
        set(handles.Experimental, 'Value', 1);
        set(handles.PlayButton,'Enable','off'); %disable Play control
        set(handles.Simulation,'Value',0);
        DisableSimControls;
        
        handles.videopositionfraction = 0; %initialize videopositionfraction - the video position fraction will be defined as the fraction of the framenumber divided by the total number of frames of the current movie
        handles.ControlBarPressed = 0; %Control bar hasn't been pressed yet
        
        NumAzBins = 120;
        SSPlotSwitch = 0;       % Switch for updating the spoke surface (the contour plot)
        MFPlotSwitch = 0;       % Switch for updating movie frames
        LFPlotSwitch = 0;       % Switch for updating line frames
        FFTPlotSwitch = 0;      % Switch for updating the FFT surface plot
        ColorBitDepthSwitch = 0;
        PlaybackMode = 'Experimental';
        
        CircleCenter = [];
        CircleRadius = [];
        RFilter = [];
        Mask2D = [];
        Mask3D = [];
        handles.FileName = [];
        ImageMatrix = [];
        MeanImage = [];
        MaxImage = [];
        MinImage = [];
        FrameIntensity = [];
        MaxPixel = 0;
        MinPixel = 0;
        NumFrames = 0;
        Normalization = [];
        SliderFraction = 0;
        IsolateChannel = 0;
        set(handles.IsolateChannelCheckbox,'Value',0);
        ACCouple = 0;
        set(handles.ACCoupleCheckbox,'Value',0);
        NormalizeVideo = 0;
        set(handles.NormalizeVideoCheckbox,'Value',0);
        CalibrateImage = 0;
        set(handles.CalibrateImageCheckbox,'Value',0); % Calibrate Image Checkbox
        SpokeSurfaceHandle=[];
        MovieFrameHandle=[];
        LineFrameHandle=[];
        cores = feature('Numcores');
        ColormapSelection = handles.ColorButton;
        
        FrameRate = 15;
        Duration = 10;
        analysisfreq1 = 25;
        FFT2Dfreq = 100;
        Nchop = 3;      % Nchop necessary for FFT surf plot
        set(handles.RateFreqDisplay,'String',num2str(FrameRate))
        set(handles.RateFreqSlider,'Value',FrameRate)
        set(handles.DurationDisplay,'String',num2str(Duration))
        set(handles.DurationSlider,'Value',Duration)
        set(handles.AnalysisDisplay1, 'String', num2str(analysisfreq1)); % Reset value of surf plot analysis frequency display
        set(handles.AnalysisSlider1, 'Value', analysisfreq1); % Reset value of surf plot analysis frequency slider
%         set(handles.AnalysisDisplay2, 'String', num2str(FFT2Dfreq)); % Reset value of FFT analysis frequency display
%         set(handles.AnalysisSlider2, 'Value', FFT2Dfreq); % Reset value of FFT analysis frequency slider
        AzViewAngle = 330;
        ElViewAngle = 60;
        AzViewAngle1 = 330; % Initial azimuth angle for FFT plot
        ElViewAngle1 = 60;  % Initial elevation angle for FFT plot
        set(handles.AzimuthSlider,'Value',AzViewAngle)
        set(handles.ElevationSlider,'Value',ElViewAngle)
%         set(handles.AzimuthSlider1,'Value',AzViewAngle1)     % Azimuth slider for FFT plot
%         set(handles.ElevationSlider1,'Value',ElViewAngle1)   % Elevation slider for FFT plot
        CurrentFrame = 0;
        DisplayColorBitDepth = 6;
        MaxRecentFiles = 5;
        
    end

%% Helper Functions
% These are just commonly called sets of commands that enable or disable
% various controls depending on the mode.  To improve code readability,
% they are renamed and collected here.

    function DisableSimControls
        set(handles.DurationDisplay,'Enable','off');%disable Duration Display control
        set(handles.DurationSlider,'Enable','off');%disable Duration slider
        %         set(handles.Simulation,'Enable','off');%disable simulation controla
    end

    function EnableSimControls
        set(handles.DurationDisplay,'Enable','on');%disable Duration Display control
        set(handles.DurationSlider,'Enable','on');%disable Duration slider
    end

    function ResetSwitches
        SSPlotSwitch = 0;
        MFPlotSwitch = 0;
        LFPlotSwitch = 0;
        FFTPlotSwitch = 0;  % ResetSwitch for FFT plot
    end

    function SpokeSurfPlot(data)
        if SSPlotSwitch  % refresh spoke surface plot
            set(SpokeSurfaceHandle,'ZData',data)
        else  % create spoke surface plot for first time
            SpokeSurfaceHandle = surf(handles.Axis2,AzBinAngles,(size(data,1):-1:1)/FASTCAMFrameRate*1000,data);
            view(handles.Axis2, [AzViewAngle,ElViewAngle])
            % set(SpokeSurfaceHandle,'XDataMode','manual','YDataMode','manual','ZDataMode','manual')
            shading(handles.Axis2,'interp')
            SSPlotSwitch = 1;
            
            xlabel(handles.Axis2,'Azimuthal Angle, deg')
            ylabel(handles.Axis2,'Elapsed Time, msec')
            zlabel(handles.Axis2,'Average Pixel Intensity, Arb. Units')
            set(handles.Axis2,'XTick',[-180 -90 0 90 180])
            xlim(handles.Axis2,[-180 180])
            
        end
    end

% Function for FFT plot
%     function FFTSurfPlot(cmavg,Nchop)
%         
%         numthetabins = size(cmavg,2);
%         Nfiles = size(cmavg,1);
%         
%         subset = floor(Nfiles/Nchop);
%         Z=zeros(subset,numthetabins,Nchop);
%         
%         %for i =1:Nchop; Z(:,:,i) = fft2(cmavg((i-1)*subset+1:i*subset,:)); end
%         for i =1:Nchop; Z(:,:,i) = fft2(cmavg((i-1)*subset+1:i*subset,:)); end
%         Z = mean(Z,3);
%         
%         %eliminate redundant second half of FFT matrix
%         Z = Z(1:floor(subset/2),1:floor(numthetabins/2));
%         mrange = 7;
%         Z = abs(real(Z(:,1:mrange)));
%         
%         if FFTPlotSwitch
%             
%             %             set(h
%             % If range of current line plot exceeds axes, adjust axes range
%             zmax1 = max(Z(:));
%             zmin1 = min(Z(:));
%             
%             % If max of data is not equal to
%             if zmax1>FFTZLims(2) || zmax1<0.1*FFTZLims(2)
%                 FFTZLims(2) = zmax1;
%                 zlim(handles.Axis8,FFTZLims)
%             end
%             if zmin1<FFTZLims(1)
%                 FFTZLims(1) = zmin1;
%                 zlim(handles.Axis8,FFTZLims)
%             end
%             set(FFTSurfaceHandle,'ZData',Z)
%             
%         else
%             %             set(handles.BaseGUIFigure,'CurrentAxes',handles.Axis8)
%             %             cla;
%             FFTSurfaceHandle = surf(handles.Axis8,(1:mrange)-1,(1:floor(subset/2))*FASTCAMFrameRate/Nfiles,Z);
%             %             shading interp
%             % Commented out code below labels the axes for the FFT plot
%             xlabel(handles.Axis8,'Mode Number');
%             ylabel(handles.Axis8,'Spoke Local Frequency, Hz');
%             zlabel(handles.Axis8,'Mode Amplitude, Arbitrary Units');
%             xlim(handles.Axis8,[0 mrange-1])
%             ylim(handles.Axis8,[0 FASTCAMFrameRate/2])
%             FFTZLims = [0 1];
%             zlim(handles.Axis8,FFTZLims)
%             view(handles.Axis8,[AzViewAngle1,ElViewAngle1])
%             shading(handles.Axis8,'interp')
%             FFTPlotSwitch = 1;
%         end
%         
%     end

    function MovieFramePlot(data)
        set(handles.BaseGUIFigure,'CurrentAxes',handles.Axis5);
        data = MakeMeanMap(data);
        
        if MFPlotSwitch  % refresh movie frame
            set(MovieFrameHandle,'CData',data)
        else  % create movie frame for first time
            MovieFrameHandle = image(data,'CDataMapping','scaled');
            daspect([1 1 1]);
            % set(hOut,'erasemode','xor');
            % trip switch so next time image is just refreshed
            MFPlotSwitch = 1;
        end
    end

    function LineFramePlot(xdata,ydata)
        
        if LFPlotSwitch
            set(LineFrameHandle,'YData',ydata)
            
            % If range of current line plot exceeds axes, adjust axes range
            ymax1 = max(ydata);
            ymin1 = min(ydata);
            ylims = get(handles.Axis4,'YLim');
            
            % If max of data is not equal to
            if ymax1>ylims(2) || ymax1<0.1*ylims(2)
                ylims(2) = ymax1;
%                 try
%                 set(handles.Axis4,'YLim',[ylims(1) ymax])
%                 catch
%                     ylims
%                     ymax1
%                     ymin1
%                     ymax
%                 end
%                 ylims = get(handles.Axis4,'YLim');
            end
            
            if ymin1<ylims(1) || ymin1>0.1*ylims(1)
                ylims(1) = ymin1;
%                 set(handles.Axis4,'YLim',[ymin ylims(2)])
            end
            % Check for pathological blank screen cases (all zero pixels)
            if range(ylims)==0
                ylims(2) = ylims(2)+1;
            end
            set(handles.Axis4,'YLim',ylims)
        else
            set(handles.BaseGUIFigure,'CurrentAxes',handles.Axis4)
            cla;
            LineFrameHandle = line(xdata,ydata);
            set(LineFrameHandle,'Marker','.');
            xlim([-180 180]);
            ylim([0 5]);
            LFPlotSwitch = 1;
            xlabel('Azimuthal Angle, deg')
            ylabel('Average Pixel Intensity, Arb. Units')
        end
        
    end

    function data = MakeMeanMap(data)
        switch ColormapSelection
            
            case handles.ColorButton
                % do nothing
                
            case handles.RBWButton
                % Normalize to largest pixel value (positive or negative)
                data = data/max(abs([MinPixel MaxPixel]));
                % Create a mask of pixels the same size as 'data' with 1's
                % where the data is positive and zeros otherwise
                posmask = data>0;
                % Create a white background in RGB space
                rMap = ones(size(data));
                gMap = ones(size(data));
                bMap = ones(size(data));
                % Create the red/blue on white colormap using the 
                gMap(posmask) = 1-data(posmask);
                bMap(posmask) = 1-data(posmask);
                gMap(~posmask) = 1+data(~posmask);
                rMap(~posmask) = 1+data(~posmask);
                data = cat(3,rMap,gMap,bMap);
                
            case handles.BWButton
                % change to black and white by making a 3-page matrix of
                % all same values
                data = data/max(abs([MinPixel MaxPixel]));
                data = cat(3,data,data,data);
             
            % This option is for red/blue on a black background.  It has
            % been removed because it was buggy.
            case handles.RBKButton
                                % Normalize to largest pixel value (positive or negative)
                data = data/max(abs([MinPixel MaxPixel]));
                posmask = data>0;
                rMap = zeros(size(data));
                gMap = zeros(size(data));
                bMap = zeros(size(data));
                rMap(posmask) = data(posmask);
                bMap(~posmask) = -data(~posmask);
                data = cat(3,rMap,gMap,bMap);
                
        end
        
    end

    function c = elementdivide(a,b)
        c = a./b;
    end

    function c = elementmultiply(a,b)
        c = a.*b;
    end

    function RefreshGlobalImageParameters
        %         This function needs to be able to be called from within any image
        %         manipulation function (ACCouple, IsolateChannel, etc.) and change
        %         the MinPixel, MaxPixel, MinImage, MaxImage and anything else
        %         necessary around so that nothing breaks.  It should also de-check
        %         any checkboxes necessary for processes that should only go in one
        %         order, and in general be a handy way to simplify the other image
        %         functions.
        MeanImage = mean(ImageMatrix,3);
        MaxImage = max(ImageMatrix,[],3);
        MinImage = min(ImageMatrix,[],3);
        PixelRange = MaxImage - MinImage;

    end

    function ResetMinMaxPixels
        % For de-AC Coupling when channel is isolated the Min/Max pixels must be recalculated for color scaling
                MeanImage = mean(ImageMatrix,3);
        MaxImage = max(ImageMatrix,[],3);
        MinImage = min(ImageMatrix,[],3);
        
        if IsolateChannel
            if ACCouple
                NewMinPixel = min(MinImage(RFilter)-MeanImage(RFilter))
                NewMaxPixel = max(MaxImage(RFilter)-MeanImage(RFilter))
            else
                NewMinPixel = min(MinImage(RFilter))
                NewMaxPixel = max(MaxImage(RFilter))
            end
        else
            if ACCouple
                NewMaxPixel = max(MaxImage(:)-MeanImage(:))
                NewMinPixel = min(MinImage(:)-MeanImage(:))
            else
                NewMinPixel = min(MinImage(:))
                NewMaxPixel = max(MaxImage(:))
            end
        end
        
    end

end