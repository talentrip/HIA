function [center R InnerRadius OuterRadius IsoMeanImage] = CircleFitGui(meanimg,GUIdir)

%% Comments
% clc
% debug = 0;
if nargin==0
disp('Error, no inputs received')
return
end
% return
%% Open GUI Figure and assign Callback Functions
% To reveal the tags of all the objects in the GUI figure window,
% un-suppress the output of the handles = guihandles(fig) command.
% Note: the '@' indicates a so-called "nested function" that is explicitly
% located later in this m-file (see "Callbacks" code section)
if isdeployed
    GUIFigFile = 'Circle_Translation_Tool.fig';
else
    homedir = userpath;
% homedir = homedir(1:end-1);
homedir = homedir(1:end);
GUIdir = fullfile(homedir, 'High Speed Image Analysis\GUI');
    GUIFigFile = fullfile(GUIdir, 'Circle_Translation_Tool.fig');
end

fig = openfig(GUIFigFile);
handles = guihandles(fig);

set(handles.Done, 'CallBack', @DoneCallback);
set(handles.Up, 'CallBack', @UpCallback);
set(handles.Down, 'CallBack', @DownCallback);
set(handles.Right, 'CallBack', @RightCallback);
set(handles.Left, 'CallBack', @LeftCallback);
set(handles.TextBox, 'CallBack', @TextBoxCallback);
set(handles.RadiusSlider, 'CallBack', @RadiusSliderCallback);

set(handles.Up, 'enable', 'on');
set(handles.Down, 'enable', 'on');
set(handles.Right, 'enable', 'on');
set(handles.Left, 'enable', 'on');
set(handles.RadiusSlider, 'enable', 'off');

%% Set Default Properties and Global Variables

set(handles.RadiusSlider,'Min',0);
set(handles.RadiusSlider,'Max',max(size(meanimg)));

FittingStep = 1;
InnerRadius = 0;
OuterRadius = 0;
PermitExit = 0;

%% Program Main

[a b R] = TaubinSVDThresholdFit(meanimg);
center = [a b];
% pause
% plot(handles.Axis1,1:5,1:5)
% pause
PlotCirc(meanimg,a,b,R,handles.Axis1,'k');
while PermitExit == 0
    pause(0.1)
end
delete(handles.figure1);
return;

%% Callback Functions

    function UpCallback(h, dummy)
        b = b + 1;
        PlotCirc(meanimg,a,b,R,handles.Axis1,'k');
    end

    function DownCallback(h, dummy)
        b = b - 1;
        PlotCirc(meanimg,a,b,R,handles.Axis1,'k');
    end

    function RightCallback(h,dummy)
        a = a + 1;
        PlotCirc(meanimg,a,b,R,handles.Axis1,'k');
    end

    function LeftCallback(h, dummy)
        a = a - 1;
        PlotCirc(meanimg,a,b,R,handles.Axis1,'k');
    end

    function RadiusSliderCallback(h, dummy)
        SliderValue = get(handles.RadiusSlider,'Value');
        if (FittingStep == 2)
            InnerRadius = SliderValue;
            PlotCirc(meanimg,a,b,InnerRadius,handles.Axis1,'k');
        elseif (FittingStep == 3)
            OuterRadius = SliderValue;
            PlotCirc(meanimg,a,b,OuterRadius,handles.Axis1,'k');
        end
    end

    function DoneCallback(h, dummy)
        if (FittingStep == 1)
            set(handles.RadiusSlider,'Value',0.75*R);
            PlotCirc(meanimg,a,b,0.75*R,handles.Axis1,'k');
            set(handles.TextBox,'String','2) Use slider to contract the inner radius bounds until satisfied, then press Enter');
            set(handles.RadiusSlider, 'enable', 'on');
            % Disable use of directional buttons after the position of the
            % center of the circle is chosen
            set(handles.Up, 'enable', 'off');
            set(handles.Down, 'enable', 'off');
            set(handles.Right, 'enable', 'off');
            set(handles.Left, 'enable', 'off');
            SliderValue = get(handles.RadiusSlider,'Value');
            InnerRadius = SliderValue;
            PlotCirc(meanimg,a,b,InnerRadius,handles.Axis1,'k');
            
        elseif (FittingStep == 2)
            set(handles.RadiusSlider,'Value',1.25*R);
            PlotCirc(meanimg,a,b,1.25*R,handles.Axis1,'k');
            set(handles.TextBox, 'String','3) Use slider to expand the outer radius bounds until satisfied, then press Enter');
            SliderValue = get(handles.RadiusSlider,'Value');
            OuterRadius = SliderValue;
            PlotCirc(meanimg,a,b,OuterRadius,handles.Axis1,'k');
            
            
        elseif (FittingStep == 3)
            set(handles.Done, 'String', 'Done');
            set(handles.TextBox, 'String','4) Thruster circle fit and cropping complete, press Done.');
            set(handles.RadiusSlider,'enable','off')
            r = PixelRadiiFromOrigin(meanimg,a,b);
            rfilter = r<InnerRadius | r>OuterRadius;
            meanimg(rfilter) = 0;
            [a,b,R] = TaubinSVDThresholdFit(meanimg);
            PlotCirc(meanimg,a,b,R,handles.Axis1,'k');
            center = [a b];
            IsoMeanImage = meanimg;
            
        end
        
        if (FittingStep == 4)
            PermitExit = 1;
        end
        
        FittingStep = FittingStep + 1;
    end

end


