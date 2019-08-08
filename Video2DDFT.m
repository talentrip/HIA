function OutputMatrix = Video2DDFT(fn,reverse,ModeSelect,FreqSelect,RangeSelect,Removal,MaxFitWidth)

if nargin==0
    [filename, pathname] = uigetfile( ...
        {'*.mat','MATLAB MAT-files (*.mat)'}, ...
        'Pick a file');
    fn = [pathname filename]
    
    reverse = 0;
    close all; clc
end

if reverse
    ReverseString = ' - Reverse';
else
    ReverseString = '';
end

%% Load in data

[filepath,filename,fileext] = fileparts(fn);
load(fn);
fps = ImageInfo.Fps;

% The reverse command is intended to look for spokes propagating in the
% (clockwise? counter-clockwise?  Opposite from normal, anyway) direction.
% It is a little easier to just flip the spokeSurface matrix and analyze as normal
% than to try to pull out the relevant pieces from the 2D DFT and code all
% the confidence interval stuff for another case.
if reverse
    spokeSurface = fliplr(spokeSurface);
end

%% Initialize constants

% Set x-axis limits for pwoer spectrum plot
xlimskHz = [0 fps/2e3];

% Define an output matrix for FWHMs, Q-values, etc.  Will not be filled
% unless lorentzian fits are performed.
OutputMatrix = [];

% Switches
showplots = 1;
showCircleFit = 0;
verbose = 1;

% These empty variables are declared since this function has a nested
% helper function MakeConfidPlot.  These particular variables only appear
% once they are loaded from an external file, but if they aren't declared
% here the memory sharing can get tricky.
fullimgfn = '';
cmavg = [];
AzBinPixels = [];
IsolatedImage = [];
MeanImage = [];

% RGB Color Definitions
white   = [1 1 1];          grey   = 0.75*white;    red     = [1 0 0];
orange  = [1 165/255 0];    yellow = [1 1 0];       yellowgreen = [.5 1 0];
green   = [0 1 0];          cyan   = [0 1 1];       blue    = [0 0 1];
blueviolet = [.5 0 1];      violet = [1 0 1];       black = [0 0 0];
% Choose colors to label each mode, 0-8
ModeColors = {black,red,orange,yellow,green,cyan,blue,violet,0.5*violet};


%% Do FFTs

% 2D FFT raw (unchopped)
[mrange,freqsraw,Zraw] = do_FFTs(spokeSurface,1,fps,0);
freqsraw = freqsraw/1e3;    % convert to kHz

% Choose centers for frequency bins when smoothing spectrum. We desire a
% good resolution, but with a minimum number of points in each bin.  This
% is a minimum resolution of the confidence interval plot (too small
% and no FFT spectrum datapoints will fall in the frequency bin) that
% depends on the video length.

% Choose a desired frequency resolution; 50 Hz would be nice
DesiredDf = 5e1;
% Compute the min. frequency resolution of full video (1 / video duration)
% and double it.
MinDf = 2/ImageInfo.Duration;
% Now set the frequency resolution as the larger of these two numbers
df = max(MinDf,DesiredDf);
% And compute the boundaries between each frequency bin
freqbins = [0:df:xlimskHz(2)*1e3]/1e3;

if showplots
    %% Video DFT for Breathing and Spoke Modes
    
    % Show fit of discharge channel used to generate spoke surface
    if showCircleFit
        VisualizeCircleFit(fn);
    end
    
    Nmodes = length(ModeColors); % there are 9 colors above for modes m=0-8
    
    % Crop the DFT to the range of frequencies given by xlimskHz above
    [dum1,dum2,FrequencyRange] = ConfidenceLevelCrop(Zraw(:,1),freqsraw,freqbins,xlimskHz);
    
    % Eliminate frequencies > Nyquist frequency of the video framerate
    FrequencyRange(FrequencyRange>(fps/(2*1e3))) = [];
    
    % Plot modes
    if nargin<3
        % Create handles for each mode plot
        plotHandles = zeros(1,Nmodes);
        
        % Label figure for all spoke modes and whether it is computed for
        % forward or reverse (cw or ccw) propagation
        figure('Name',['All Spoke Modes' ReverseString],'NumberTitle','off')
        
        for i=1:Nmodes
            % Make plot
            [h{i} Fconf Pconf] = MakeConfidPlot(Zraw,i,freqsraw,freqbins,xlimskHz,ModeColors{i});
            
            % acquire the handle for the midline plot
            plotHandles(i) = h{i}.mid;
            
            % construct the legend label, and note it as a "breathing" or
            % "spoke" mode
            plotLabels{i} = ['m=' num2str(i-1)];
            if i==1; plotLabels{i} = [plotLabels{i} ' Breathing Mode']; end
            if i>=2; plotLabels{i} = [plotLabels{i} ' Spoke Mode']; end
        end
        
        
    else
        % In the case with >3 arguments, the function is being called to
        % fit Lorentzian (damped driven oscillator) functions to the power
        % spectra.  This requires input guesses for a nonlinear solver.
        plotHandles = zeros(1,numel(ModeSelect));
        OutputMatrix = zeros(5,numel(ModeSelect));
        A_scaled = zeros(1,numel(ModeSelect));
        
        % Set figure window label
        figure('Name',['Chosen Spoke Modes' ReverseString],'NumberTitle','off')
        
        for i=1:numel(ModeSelect)
            % Create DFT plot with confidence intervals
            [h{i} Fconf Pconf] = MakeConfidPlot(Zraw,ModeSelect(i)+1,freqsraw,freqbins,xlimskHz,ModeColors{ModeSelect(i)+1});
            
            % acquire the handle for the midline plot
            plotHandles(i) = h{i}.mid;
            
            % construct the legend label, and note it as a "breathing" or
            % "spoke" mode
            plotLabels{i} = ['m=' num2str(ModeSelect(i))];
            if ModeSelect(i)==0;
                plotLabels{i} = [plotLabels{i} ' Breathing Mode'];
            else
                plotLabels{i} = [plotLabels{i} ' Spoke Mode'];
            end
            
            % Remove problematic sections of the power spectrum for the
            % given mode.  This is most often a problem when a large
            % breathing mode bleeds into the spoke mode spectrum and must
            % be excised to fit a Lorentzian to the spoke response.
            Pconf.mean(Removal{i}) = [];
            Fconf.mean(Removal{i}) = [];
            
            % Create copies of the mean (i.e., the smoothed) power spectrum
            % vectors and crop them to the region of interest where the
            % peak has been manually identified by FreqSelect and
            % RangeSelect
            f = Fconf.mean;
            p = Pconf.mean;
            inrange = abs(f-FreqSelect(i))<RangeSelect(i);
            f = f(inrange);
            p = p(inrange);
            
            % Display the current mode to the standard out on the command
            % window
            if verbose
                disp(plotLabels{i})
            end
            
            % And finally fit the Lorentzian (lorentzianfit will construct
            % a text output to the command window)
            
            [F0,FWHM0,A0,Q0] = lorentzianfit(f,p,0,0);
            
            % Now re-compute the Lorentzian using these initial values
            % to refine the fit:
            % First, re-copy the power spectrum
            f = Fconf.mean;
            p = Pconf.mean;
            
            % Choose a fitting width for the Lorentzian
            % function of at least 3 kHz and at most the
            % smaller of MaxFitWidth (a manual input) or 4
            % times the FWHM determined above.
            FitWidth = max(min(MaxFitWidth,4*FWHM0),3);
            
            % Crop the power spectrum to the computed frequency range
            inrange = abs(f-F0)<FitWidth/2;
            f = f(inrange);
            p = p(inrange);
            
            % Compute the new Lorentzian
            [F,FWHM,A,Q] = lorentzianfit(f,p,0,0);
            
            % Check frequency F against initial guess; if off, repeat with
            % half the FitWidth
            FreqErrorThreshold = 1; %kHz
            FreqGuessError = abs(F-FreqSelect(i))
            while FreqGuessError>FreqErrorThreshold
                f = Fconf.mean;
                p = Pconf.mean;
                
                FitWidth = FitWidth/2;
                inrange = abs(f-F0)<FitWidth/2;
                f = f(inrange);
                p = p(inrange);
                
                [F,FWHM,A,Q] = lorentzianfit(f,p,0,0);
                
                % Check that the new frequency peak identification error is
                % less than the old -- otherwise, this loop is prone toward
                % runaway failure.
                if abs(F-FreqSelect(i))<FreqGuessError
                    FreqGuessError = abs(F-FreqSelect(i))                    
                else
                    disp(['m=' num2str(ModeSelect(i))' ' frequency guess unable to converge near guess.  Refine guess.'])
                    % now undo the damage done in this loop and restore the
                    % initial fit results from before the loop started
                    FreqGuessError = 0;
                    FitWidth = FitWidth*2;
                    
                    f = Fconf.mean;
                    p = Pconf.mean;
                    
                    inrange = abs(f-F0)<FitWidth/2;
                    f = f(inrange);
                    p = p(inrange);
                end                                
            end
            
            % This final duplication is with the 'showplots' switch
            % enabled, to illustrate the fit to the user
            [F,FWHM,A,Q] = lorentzianfit(f,p,0,1);
            
            
            A_scaled(i) = sqrt(A)/180;
            
            % And compute the error between the initial and refined
            % estimates
            Error = ([F,FWHM,A,Q]-[F0,FWHM0,A0,Q0])./[F,FWHM,A,Q];
            
            OutputMatrix(:,i) = [ModeSelect(i); F; FWHM; A; Q];
            
        end
        A_scaled;
    end
    
    
    ylabel({'Power Spectral Density';'Pixel Amplitude^2 / Hz'})
    xlabel('Frequency, kHz')
    % assign the iteratively constructed legend labels to the various
    % midline plots
    legend(plotHandles,plotLabels,'Location','NorthEast')
    grid on
    set(gca,'YScale','log','XScale','linear','YMinorGrid','off')
    title(filename)
    
    
    OutputMatrix = OutputMatrix';
end

%% Helper Functions
function [h F P] = MakeConfidPlot(Z,Zcol,f,fbins,xlimskHz,rgbcolor)

[f,Z,fbins] = ConfidenceLevelCrop(Z(:,Zcol),f,fbins,xlimskHz);
[F,P] = ConfidenceLevelInterval(Z,f,fbins);
h = ConfidenceLevelPlot(F,P,rgbcolor);

