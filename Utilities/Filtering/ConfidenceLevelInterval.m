function [Fconf, Pconf] = ConfidenceLevelInterval(signal,freq,freqbins,hilo,RGBColor,lnstyle,sclfactor)
%% Comments
%
% hmid = conflevelint(signal,freq,freqbins,filtwin,hilo,RGBcolor,lnstyle,sclfactor)
%
% This function generates a confidence interval band on an FFT by assuming
% that within each frequency band, of width given by the number of elements
% in signal divided by the number of frequency bins given by BINS, the points
% in the FFT are actually a normal distribution about some true mean value.
% The upper and lower bounds of the 95% confidence interval are then
% returned in LOWLIMIT and HILIMIT.

%% Initialization
showplots = 0;

if nargin<4;     hilo = 0; end
if nargin<5;     RGBColor = [1 0 0]; end
if nargin<6;     lnstyle = '-'; end
if nargin<7;     sclfactor = 1; end



lo95 = zeros(1,length(freqbins)-1);
hi95 = lo95;
mid = lo95;
fout = lo95;


%% Confidence Level Calculation
for i=1: (numel(freqbins)-1)
    % Choose range of frequencies this bin will cover
    lofreq = freqbins(i);
    hifreq = freqbins(i+1);
    inrange = freq > lofreq & freq < hifreq;
    % Return center frequency in this range
        fout(i) = (lofreq + hifreq)/2;

    
    % Idenfify FFT values in this range of frequencies
    s = signal(inrange);
%     numel(s)
    
    % Take statistics in this range (mean and upper/lower 95% confidence
    % intervals based on standard distribution assumption)
    sigma = std(s); % note this goes to zero if only one element in 's'
    mu = mean(s);
    stderror = sigma/sqrt(numel(s));
    factor95 = 1.96;    % This factor is from the standard error math, I think it's an output of the erf() function
    lo95(i) = mu - factor95*stderror;
    hi95(i) = mu + factor95*stderror;
    mid(i) = mu;
    
end

% Sometimes confidence interval goes into negatives; this is bad on a
% log-log plot.  Remove these elements to enable good plotting
foutlo = fout(lo95>0);        lo95 = lo95(lo95>0);
fouthi = fout(hi95>0);        hi95 = hi95(hi95>0);
foutmid = fout(mid>0);      mid = mid(mid>0);

lo95 = lo95*sclfactor;
mid = mid*sclfactor;
hi95 = hi95*sclfactor;

% Store frequencies and power densities in a structure so they can be
% easily passed back out to calling routines
Fconf.lo = foutlo;
Fconf.mean = foutmid;
Fconf.hi = fouthi;

Pconf.lo = lo95;
Pconf.mean = mid;
Pconf.hi = hi95;

if showplots
hmid = line(foutmid,mid);       set(hmid,'LineWidth',2,'Color',RGBColor,'LineStyle',lnstyle)
if hilo
    hlo  = line(foutlo, lo95 );       set(hlo,'LineWidth',1,'Color',.75*RGBColor,'LineStyle',lnstyle)
    hhi  = line(fouthi, hi95 );       set(hhi,'LineWidth',1,'Color',.75*RGBColor,'LineStyle',lnstyle)
end
end
