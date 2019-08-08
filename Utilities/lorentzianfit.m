function [F,FWHM,A,Q] = lorentzianfit(x,y,DHO,showplots)
%% Comments
%  This function is a quick way to quantify spoke oscillation frequency,
%  amplitude, and full width at half maximum (FWHM).  It fits a Lorentzian
%  function (a forced, damped oscillator model) to a spoke's peak in the
%  FFT, and pulls out the peak frequency F, full-width at half-maximum
%  FWHM, amplitude A, and quality factor Q.
%
%  The code includes an option to fit the full damped harmonic oscillator
%  function to the curve if desired, first using the results from a
%  Lorentzian fit to seed the solver for the DHO model.  To see this fit,
%  set DHO = 1.
%
%  Note: The amplitude of the peak is a better measure of spoke strength
%  than the area.  A large amplitude and small FWHM sometimes give the same
%  area as a small amplitude and large FWHM, but the large amplitude case
%  is the one you care about.  Alternately, 1/FWHM can also be used as a
%  measure of spoke strength (narrow peaks are usually strong).
%  Mathematically, 1/FWHM is related to the Q (quality factor) of the
%  oscillator in the forced oscillator model.
%
%  However, this function should mainly be used to quantify the frequency
%  peaks in the DFT, rather than spoke strength.  The ratio of RMS values
%  in the frequency space discussed in the RSI paper is a more robust
%  method for evaluating spoke strength (it doesn't require initial
%  guesses, for instance).

%% Initialization

% Define the functional form of the resonance, using a 3-element vector 'a'
% and a variable 'x', where 'x' is the independent variable for frequency
% a = [Amplitude, Freq of peak,FWHM]
% b = [a(1)*(2*a(2)*a(3))^2, a(2), 2*a(3)];
Lorentzian = @(a, x) a(1)./(1 + (a(2) - x).^2/a(3)^2);
DampHarmOsc = @(b, x) b(1)./((b(2)^2 - x.^2).^2 + b(3)^2*x.^2);

if nargin<4
    showplots = 0;
end
if nargin<3
    DHO = 0;
end

% Illustrate the different fit types if nargin==0
if nargin==0
    close all; clc
    figure()
    showplots = 1;
    DHO = 1;
    
    % Choose parameters for a sample damped harmonic oscillator (DHO)
    % signal
    omega_0 = 5;    % undamped oscillation frequency
    gam = 2;      % full-width at half-maximum (FWHM) of the resonance
    % note that gamma = omega_0 / Q, where Q = the quality factor of the
    % resonance
    FoverMsqrd = 100;   % (Fo/m)^2
    
    % Create sample DHO signal
    x = linspace(0, 2*omega_0, 201);
    SimpleLorentzian = (FoverMsqrd*1/omega_0^2*1/gam^2) ./ (1+(omega_0-x).^2/(gam/2)^2);
    SimpleDHO = FoverMsqrd ./ ((omega_0^2-x.^2).^2 + gam^2*x.^2);
    
    % Add noise proportional to local value of signal to corrupt it but
    % leave it similar to power spectra obtained in thruster discharge
    % current traces and high speed videos
    %
    % Note: the default random number seed is called so the noise is the
    % same on every run. This is useful to check that the least squares
    % fit is properly converging on the same value each time.
    randSeed = RandStream.create('mrg32k3a');
    
    noise = .25*randn(randSeed,size(SimpleDHO)).*SimpleDHO;
    
    % use abs. value so log plots don't get angry about plotting tiny
    % negative values
    y = abs(SimpleDHO + noise);
end

% Force input vectors to columns
x = x(:);
y = y(:);

%% Make Initial Guess for Nonlinear Fit

% Make an initial guess for the Lorentzian parameters to seed the solving
% routine.  Use the location of the peak in the spectrum to gauge amplitude
% and frequency, and guess a FWHM of unity (in kHz this turns out to be a
% good guess for spokes and breathing modes)
[m, ii] = max(y);
gamma_0 = 1;
a0 = [m, x(ii), ( gamma_0 /2)];

%% Perform nonlinear fit for Lorentzian
% do a non-linear least-squares fit of the Lorentzian to the measured data
opts = optimset('lsqcurvefit');
opts.Display = 'off';
a = lsqcurvefit(Lorentzian, a0, x, y,[],[],opts);

% Compute oscillator parameters from fit parameters
a=abs(a);
F = a(2);
FWHM = 2*a(3);
A = a(1);
Q = F/FWHM;

%% Use Lorentzian fit to perform DHO fit, if desired
if DHO
    % Now do an initial guess for the full damped harmonic oscillator model
    % fit, using the Lorentzian fit to compute an initial guess
    b0 = [a(1)*(2*a(2)*a(3))^2, a(2), 2*a(3)];
    opts = optimset('lsqcurvefit');
    opts.Display = 'off';
    b = lsqcurvefit(DampHarmOsc, b0, x, y,[],[],opts);
    
    b=abs(b);
    
    % Compute the DHO parameters from the fit
    Q = b(2)/b(3);
    A0_sqrd = b(1)/b(2)^4;
    Qfactor_sqrd = Q^2/(1-1/(4*Q^2));
    Am_sqrd = A0_sqrd*Qfactor_sqrd;
    omega_m = b(2)*sqrt(1-1/(2*Q^2));
    
    % Assign to convenient variables to return to end user
    F = omega_m;  % F for frequency of the peak
    FWHM = b(3);   % the Full-Width at Half-Max
    A = Am_sqrd;    % A for amplitude (really amplitude squared since it's the power spectral density
    Q = Q;          % redundant but just for completeness
end

if nargin==0
    % plot the results for the demo case
    subplot(1,2,1)
    plot(x,SimpleLorentzian,'ro-')
    hold on
    plot(x,SimpleDHO,'bo-')
    hold off
    set(gca,'YScale','log')
    legend('pure Lorentzian', 'pure DHO', 'Location', 'South');
    title({'Pure signals for Lorentzian and DHO';'computed using input frequency,';'FWHM and amplitude'})
    
    subplot(1,2,2)
    plot(x, y, 'bo-');
    hold all
    % plot(x, signal, '--');
    plot(x, Lorentzian(a, x), 'ro-');
    plot(x,DampHarmOsc(b,x), 'go-');
    plot(x,SimpleDHO,'k-', 'LineWidth',2)
    hold off
    set(gca,'YScale','log')
    legend('noisy signal', 'Lorentzian fit', 'DHO fit', 'pure signal', ...
        'Location', 'South');
    title({'DHO and Lorentzian fits to noisy DHO signal,';'with pure DHO signal also shown';'for comparison'})
    
    % Compute what the values *should* have been for the pure DHO signal
    % before it was corrupted by noise
    
    Q0 = omega_0/gam;
    A00_sqrd = b0(1)/omega_0^4;
    Q0factor_sqrd = Q0^2/(1-1/(4*Q0^2));
    Am0_sqrd = A00_sqrd*Q0factor_sqrd;
    omega_m0 = omega_0*sqrt(1-1/(2*Q0^2));
    fprintf('Pure values: F = %3.4G,    Gamma = %3.3G,    Amp = %3.3G,    Q = %3.3G \r',omega_m0,gam,Am0_sqrd,Q0)
else
    if showplots
        if DHO            
            h = line(x,DampHarmOsc(b,x));
        else
            h = line(x,Lorentzian(a,x));            
        end
        set(h,'Marker','o','Color',[1 0 0],'LineStyle','-')
    end
end


%% Print results to user
if showplots
fprintf('Fit values: F = %3.4G,    Gamma (FWHM) = %3.3G,    Amp = %3.3G,    Q = %3.3G \r',F,FWHM,A,Q)
end