function [F,FWHM,A,Area] = lorentzianfit(x,y)
%% Comments
%  This function is a quick-and-dirty way to quantify spoke oscillation
%  strength. It fits a Lorentzian function (a forced, damped oscillator
%  model) to a spoke's peak in the FFT, and pulls out the peak frequency F,
%  full-width at half-maximum FWHM, amplitude A, and area under the curve
%  AREA.
%
%  Note: The amplitude of the peak is a better measure of spoke strength
%  than the area.  A large amplitude and small FWHM sometimes give the same
%  area as a small amplitude and large FWHM, but the large amplitude case
%  is the one you care about.  Alternately, 1/FWHM can also be used as a
%  measure of spoke strength (narrow peaks are usually strong).
%  Mathematically, 1/FWHM is related to the Q (quality factor) of the
%  oscillator in the forced oscillator model.

%% Code

% define the functional form of the resonance, using a 3-element vector a
% and a variable x.  
% a = 
% b = 
% and x is the independent variable for frequency

DampHarmOsc = @(b, x) b(1)./((b(2)^2 - x.^2).^2 + b(3)^2*x.^2);

if nargin==0
    close all
    figure()
    % make a signal and noise:
    x = linspace(60, 80, 101);
    omega_0 = 70;
    gamma = 3;
    MaxAmp = 20;
    b0 = [MaxAmp, omega_0, gamma]
    

    signal2 = DampHarmOsc(b0,x);
    

    max(signal2(:))

    
    
    noise = .1*max(signal2(:))*randn(size(signal2));
    y = signal2 + noise;
end


% plot(x,DampHarmOsc([70 3 20],x))
% return
x = x(:);
y = y(:);
% Make an initial guess for the Lorentzian parameters to seed the solving
% routine.  Use the location of the peak in the spectrum to gauge amplitude
% and frequency, and guess a FWHM of unity (in kHz this turns out to be a
% good guess for spokes and breathing modes)
[m, ii] = max(y);

b0 = [m, x(ii), 1]

% do a least-squares fitting of the Lorentzian to the measured data
opts = optimset('lsqcurvefit');
opts.Display = 'off';
b = lsqcurvefit(DampHarmOsc, b0, x, y,[],[],opts)

if nargin==0
        % plot the results
    subplot(2,1,1)
    % plot(x, signal, '--');
    plot(x, signal2, 'ro-');   
    legend('pure DHO','Location', 'Best');

    % plot the results
    subplot(2,1,2)
    plot(x, y, 'bo-');
    hold all
    % plot(x, signal, '--');
    plot(x, DampHarmOsc(b, x), 'ro-');
    plot(x,signal2,'g+')
    hold off
    
    legend('noisy signal', 'fitted curve', 'pure signal', ...
        'Location', 'Best');
    
else
    
    h = line(x,Lorentzian(a,x));
    set(h,'Marker','o','Color',[1 0 0],'LineStyle','-')
    a=abs(a);
    F = a(1);
    FWHM = a(2);
    A = a(3);
    Area = A*pi*FWHM;
    fprintf('F = %3.4G,    Gamma = %3.3G,    Amp = %3.3G, Area = %3.3G \r',F,FWHM,A,Area)
    %     disp('')
end