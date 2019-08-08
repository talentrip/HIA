function h = ConfidenceLevelPlot(F,P,rgbcolor,hilo,lnstyle)

% F = frequency vector in struct format with fields mean, hi and lo
% P = power spectral density vector, same fields
% rgbcolor = 3-element row vector in [R G B] format with values on (0,1)
% lnstyle = LineStyle argument, i.e., '-', '.', '.-', etc. 
if nargin<5; lnstyle = '-'; end
if nargin<4; hilo = 1; end
if nargin<3; rgbcolor = [1 0 0]; end

% Split-core Hall probe discharge current
h.mid = line(F.mean,P.mean);       set(h.mid,'LineWidth',2,'Color',rgbcolor,'LineStyle',lnstyle)
if hilo
    hilolnstyle = '-';
    h.lo  = line(F.lo,P.lo);       set(h.lo,'LineWidth',1,'Color',.75*rgbcolor,'LineStyle',hilolnstyle)
    h.hi  = line(F.hi,P.hi);       set(h.hi,'LineWidth',1,'Color',.75*rgbcolor,'LineStyle',hilolnstyle)
end