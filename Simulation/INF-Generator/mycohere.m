function [Cxy, f] = mycohere(varargin)
%COHERE Coherence function estimate.
%   Cxy = COHERE(X,Y,NFFT,Fs,WINDOW) estimates the coherence of X and Y
%   using Welch's averaged periodogram method.  Coherence is a function
%   of frequency with values between 0 and 1 that indicate how well the
%   input X corresponds to the output Y at each frequency.  X and Y are 
%   divided into overlapping sections, each of which is detrended, then 
%   windowed by the WINDOW parameter, then zero-padded to length NFFT.  
%   The magnitude squared of the length NFFT DFTs of the sections of X and 
%   the sections of Y are averaged to form Pxx and Pyy, the Power Spectral
%   Densities of X and Y respectively. The products of the length NFFT DFTs
%   of the sections of X and Y are averaged to form Pxy, the Cross Spectral
%   Density of X and Y. The coherence Cxy is given by
%       Cxy = (abs(Pxy).^2)./(Pxx.*Pyy)
%   Cxy has length NFFT/2+1 for NFFT even, (NFFT+1)/2 for NFFT odd, or NFFT
%   if X or Y is complex. If you specify a scalar for WINDOW, a Hanning 
%   window of that length is used.  Fs is the sampling frequency which does
%   not effect the cross spectrum estimate but is used for scaling of plots.
%
%   [Cxy,F] = COHERE(X,Y,NFFT,Fs,WINDOW,NOVERLAP) returns a vector of freq-
%   uencies the same size as Cxy at which the coherence is computed, and 
%   overlaps the sections of X and Y by NOVERLAP samples.
%
%   COHERE(X,Y,...,DFLAG), where DFLAG can be 'linear', 'mean' or 'none', 
%   specifies a detrending mode for the prewindowed sections of X and Y.
%   DFLAG can take the place of any parameter in the parameter list
%   (besides X and Y) as long as it is last, e.g. COHERE(X,Y,'mean');
%   
%   COHERE with no output arguments plots the coherence in the current 
%   figure window.
%
%   The default values for the parameters are NFFT = 256 (or LENGTH(X),
%   whichever is smaller), NOVERLAP = 0, WINDOW = HANNING(NFFT), Fs = 2, 
%   P = .95, and DFLAG = 'none'.  You can obtain a default parameter by 
%   leaving it off or inserting an empty matrix [], e.g. 
%   COHERE(X,Y,[],10000).
%
%   See also PSD, CSD, TFE.
%   ETFE, SPA, and ARX in the Identification Toolbox.

%   Author(s): T. Krauss, 3-31-93
%   Copyright (c) 1988-98 by The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 1998/06/03 14:42:19 $

narginchk(2,7);
x = varargin{1};
y = varargin{2};
%[msg,nfft,Fs,window,noverlap,p,dflag]=psdchk(varargin(3:end),x,y);
%error(msg)
nfft = 256;
Fs =8000;
window=hanning(nfft);
noverlap = 0.75*nfft;
p = .95;
dflag = 'none';

% compute PSD and CSD
window = window(:);
n = length(x);		% Number of data points
nwind = length(window); % length of window
if n < nwind    % zero-pad x , y if length is less than the window length
    x(nwind)=0;
    y(nwind)=0;  
    n=nwind;
end
x = x(:);		% Make sure x is a column vector
y = y(:);		% Make sure y is a column vector
k = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
					% (k = fix(n/nwind) for noverlap=0)
index = 1:nwind;

Pxx = zeros(nfft,1); Pxx2 = zeros(nfft,1);
Pyy = zeros(nfft,1); Pyy2 = zeros(nfft,1);
Pxy = zeros(nfft,1); Pxy2 = zeros(nfft,1);
for i=1:k
    if strcmp(dflag,'none')
        xw = window.*x(index);
        yw = window.*y(index);
    elseif strcmp(dflag,'linear')
        xw = window.*detrend(x(index));
        yw = window.*detrend(y(index));
    else
        xw = window.*detrend(x(index),0);
        yw = window.*detrend(y(index),0);
    end
    index = index + (nwind - noverlap);
    Xx = fft(xw,nfft);
    Yy = fft(yw,nfft);
    Xx2 = abs(Xx).^2;
    Yy2 = abs(Yy).^2;
    Xy2 = Yy.*conj(Xx);
    Pxx = Pxx + Xx2;
    Pxx2 = Pxx2 + abs(Xx2).^2;
    Pyy = Pyy + Yy2;
    Pyy2 = Pyy2 + abs(Yy2).^2;
    Pxy = Pxy + Xy2;
    Pxy2 = Pxy2 + Xy2.*conj(Xy2);
end

% Select first half
if ~any(any(imag([x y])~=0)),   % if x and y are not complex
    if rem(nfft,2),    % nfft odd
        select = [1:(nfft+1)/2];
    else
        select = [1:nfft/2+1];   % include DC AND Nyquist
    end
    Pxx = Pxx(select);
    Pxx2 = Pxx2(select);
    Pyy = Pyy(select);
    Pyy2 = Pyy2(select);
    Pxy = Pxy(select);
    Pxy2 = Pxy2(select);
else
    select = 1:nfft;
end
%Coh = (abs(Pxy).^2)./(Pxx.*Pyy);             % coherence function estimate 
Coh = Pxy./sqrt(Pxx.*Pyy);
freq_vector = (select - 1)'*Fs/nfft;

% set up output parameters
if (nargout == 2),
   Cxy = Coh;
   f = freq_vector;
elseif (nargout == 1),
   Cxy = Coh;
elseif (nargout == 0),   % do a plot
   newplot;
   plot(freq_vector,Coh), grid on
   xlabel('Frequency'), ylabel('Coherence Function Estimate');
end