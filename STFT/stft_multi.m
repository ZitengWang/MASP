function X=stft_multi(x,wlen)

% STFT_MULTI Multichannel short-time Fourier transform (STFT) using
% half-overlapping sine windows.
%
% X=stft_multi(x)
% X=stft_multi(x,wlen)
%
% Inputs:
% x: nchan x nsampl matrix containing nchan time-domain mixture signals
% with nsampl samples
% wlen: window length (default: 1024 samples or 64ms at 16 kHz, which is
% optimal for speech source separation via binary time-frequency masking)
%
% Output:
% X: nbin x nfram x nchan matrix containing the STFT coefficients with nbin
% frequency bins and nfram time frames
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2008 Emmanuel Vincent
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Errors and warnings %%%
if nargin<1, error('Not enough input arguments.'); end
if nargin<2, wlen=1024; end
[nchan,nsampl]=size(x);
if nchan>nsampl, error('The signals must be within rows.'); end
if wlen~=4*floor(wlen/4), error('The window length must be a multiple of 4.'); end

%%% Computing STFT coefficients %%%
% Defining sine window
win=sin((.5:wlen-.5)/wlen*pi).';
% Zero-padding
nfram=ceil(nsampl/wlen*2);
x=[x,zeros(nchan,nfram*wlen/2-nsampl)];
% Pre-processing for edges
x=[zeros(nchan,wlen/4),x,zeros(nchan,wlen/4)];
swin=zeros((nfram+1)*wlen/2,1);
for t=0:nfram-1,
    swin(t*wlen/2+1:t*wlen/2+wlen)=swin(t*wlen/2+1:t*wlen/2+wlen)+win.^2;
end
swin=sqrt(wlen*swin);
nbin=wlen/2+1;
X=zeros(nbin,nfram,nchan);
for i=1:nchan,
    for t=0:nfram-1,
        % Framing
        frame=x(i,t*wlen/2+1:t*wlen/2+wlen).'.*win./swin(t*wlen/2+1:t*wlen/2+wlen);
        % FFT
        fframe=fft(frame);
        X(:,t+1,i)=fframe(1:nbin);
    end
end

return;