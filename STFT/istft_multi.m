function x=istft_multi(X,nsampl)

% ISTFT_MULTI Multichannel inverse short-time Fourier transform (ISTFT)
% using half-overlapping sine windows.
%
% x=istft_multi(X,nsampl)
%
% Inputs:
% X: nbin x nfram x nsrc matrix containing STFT coefficients for nsrc
% sources with nbin frequency bins and nfram time frames or nbin x nfram x
% nsrc x nchan matrix containing the STFT coefficients of nsrc spatial
% source images over nchan channels
% nsampl: number of samples to which the corresponding time-domain signals
% are to be truncated
%
% Output:
% x: nsrc x nsampl matrix or nsrc x nsampl x nchan matrix containing the
% corresponding time-domain signals
% If x is a set of signals of length nsampl and X=stft_multi(x), then
% x=istft_multi(X,nsampl).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2008 Emmanuel Vincent
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Errors and warnings %%%
if nargin<2, error('Not enough input arguments.'); end
[nbin,nfram,nsrc,nchan]=size(X);
if nbin==2*floor(nbin/2), error('The number of frequency bins must be odd.'); end
wlen=2*(nbin-1);

%%% Computing inverse STFT signal %%%
% Defining sine window
win=sin((.5:wlen-.5)/wlen*pi);
% Pre-processing for edges
swin=zeros(1,(nfram+1)*wlen/2);
for t=0:nfram-1,
    swin(t*wlen/2+1:t*wlen/2+wlen)=swin(t*wlen/2+1:t*wlen/2+wlen)+win.^2;
end
swin=sqrt(swin/wlen);
x=zeros(nsrc,(nfram+1)*wlen/2,nchan);
for i=1:nchan,
    for j=1:nsrc,
        for t=0:nfram-1,
            % IFFT
            fframe=[X(:,t+1,j,i);conj(X(wlen/2:-1:2,t+1,j,i))];
            frame=real(ifft(fframe));
            % Overlap-add
            x(j,t*wlen/2+1:t*wlen/2+wlen,i)=x(j,t*wlen/2+1:t*wlen/2+wlen,i)+frame.'.*win./swin(t*wlen/2+1:t*wlen/2+wlen);
        end
    end
end
% Truncation
x=x(:,wlen/4+1:wlen/4+nsampl,:);

return;