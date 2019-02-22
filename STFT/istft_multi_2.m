function x=istft_multi_2(X,nsampl)

% ISTFT_MULTI Multichannel inverse short-time Fourier transform (ISTFT)
% using half-overlapping sine windows.
%
% x=istft_multi(X,nsampl)
%
% Inputs:
% X: nfram x nbin x nchan matrix containing STFT coefficients for nchan
% sources with nbin frequency bins and nfram time frames 
% nsampl: number of samples to which the corresponding time-domain signals
% are to be truncated
%
% Output:
% x:  nsampl x nchan matrix containing the corresponding time-domain signals
% If x is a set of signals of length nsampl and X=stft_multi(x), then
% x=istft_multi(X,nsampl).

%%% Errors and warnings %%%
if nargin<2, error('Not enough input arguments.'); end
[nfram,nbin,nchan]=size(X);
if nbin==2*floor(nbin/2), error('The number of frequency bins must be odd.'); end
wlen=2*(nbin-1);

%%% Computing inverse STFT signal %%%
% Defining sine window
win=sin((.5:wlen-.5)/wlen*pi).';

%%% windowing method 1:
swin=ones((nfram+1)*wlen/2,1);
% for t=0:nfram-1,
%     swin(t*wlen/2+1:t*wlen/2+wlen)=swin(t*wlen/2+1:t*wlen/2+wlen)+win.^2;
% end
% swin=sqrt(swin);
%%%
swin(1:wlen/2,1)=win(1:wlen/2);
swin(nfram*wlen/2+1:end,1)=win(wlen/2+1:wlen);

x=zeros((nfram+1)*wlen/2,nchan);
for i=1:nchan,
    for t=0:nfram-1,
        % IFFT
        fframe=[X(t+1,:,i),conj(X(t+1,wlen/2:-1:2,i))];
        frame=real(ifft(fframe)).';
        %%% Overlap-add method 1
        x(t*wlen/2+1:t*wlen/2+wlen,i)=x(t*wlen/2+1:t*wlen/2+wlen,i)+frame.*win./swin(t*wlen/2+1:t*wlen/2+wlen);
%         %%% Overlap-add method 2
%         x(t*wlen/2+1:t*wlen/2+wlen,i)=x(t*wlen/2+1:t*wlen/2+wlen,i)+frame.*win;
    end
end

%%% just for method 2: to keep the energy before and after windowing the same.
% swin=zeros((nfram+1)*wlen/2,1);
% for t=0:nfram-1,
%     swin(t*wlen/2+1:t*wlen/2+wlen)=swin(t*wlen/2+1:t*wlen/2+wlen)+win.^2;
% end
% x=x./swin;

% % Truncation
x=x(1:nsampl,:);
return;