function z = sinf_3D(P,len,params)

% Generating sensor signals for a 1D sensor array in a spherically 
% isotropic noise field [1,2]
%
%    z = sinf_3D(P,len,params)
%
% Input parameters:
%    P            : sensor positions
%    len          : desired data length
%    params.fs    : sample frequency
%    params.c     : sound velocity in m/s
%    params.N_phi : number of cylindrical angles
%
% Output parameters:
%    z            : output signal
%
% References:
%    [1] E.A.P. Habets and S. Gannot, 'Generating sensor signals
%        in isotropic noise fields', The Journal of the Acoustical 
%        Society of America, Vol. 122, Issue 6, pp. 3464-3470, Dec. 2007.
%    [2] E.A.P. Habets and S. Gannot, 'Comments on Generating 
%        sensor signals in isotropic noise fields', Technical Report,
%        Imperial College London, Sept. 2010.
%
% Authors:  E.A.P. Habets and S. Gannot
%
% History:  2007-11-02 - Initial version
%           2010-09-16 - Near-uniformly distributed points over S^2
%           2017-20-06 - Use native waitbar
%
% Copyright (C) 2007-2017 E.A.P. Habets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = size(P,2);                        % Number of sensors
NFFT = 2^ceil(log2(len));             % Number of frequency bins
X = zeros(M,NFFT/2+1);

if ~isfield(params,'fs')
    fs = 8000;                        % Default
else
    fs = params.fs;
end
if ~isfield(params,'c')               
    c = 340;                          % Default
else
    c = params.c;
end
if ~isfield(params,'N')                 
    N = 256;                          % Default 
else
    N = params.N;
end

w = 2*pi*fs*(0:NFFT/2)/NFFT;

% Generate N points that are near-uniformly distributed over S^2
for k=1:N
    h = -1 + 2*(k-1)/(N-1);
    phi(k) = acos(h);
    if k==1 || k==N
        theta(k) = 0;
    else
        theta(k) = mod((theta(k-1) + 3.6/sqrt(N*(1-h^2))),2*pi);
    end
end

% Caculate relative positions
P_rel = zeros(3,M);
for m = 1:M
    P_rel(:,m) = P(:,m) - P(:,1);
end

% Initialize waitbar
% h = waitbar(0,'Generating sensor signals...');

% Calculate sensor signals in the frequency domain
for idx = 1:N
%     waitbar(idx/N);
    
    X_prime = randn(1,NFFT/2+1) + 1i*randn(1,NFFT/2+1);
    X(1,:) = X(1,:) + X_prime;
    for m = 2:M
        v = [cos(theta(idx))*sin(phi(idx)) ; sin(theta(idx))*sin(phi(idx)) ; cos(phi(idx))];
        Delta = v.' * P_rel(:,m);
        X(m,:) = X(m,:) + X_prime.*exp(-1i*Delta*w/c);
    end
end
X = X/sqrt(N);

% Transform to time domain
X = [sqrt(NFFT)*real(X(:,1)), sqrt(NFFT/2)*X(:,2:NFFT/2),...
    sqrt(NFFT)*real(X(:,NFFT/2+1)), sqrt(NFFT/2)*conj(X(:,NFFT/2:-1:2))];
z = real(ifft(X,NFFT,2));

% Truncate output signals
z = z(:,1:len);

% Close waitbar
% close(h);