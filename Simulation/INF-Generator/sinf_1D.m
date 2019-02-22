function z = sinf_1D(d,len,params)

% Generating sensor signals for a 1D sensor array in a spherically 
% isotropic noise field [1]
%
%    z = sinf_1D(d,len,params)
%
% Input parameters:
%    d            : sensor distances
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
%
% Authors:  E.A.P. Habets and S. Gannot
%
% History:  2007-11-02 - Initial version
%           2010-09-16 - Minor corrections
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

M = length(d);                        % Number of sensors
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
if ~isfield(params,'N_phi')
    N_phi = 64;                       % Default
else
    N_phi = params.N_phi;
end

w = 2*pi*fs*(0:NFFT/2)/NFFT;
phi = acos(2*(0:1/N_phi:1)-1);

% Calculate relative sensor distances
d_rel = d - d(1);

% Initialize waitbar
% h = waitbar(0,'Generating sensor signals...');

% Calculate sensor signals in the frequency domain
for phi_idx = 1:N_phi
%     waitbar(phi_idx/N_phi);
    X_prime = randn(1,NFFT/2+1) +1i*randn(1,NFFT/2+1);    
    X(1,:) = X(1,:) + X_prime;   
    for m = 2:M
        Delta = d_rel(m)*cos(phi(phi_idx));
        X(m,:) = X(m,:) + X_prime.*exp(-1i*Delta*w/c);
    end    
end
X = X/sqrt(N_phi);

% Transform to time domain
X = [sqrt(NFFT)*real(X(:,1)), sqrt(NFFT/2)*X(:,2:NFFT/2),...
    sqrt(NFFT)*real(X(:,NFFT/2+1)), sqrt(NFFT/2)*conj(X(:,NFFT/2:-1:2))];
z = real(ifft(X,NFFT,2));

% Truncate output signals
z = z(:,1:len);

% Close waitbar
% close(h);