% the MWF beamformer
% r1MWF:    h = PhiN^-1 * PhiX * refVec / (mu + trace(PhiN^-1 * PhiX)
% SDWMWF:   h = (PhiX + mu * PhiN)^-1 * PhiX * refVec
% input:    PhiN, (Nch, Nch, Nbin) the noise covariance matrix
%           PhiX, (Nch, Nch, Nbin) the speech covariance matrix
%           mu, the speech distortion/noise reduction trade-off parameter
%           refMic, the reference microphone, default = 1
% output:   h, (Nch, Nbin)  the beamformer coefficients
% % Ziteng Wang @ 201812

function h = MWF(PhiX, PhiN, mu, refMic, choice)
if nargin < 3
    mu = 1;                 % typical value {0, 1}
elseif nargin < 4
    refMic = 1;
elseif nargin < 5
    choice = 'r1MWF';       % collected from 'r1MWF' 'SDWMWF'
end

[Nch, ~, Nbin] = size(PhiX);
refVec = zeros(Nch, 1); 
refVec(refMic) = 1;
lambda = zeros(Nbin, 1);
h = zeros(Nch, Nbin);

for bin = 1:Nbin
    if rcond(PhiN(:,:,bin)) < eps
        disp(['bin ' num2str(bin) ': Noise covariance ill-conditioned.']);
        PhiN(:,:,bin) = PhiN(:,:,bin) + 1e-10 * eye(Nch);
    end
    if strcmp(choice, 'r1MWF')
        tmp = PhiN(:,:,bin) \ PhiX(:,:,bin);
        lambda(bin) = trace(tmp);
        h(:,bin) = tmp * refVec / (mu + lambda(bin));
        %%% Gain factor to normalize the steering vector. The results are 
        %%% better. Needs further checking
        %gain_factor = norm(PhiX(:,refMic,bin) / PhiX(refMic,refMic,bin));
        %h(:,bin) = gain_factor * h(:,bin);
    elseif strcmp(choice, 'SDWMWF') % SDW-MWF is more sensitive to estimation errors.
        tmp = (PhiX(:,:,bin) + mu * PhiN(:,:,bin)) \ PhiX(:,:,bin);
        h(:,bin) = tmp * refVec;
    end 
end


