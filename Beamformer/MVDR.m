% the MVDR beamformer
% h = PhiN^-1 * steerVec / steerVec^H * PhiN^-1 * steerVec
% input:    PhiN, (Nch, Nch, Nbin) the noise covariance matrix
%           PhiX, (Nch, Nch, Nbin) the speech covariance matrix
%           refMic, the reference microphone, default = 1
% output:   h, (Nch, Nbin)  the beamformer coefficients
% % Ziteng Wang @ 201812
% Note: MPDR is more sensitive to misalignment errors than MVDR

function h = MVDR(PhiX, PhiN, refMic, choice)
if nargin < 3
    refMic = 1;
elseif nargin < 4
    choice = 'eigenVec';    % collected from 'eigenVec' 'gevd'
end

[Nch, ~, Nbin] = size(PhiX);
steerVec = zeros(Nch, Nbin);
PhiNinv = zeros(Nch, Nch, Nbin);
for bin = 1:Nbin
    if rcond(PhiN(:,:,bin)) < eps
        disp(['bin ' num2str(bin) ': Noise covariance ill-conditioned.']);
        PhiN(:,:,bin) = PhiN(:,:,bin) + 1e-10 * eye(Nch);
    end
    PhiNinv(:,:,bin) = inv(PhiN(:,:,bin));
    steerVec(:,bin) =  get_steering_vector(PhiX(:,:,bin), PhiN(:,:,bin), refMic, choice);
end
tmp = squeeze(sum(bsxfun(@times, PhiNinv, permute(steerVec, [3,1,2])), 2));
h = bsxfun(@rdivide, tmp, sum(bsxfun(@times, conj(steerVec), tmp),1));


% calculate the steering vector (or the Relative Transfer Function)
% reference: Effect of Steering Vector Estimation on MVDR Beamformer for
% Noisy Speech Recognition, DSP 2018
function steerVec = get_steering_vector(PhiXbin, PhiNbin, refMic, choice)
if strcmp(choice, 'eigenVec')
    [vv, dd] = eig(PhiXbin);
    [~, idx] = max(diag(dd));
    steerVec = vv(:,idx);
    %%% By default, this steering vector is unit normalized. This however changes 
    %%% the scale of the output signal.
elseif strcmp(choice, 'gevd')
    [vv, dd] = eig(PhiXbin, PhiNbin);
    [~, idx] = max(diag(dd));
    steerVec = PhiNbin * vv(:,idx);
    %%% By default, this steering vector is normalized to the reference channel.
    steerVec = steerVec / (steerVec(refMic));
    %%% But unit normalization brings better result. Needs further checking
    % steerVec = steerVec / norm(steerVec);
end
