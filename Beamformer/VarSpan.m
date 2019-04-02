% the Variable Span beamformer
% h = sum( gevdVec * gevdVec^H / (mu + lambda)) * PhiX * refVec
% input:    PhiN, (Nch, Nch, Nbin) the noise covariance matrix
%           PhiX, (Nch, Nch, Nbin) the speech covariance matrix
%           mu, the speech distortion/noise reduction trade-off parameter
%           span, 1 <= span <= Nch, dimension of the assumed source space
%           refMic, the reference microphone, default = 1
% output:   h, (Nch, Nbin)  the beamformer coefficients
% % Ziteng Wang @ 201812

function h = VarSpan(PhiX, PhiN, mu, span, refMic)
if nargin < 3
    mu = 1;             %%% typical value {0, 1}
elseif nargin < 4
    span = 1;           %%% default 1 for single-target
elseif nargin < 5
    refMic = 1;
end

[Nch, ~, Nbin] = size(PhiX);
refVec = zeros(Nch, 1);
refVec(refMic) = 1;
h = zeros(Nch, Nbin);

for bin = 1:Nbin
    if rcond(PhiN(:,:,bin)) < eps
        disp(['bin ' num2str(bin) ': Noise covariance ill-conditioned.']);
        PhiN(:,:,bin) = PhiN(:,:,bin) + 1e-10 * eye(Nch);
    end
    [vv, dd] = eig(PhiX(:,:,bin), PhiN(:,:,bin));
    [~, idx] = sort(diag(dd), 'descend');
    tmp = 0;
    for i = 1:span
        tmp = tmp + vv(:,idx(i)) * vv(:,idx(i))' / (mu + dd(idx(i),idx(i)));
    end
    h(:,bin) = tmp * PhiX(:,:,bin) * refVec;
end
