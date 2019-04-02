% the GEV beamformer
% GEV:    h = eigenVec(PhiN^-1 * PhiX)
% input:    PhiN, (Nch, Nch, Nbin) the noise covariance matrix
%           PhiX, (Nch, Nch, Nbin) the speech covariance matrix
%           BANflag, performe BAN post processing or not
% output:   h, (Nch, Nbin)  the beamformer coefficients
% % Ziteng Wang @ 201812

function h = GEV(PhiX, PhiN, BANflag)
if nargin < 3
    BANflag = 1;        %%% apply BAN by default
end

[Nch, ~, Nbin] = size(PhiX);
h = zeros(Nch, Nbin);

for bin = 1:Nbin
    if rcond(PhiN(:,:,bin)) < eps
        disp(['bin ' num2str(bin) ': Noise covariance ill-conditioned.']);
        PhiN(:,:,bin) = PhiN(:,:,bin) + 1e-10 * eye(Nch);
    end
    [vv, dd] = eig(PhiX(:,:,bin), PhiN(:,:,bin));
    [~, idx] = max(diag(dd));
    h(:,bin) = vv(:,idx);
    if BANflag        %%% If not, the output signal is badly scaled.
        h(:,bin) = h(:,bin) * sqrt(trace(PhiN(:,:,bin))) / sqrt(Nch);
    end
end

