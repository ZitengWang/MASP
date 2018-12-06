% use the estimated RTF to reconstruct the speech covariance matrix of 
% a single target to be a rank-1 matrix.
% Ziteng Wang @201812

function PhiX = rank1Approx(PhiX, PhiN)
[Nch,~,Nbin] = size(PhiX);
for bin = 1:Nbin
    if rcond(PhiN(:,:,bin)) < eps
        disp(['bin ' num2str(bin) ': Noise covariance ill-conditioned.']);
        PhiN(:,:,bin) = PhiN(:,:,bin) + 1e-10 * eye(Nch);
    end
    % GEVD (or Covariance Whitenning) based RTF estimation;
    % could also use the principle eigenvector (or Covariance Substraction)
    [vv,dd] = eig(PhiX(:,:,bin), PhiN(:,:,bin));
    [~,idx] = max(diag(dd));
    RTF = PhiN(:,:,bin) * vv(:,idx);
    tmp = RTF * RTF';
    % normalize
    PhiX(:,:,bin) = trace(PhiX(:,:,bin)) / trace(tmp) * tmp;
end