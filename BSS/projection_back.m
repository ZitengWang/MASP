function c = projection_back(Y, ref)
% This function computes the frequency-domain filter that minimizes
%     the squared error to a reference signal. This is commonly used
%     to solve the scale ambiguity in BSS.
%     Derivation of the projection
%     ----------------------------
%     The optimal filter `z` minimizes the squared error.
%     
%     .. math::
%         \min E[|z^* y - x|^2]
%     It should thus satsify the orthogonality condition
%     and can be derived as follows
% .. math::
%         0 & = E[y^*\\, (z^* y - x)]
%         0 & = z^*\\, E[|y|^2] - E[y^* x]
%         z^* & = \\frac{E[y^* x]}{E[|y|^2]}
%         z & = \\frac{E[y x^*]}{E[|y|^2]}
%     In practice, the expectations are replaced by the sample
%     mean.
% Parameters:
% Y : T x F x N
% ref : T x F
% c : F x N
num = squeeze(sum(bsxfun(@times, conj(ref), Y), 1));
denom = squeeze(sum(abs(Y).^2, 1));
I = (denom > 0);
c = ones(size(num));
c(I) = num(I) ./ denom(I);
