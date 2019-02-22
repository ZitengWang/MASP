function c = projection_back(Y, ref)
% Parameters:
% Y : T x F x N
% ref : T x F
% c : F x N
num = squeeze(sum(bsxfun(@times, conj(ref), Y), 1));
denom = squeeze(sum(abs(Y).^2, 1));
I = (denom > 0);
c = ones(size(num));
c(I) = num(I) ./ denom(I);
