
clear all
close all
addpath('..\STFT')

Nfft = 1024;
Nshift=Nfft/2;

% input data
[x, fs]=audioread('Data\reverb_RT0_interf_SIR0.wav'); 
x = x(:,1:2);

% fft
X = stft_multi_2(x, Nfft);
[Nframe,Nbin,Nch] = size(X);

%% AuxIVA
% demixing matrix, Nbin x Nch x Nch 
W = repmat(eye(Nch),1,1,Nbin); W = permute(complex(W,0), [3,1,2]);
I = eye(Nch);
% Y: Nframe x Nbin x Nch
Y = squeeze(sum(bsxfun(@times, X, permute(conj(W), [4,1,2,3])) ,3));

epochs = 30;
for ep=1:epochs
    % Nframe x Nch
    R = sqrt(squeeze(sum(abs(Y).^2, 2)));
    Gr = 1 ./ (R+eps);
    for bin = 1:Nbin
        for ch = 1:Nch
            % compute V
            Vtmp = bsxfun(@times, permute(X(:,bin,:),[1,3,2]), conj(X(:,bin,:)));
            V = squeeze(sum(bsxfun(@times, Gr(:,ch), Vtmp), 1)) / Nframe;
            % update W
            WV = squeeze(W(bin,:,:))' * V;
            Wtmp = WV \ I(:,ch);
            Wtmp = Wtmp / sqrt(Wtmp' * V * Wtmp);
            W(bin,:,ch) = Wtmp;
        end
    end
    %¡¡demixing
    Y = squeeze(sum(bsxfun(@times, X, permute(conj(W), [4,1,2,3])) ,3));
end

% project back
z = projection_back(Y, X(:,:,1));
Y = bsxfun(@times, Y, permute(z, [3,1,2]));


%% output
y = istft_multi_2(Y, size(x,1));
saveDir = 'Data\';
audiowrite([saveDir 'AuxIVA_src1.wav'], y(:,1), fs);
audiowrite([saveDir 'AuxIVA_src2.wav'], y(:,2), fs);