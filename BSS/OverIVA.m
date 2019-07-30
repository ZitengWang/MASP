%OverIVA, reference:
% Robin Scheibler and Nobutaka Ono "Independent Vector Analysis with more
% microphones than sources"
% ZitengWANG@201907

clear all
close all
addpath('..\STFT')

Nfft = 1024;
Nshift=Nfft/2;

% input data
[x, fs]=audioread('GeneratedData\reverb_RT0.3_white_SNR50_interf_SIR0.wav'); 
x = x(:,1:4);

% fft
X = stft_multi_2(x, Nfft);
[Nframe,Nbin,Nch] = size(X);
Nsrc = 2;

%% OverIVA
% W_hat: Nbin x Nch x Nch 
W_hat = repmat(zeros(Nch, Nch),1,1,Nbin); W_hat = permute(complex(W_hat,0), [3,1,2]);
I = eye(Nch);
% init W_hat
for bin=1:Nbin
    W_hat(bin,1:Nsrc,1:Nsrc) = permute(eye(Nsrc),[3,1,2]);
    W_hat(bin,Nsrc+1:Nch,Nsrc+1:Nch) = -permute(eye(Nch-Nsrc), [3,1,2]);
end
% demix Y: Nframe x Nbin x Nsrc
W = W_hat(:,:,1:Nsrc);
Y = squeeze(sum(bsxfun(@times, X, permute(conj(W), [4,1,2,3])) ,3));
% Covariance: Nbin x Nch x Nch
C = squeeze(mean(bsxfun(@times, X, permute(conj(X), [1,2,4,3])), 1));

epochs = 30;
for ep=1:epochs
    % Laplacian: Nframe x Nsrc
    R = sqrt(squeeze(sum(abs(Y).^2, 2)));
    Gr = 1 ./ (R+0.01);
    for bin = 1:Nbin
        for src = 1:Nsrc
            % compute V
            Vtmp = bsxfun(@times, permute(X(:,bin,:),[1,3,2]), conj(X(:,bin,:)));
            V = squeeze(sum(bsxfun(@times, Gr(:,src), Vtmp), 1)) / Nframe;
            % update W
            WV = squeeze(W_hat(bin,:,:))' * V;
            Wtmp = WV \ I(:,src);
            Wtmp = Wtmp / sqrt(Wtmp' * V * Wtmp);
            W_hat(bin,:,src) = Wtmp;
            % update J 
            CW = squeeze(C(bin,:,:))*squeeze(W_hat(bin,:,1:Nsrc));
            J = CW(Nsrc+1:Nch,:)/CW(1:Nsrc,:);
            W_hat(bin,1:Nsrc,Nsrc+1:Nch) = J';
        end
    end
    %¡¡demixing
    W = W_hat(:,:,1:Nsrc);
    Y = squeeze(sum(bsxfun(@times, X, permute(conj(W), [4,1,2,3])) ,3));
end

% project back
z = projection_back(Y, X(:,:,1));
Y = bsxfun(@times, Y, permute(z, [3,1,2]));


%% output
y = istft_multi_2(Y, size(x,1));
saveDir = 'GeneratedData\';
audiowrite([saveDir 'OverIVA_src1.wav'], y(:,1), fs);
audiowrite([saveDir 'OverIVA_src2.wav'], y(:,2), fs);
