%  the GSC filter in a two source scenario
% reference: Multichannel Eigenspace Beamforming in a Reverberant Noisy
% Environment With Multiple Interfering Speech Signals, 2009
% % Ziteng Wang @ 201812
% % Note: the smoothing parameter mu affects the performance

clear all
close all
warning('off')
addpath('..\STFT\')
addpath('..\Simulation\')
addpath('..\Simulation\RIR-Generator\')

%% simulation start
flatStart = 1;
prefix = '';   % for saving file

speechDir = '..\Simulation\Data\';
speechFile = 'fajw0_sa1.wav';
saveDir = 'GeneratedData2\';
if ~exist(saveDir)
    mkdir(saveDir)
end
speech = audioread([speechDir speechFile]);

% configuration
cfg = [];
cfg.fs = 16000;                     % sampling rate
cfg.room = [6 5 3];                 % room dimension (m)
cfg.T60 = 0.3;                      % reverberation time (s)

cfg.Nch = 6;
cfg.micCenter = [2 3 1.5];          % array center (m)
cfg.micCoordinate =[0.0425,0.0,0.0;
        0.02125,0.03680608,0.0;
        -0.02125,0.03680608,0.0;
        -0.0425,0.0,0.0;
        -0.02125,-0.03680608,0.0;
        0.02125,-0.03680608,0.0;];  % microphone array coordinates

cfg.az = 180;
cfg.el = 0;
cfg.dist = 2;

cfg.SNR = 30;
if cfg.SNR ~= inf                   % inf means no noise
    cfg.noiseType = 'diffuse';      % choice {'white' 'diffuse' 'recorded'}
    if strcmp(cfg.noiseType, 'recorded')
        % check first the noise is longer than speech!
        cfg.noiseFile = '';     
    end
end

cfg.SIR = 0;              
if cfg.SIR ~= inf
    % check the interference is longer than speech!
    cfg.interfFile =[speechDir '\fgjd0_si818.wav'];        
    cfg.azITF = 120;
    cfg.elITF = 0;
    cfg.distITF = 2;
end

% setup room and collect data
setup_room
setup_noise


%% processing start
% the following parts are specific to the algorithm
Nfft = 1024;
Nshift = Nfft / 2;
Nbin = Nfft / 2 + 1;
Nch = cfg.Nch;
win = sqrt(hanning(Nfft));

Nframe = floor(size(y,1)/Nshift) - 1;

% constraint set
C = zeros(Nch, 2, Nbin);

% % free field
for bin=1:Nbin
    C(:,1,bin) = exp(-2*1i*pi*(bin-1)/Nfft*fs*TDOA)/Nch;
    C(:,2,bin) = exp(-2*1i*pi*(bin-1)/Nfft*fs*TDOAinterf)/Nch;
end
postfix = [postfix '_freeC'];

% % ATF
% for ch=1:Nch
%     tmp = fft(rirSimu(ch,:), Nfft);
%     C(ch,1,:) = tmp(1:Nbin);
%     tmp = fft(rirITFSimu(ch,:), Nfft);
%     C(ch,2,:) = tmp(1:Nbin);
% end
% postfix = [postfix '_atfC'];


%% GSC lms
g = [0; 1];

q = zeros(Nch, Nbin);
pest = zeros(Nbin,1);
mu = 0.05;
alphaP = 0.9;
% resume
st = 1;
Y = zeros(Nch, Nbin);
Yout = zeros(Nbin,1);
yout = zeros(size(y,1), 1);

w0 = zeros(Nch, Nbin);
B = zeros(Nch, Nch, Nbin);
for bin = 1:Nbin
    tmp = C(:,:,bin)' * C(:,:,bin);
    if rcond(tmp) < eps
        tmp = tmp + 1e-10 * eye(size(tmp,1));
    end
    w0(:,bin) = C(:,:,bin) / tmp * g;
    B(:,:,bin) = eye(Nch) - C(:,:,bin) / tmp * C(:,:,bin)';
end
for frm=1:Nframe
    for ch=1:Nch
        tmp = fft(win .* y(st:st+Nfft-1,ch), Nfft);
        Y(ch,:) = tmp(1:Nbin);
    end
    for bin=1:Nbin  
        w = w0(:,bin) - B(:,:,bin)' * q(:,bin);
        Yout(bin) = w' * Y(:,bin);
        
        u = B(:,:,bin) * Y(:,bin);
        pest(bin) = alphaP * pest(bin) + (1-alphaP) * sum(abs(u).^2);
        q(:,bin) = q(:,bin) + mu * u * conj(Yout(bin)) / pest(bin);
    end
    yout(st:st+Nfft-1) = yout(st:st+Nfft-1) + win .* real(ifft([Yout; conj(Yout(end-1:-1:2))]));
    st = st + Nshift;
end

audiowrite([saveDir prefix 'GSC_lms' postfix '.wav'], yout, fs);

%% GSC wiener (closed form solution) is LCMV
% % resume
% st = 1;
% Y = zeros(Nch, Nbin);
% PhiN = zeros(Nch, Nch, Nbin);
% alpha = 0.98;
% Yout = zeros(Nbin,1);
% yout = zeros(size(y,1), 1);
% for frm=1:Nframe
%     for ch=1:Nch
%         tmp = fft(win .* y(st:st+Nfft-1,ch), Nfft);
%         Y(ch,:) = tmp(1:Nbin);
%         tmp = fft(win .* noise(st:st+Nfft-1,ch), Nfft);
%         N(ch,:) = tmp(1:Nbin);
%     end
%     for bin=1:Nbin
%         if frm == 1
%             PhiN(:,:,bin) = N(:,bin) * N(:,bin)';
%         else
%             PhiN(:,:,bin) = alpha * PhiN(:,:,bin) + (1-alpha) * N(:,bin) * N(:,bin)';
%         end
%         if frm > 5
%             tmp = B(:,:,bin) * PhiN(:,:,bin) * B(:,:,bin)';
%             if rcond(tmp) < eps;
%                 tmp = tmp + 1e-10 * eye(Nch);
%             end
%             q = tmp \ B(:,:,bin) * PhiN(:,:,bin) * w0(:,bin);
%             w = w0(:,bin) - B(:,:,bin)'* q;
%             Yout(bin) = w'* Y(:,bin);
%         end
%     end
%     yout(st:st+Nfft-1) = yout(st:st+Nfft-1) + win .* real(ifft([Yout; conj(Yout(end-1:-1:2))]));
%     st = st + Nshift;
% end
% 
% audiowrite([saveDir prefix 'GSC_wiener' postfix '.wav'], yout, fs);