%  the LCMV filter in a two source scenario
% reference: Multichannel Eigenspace Beamforming in a Reverberant Noisy
% Environment With Multiple Interfering Speech Signals, 2009
% % Ziteng Wang @ 201812
% % Note: by setting g=[1; 0], LCMV attenuates the second source to some
% extent. 
% When using MVDR, the interferece can be well cancelled if its statitics 
% is correctly estimated.

clear all
close all
warning('off')
addpath('..\STFT\')
addpath('..\Simulation\RIR-Generator\')

%% simulation start
flatStart = 1;
postfix = '';   % for saving file

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
    cfg.interfFile = [speechDir '\fgjd0_si818.wav'];         
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

%%% free field
% for bin=1:Nbin
%     C(:,1,bin) = exp(-2*1i*pi*(bin-1)/Nfft*fs*TDOA)/Nch;
%     C(:,2,bin) = exp(-2*1i*pi*(bin-1)/Nfft*fs*TDOAinterf)/Nch;
% end
% postfix = [postfix '_freeC1'];

%%% ATF
for ch=1:Nch
    tmp = fft(rirSimu(ch,:), Nfft);
    C(ch,1,:) = tmp(1:Nbin);
    tmp = fft(rirITFSimu(ch,:), Nfft);
    C(ch,2,:) = tmp(1:Nbin);
end
postfix = [postfix '_atfC'];

%% LCMV
g = [0; 1];

% start processing

st = 1;
Y = zeros(Nch, Nbin);
PhiY = zeros(Nch, Nch, Nbin);
N = zeros(Nch, Nbin);
PhiN = zeros(Nch, Nch, Nbin);
alpha = 0.98;
Yout = zeros(Nbin,1);
yout = zeros(size(y,1), 1);
for frm=1:Nframe
    for ch=1:Nch
        tmp = fft(win.*y(st:st+Nfft-1,ch), Nfft);
        Y(ch,:) = tmp(1:Nbin);
        tmp = fft(win .* noise(st:st+Nfft-1,ch), Nfft);
        N(ch,:) = tmp(1:Nbin);
    end
    for bin=1:Nbin
        if frm == 1
            PhiN(:,:,bin) = N(:,bin) * N(:,bin)';
        else
            PhiN(:,:,bin) = alpha * PhiN(:,:,bin) + (1-alpha) * N(:,bin) * N(:,bin)';
        end
        if frm > 5
            wtmp = PhiN(:,:,bin) \ C(:,:,bin);
            wtmp2 = C(:,:,bin)' * wtmp;
            if rcond(wtmp2) < eps
                wtmp2 = wtmp2 + 1e-10*eye(2);
            end
            w = wtmp / (wtmp2) * g;
            Yout(bin) = w' * Y(:,bin);
        end
    end
    yout(st:st+Nfft-1) = yout(st:st+Nfft-1) + win .* real(ifft([Yout; conj(Yout(end-1:-1:2))]));
    st = st + Nshift;
end

audiowrite([saveDir 'LCMV' postfix '.wav'], yout / max(abs(yout)) * 0.8, cfg.fs);

