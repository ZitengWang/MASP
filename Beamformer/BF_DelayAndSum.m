% % multichannel delay and sum beamformer
% % Ziteng Wang @ 201811

clear all
warning('off')
addpath('..\STFT\')
addpath('..\Simulation\')
addpath('..\Simulation\RIR-Generator\')

%% simulation start
flatStart = 1;
prefix = '';   % for saving file

speechDir = '..\Simulation\Data\';
speechFile = 'fajw0_sa1.wav';
saveDir = 'GeneratedData\';
if ~exist(saveDir)
    mkdir(saveDir)
end
speech = audioread([speechDir speechFile]);

% configuration
cfg = [];
cfg.fs = 16000;                     % sampling rate
cfg.room = [6 5 3];                 % room dimension (m)
cfg.T60 = 0.6;                      % reverberation time (s)

cfg.Nch = 6;
cfg.micCenter = [2 3 1.5];          % array center (m)
cfg.micCoordinate =[0.0425,0.0,0.0;
        0.02125,0.03680608,0.0;
        -0.02125,0.03680608,0.0;
        -0.0425,0.0,0.0;
        -0.02125,-0.03680608,0.0;
        0.02125,-0.03680608,0.0;];  % microphone array coordinates

cfg.az = 45;
cfg.el = 0;
cfg.dist = 2;

cfg.SNR = 10;
if cfg.SNR ~= inf               % inf means no noise
    cfg.noiseType = 'white';    % choice {'white' 'diffuse' 'recorded'}
    if strcmp(cfg.noiseType, 'recorded')
        % check first the noise is longer than speech!
        cfg.noiseFile = '';     
    end
end

cfg.SIR = inf;              
if cfg.SIR ~= inf
    % check the interference is longer than speech!
    cfg.interfFile = '';    
    cfg.azITF = 180;
    cfg.elITF = 0;
    cfg.distITF = 2;
end

% setup room and collect data
setup_room
setup_noise


%% processing start
% the following parts are specific to the algorithm
Nfft = 1024;

Y = stft_multi_2(y, Nfft);
[Nframe, Nbin, Nch] = size(Y);

%% delay and sum 
% h = 1/Nch * dirctVec
hDS = zeros(Nch, Nbin);
for bin=1:Nbin
    hDS(:,bin) = 1/Nch * exp(-2*1i*pi*(bin-1)/Nfft*fs*TDOA);
end
% apply the filter
Xest = sum(bsxfun(@times, conj(permute(hDS, [3,2,1])), Y), 3);
xest = istft_multi_2(Xest, length(speech));


%% record and plot
audiowrite([saveDir prefix 'DelayAndSum' postfix '.wav'], xest, fs);

% plot the beampattern
dist = 2;
el = 0;
for az=1:360
    sourcePose = dist * [cos(el/180*pi)*cos(az/180*pi) cos(el/180*pi)*sin(az/180*pi) sin(el/180*pi)] + micCenter;
    TDOA = sqrt(sum((bsxfun(@minus, sourcePose, micPose)).^2, 2))/c;
    for bin=1:Nbin
        azVec = exp(-2*1i*pi*(bin-1)/Nfft*fs*TDOA);
        beamPattern(az,bin) = abs(sum(conj(hDS(:,bin)) .* azVec));
        beamPattern(az,bin) = max(20*log10(beamPattern(az,bin)), -30); % in dB
    end
end
figure; imagesc(beamPattern); axis xy; colorbar;
title(['BeamPattern of DS steered towards azimuth ' num2str(cfg.az)]);
xlabel('frequency index'); ylabel('azimuth angle');

% white noise gain
for bin=1:Nbin
    whiteNoiseGain(bin) = 10*log10(1/(hDS(:,bin)'*hDS(:,bin)));
end
figure; plot(whiteNoiseGain); ylim([0,10]);
title('WNG of DS');
xlabel('frequency index'); ylabel('WNG in dB');