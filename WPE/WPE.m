% Weighted Prediction Error for dereverberation
% reference: Nakatani, Tomohiro, et al. "Speech dereverberation based on
% variance-normalized delayed linear prediction." IEEE TASLP 18.7 (2010)
% ZitengWANG@201903


clear all
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
cfg.T60 = 0;                      % reverberation time (s)

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
cfg.dist = 3;

cfg.SNR = 20;
if cfg.SNR ~= inf               % inf means no noise
    cfg.noiseType = 'white';    % choice {'white' 'diffuse' 'recorded'}
    if strcmp(cfg.noiseType, 'recorded')
        % check first the noise is longer than speech!
        cfg.noiseFile = '';     
    end
end

cfg.SIR = inf;

% setup room and collect data
setup_room
setup_noise


%% offline processing start
% the following parts are specific to the algorithm
Nfft = 512;

Y = stft_multi_2(y, Nfft);
[Nframe, Nbin, Nch] = size(Y);

iterMax = 5;
wpe.delay = 3;
wpe.taps = 10;
XeEst = Y;
for bin = 1:Nbin
    % get audio context first
    YbarFrm = cell(Nframe, 1);
    for frm = (wpe.taps + wpe.delay):Nframe
        Ybar = squeeze(Y(frm - wpe.delay,bin,:));
        for tap=1:wpe.taps-1
            Ybar = [Ybar; squeeze(Y(frm - wpe.delay - tap,bin ,:))];
        end
        YbarFrm{frm} = Ybar;
    end
    
    % iterate
    X = squeeze(Y(:,bin,:));
    for iter = 1:iterMax
        % calculate mean signal power
        Xpow = mean(abs(X).^2, 2);
        XpowMax = max(Xpow);
        % calculate correlation matrix and correlation vector
        R = 0;
        P = 0;
        for frm = (wpe.taps + wpe.delay):Nframe
            Ybar = YbarFrm{frm};
            YbarTmp = Ybar / max(Xpow(frm), 1e-10*XpowMax);
            R = R + YbarTmp * Ybar' ;
            P = P + YbarTmp * squeeze(Y(frm,bin,:))';
        end
        hWPE = R \ P;
        
        % apply the filter
        for frm = (wpe.taps + wpe.delay):Nframe
            X(frm, :) = squeeze(Y(frm,bin,:)) - hWPE'*YbarFrm{frm};
        end
    end
    XeEst(:,bin,:) = permute(X,[1,3,2]);
end

xeEst = istft_multi_2(XeEst, length(speech));
audiowrite([saveDir prefix 'WPE' postfix '.wav'], xeEst, fs);




