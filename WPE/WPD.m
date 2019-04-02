% Weighted Power minimization Distortionless response
% reference: Nakatani, Tomohiro, and Keisuke Kinoshita. 
% "A unified convolutional beamformer for simultaneous denoising and 
%  dereverberation." arXiv preprint arXiv:1812.08400 (2018).
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
cfg.T60 = 0.6;                      % reverberation time (s)

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
Xe = stft_multi_2(xe, Nfft);
Xl = stft_multi_2(xl, Nfft);
N = stft_multi_2(noise, Nfft);

% covariance matrix (Nch, Nch, Nbin)
PhiXe = mean(bsxfun(@times, permute(Xe,[3,4,2,1]), permute(conj(Xe), [4,3,2,1])), 4);
PhiXlN = mean(bsxfun(@times, permute(Xl+N,[3,4,2,1]), permute(conj(Xl+N), [4,3,2,1])), 4);

wpd.b = 4;  % prediction delay
% {12,10, 8, 6} prediction filter for 0-0.8kHz, 0.8-1.5kHz, 1.5-3kHz, 3-8kHz
wpd.Lw = zeros(Nbin, 1);
wpd.Lw(1 : floor(800/cfg.fs*Nfft + 1)) = 12;
wpd.Lw(floor(800/cfg.fs*Nfft + 1)+1 : floor(1500/cfg.fs*Nfft + 1)) = 10;
wpd.Lw(floor(1500/cfg.fs*Nfft + 1)+1 : floor(3000/cfg.fs*Nfft + 1)) = 8;
wpd.Lw(floor(3000/cfg.fs*Nfft + 1)+1 : end) = 6;


refMic = 1;
XeEst = Y(:,:,refMic);
for bin=1:Nbin
    % calculate covariance
    R = 0;
    for frm = max(wpd.Lw)+1:Nframe
        Ybar = squeeze(Y(frm, bin, :));
        for dly = wpd.b:wpd.Lw(bin)
            Ybar = [Ybar;squeeze(Y(frm - dly, bin , :))];
        end
        PowXe = mean(abs(squeeze(Xe(frm,bin,:))).^2);
        R = R + Ybar*Ybar' / PowXe;
    end
    % calculate steering vector
    [vv,dd] = eig(PhiXe(:,:,bin), PhiXlN(:,:,bin));
    [~,idx] = max(diag(dd));
    RTFvec = PhiXlN(:,:,bin) * vv(:,idx);
    RTF = RTFvec / (RTFvec(refMic));    %%% reference normalization 
    RTF = RTF / norm(RTF);              %%% unit normalization
    steerVec = [RTF; zeros(cfg.Nch*(wpd.Lw(bin)-wpd.b+1),1)];
    % calculate the WPD beamformer
    tmp = R \ steerVec;
    hWPD = tmp / (steerVec'*tmp);
    
    % apply beamformer
    for frm = max(wpd.Lw)+1:Nframe
        Ybar = squeeze(Y(frm, bin, :));
        for dly = wpd.b:wpd.Lw(bin)
            Ybar = [Ybar;squeeze(Y(frm - dly, bin , :))];
        end
        XeEst(frm, bin) = hWPD'*Ybar;
    end
end
% check output 
xeEst = istft_multi_2(XeEst, length(speech));
audiowrite([saveDir prefix 'WPD' postfix '.wav'], xeEst, fs);

%%% impressive results
%%% to apply WPE+MPDR for comparison

