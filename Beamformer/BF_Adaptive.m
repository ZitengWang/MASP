% % multichannel adaptive beamformer for single target enhancement
% % Ziteng Wang @ 201812


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
saveDir = 'GeneratedData\';
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

cfg.SNR = 10;
if cfg.SNR ~= inf                   % inf means no signal
    cfg.noiseType = 'diffuse';      % choice {'white' 'diffuse' 'recorded'}
    if strcmp(cfg.noiseType, 'recorded')
        % check first the noise is longer than speech!
        cfg.noiseFile = '';     
    end
end

cfg.SIR = inf;              
if cfg.SIR ~= inf
    % check the interference is longer than speech!
    cfg.interfFile = [speechDir '\mrgg0_si1829.wav'];        
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
X = stft_multi_2(x, Nfft);
N = stft_multi_2(noise, Nfft);

PhiX = mean(bsxfun(@times, permute(X,[3,4,2,1]), permute(conj(X), [4,3,2,1])), 4);
PhiN = mean(bsxfun(@times, permute(N,[3,4,2,1]), permute(conj(N), [4,3,2,1])), 4);

%% Adaptive beamformer
% collected from {'MVDR','MWF','GEV','VS'}
BFtype = categorical({'MVDR', 'MWF', 'GEV', 'VS'});  

% optionally perform rank-1 approximation
% PhiX = rank1Approx(PhiX, PhiN);

if iscategory(BFtype, 'MVDR')
    % the MVDR filter
    refMic = 1;
    choice = 'eigenVec';    %%% collected from 'eigenVec' 'gevd'
    hMVDR = MVDR(PhiX, PhiN, refMic, choice);
    
    % apply the filter
    Xest = sum(bsxfun(@times, conj(permute(hMVDR, [3,2,1])), Y), 3);
    xest = istft_multi_2(Xest, length(speech));
    audiowrite([saveDir prefix 'MVDR_' choice postfix '.wav'], xest, fs);
end  
if iscategory(BFtype, 'MWF')
    % the MWF filter
    mu = 1;
    refMic = 1;
    choice = 'r1MWF';   %%% collected from 'r1MWF' 'SDWMWF'
    hMWF = MWF(PhiX, PhiN, mu, refMic, choice);
    
    % apply the filter
    Xest = sum(bsxfun(@times, conj(permute(hMWF, [3,2,1])), Y), 3);
    xest = istft_multi_2(Xest, length(speech));
    audiowrite([saveDir prefix 'MWF_' choice postfix '.wav'], xest, fs);
end    
if iscategory(BFtype, 'GEV')
    % the GEV filter
    BANflag = 1;
    hGEV = GEV(PhiX, PhiN, BANflag);
    
    % apply the filter
    Xest = sum(bsxfun(@times, conj(permute(hGEV, [3,2,1])), Y), 3);
    xest = istft_multi_2(Xest, length(speech));
    audiowrite([saveDir prefix 'GEV' postfix '.wav'], xest, fs);
end    
if iscategory(BFtype, 'VS')
    % the Variable Span filter
    mu = 1;
    span = 1;
    refMic = 1;
    hVS = VarSpan(PhiX, PhiN, mu, span, refMic);
    
    % apply the filter
    Xest = sum(bsxfun(@times, conj(permute(hVS, [3,2,1])), Y), 3);
    xest = istft_multi_2(Xest, length(speech));
    audiowrite([saveDir prefix 'VS' postfix '.wav'], xest, fs);
end


%% record and plot
% The Signal Distortion / Noise Reduction metrics are affected by the scale
% of the output signal.
Xout = sum(bsxfun(@times, conj(permute(hMVDR, [3,2,1])), X), 3);
xout = istft_multi_2(Xout, length(speech));
Nout = sum(bsxfun(@times, conj(permute(hMVDR, [3,2,1])), N), 3);
nout = istft_multi_2(Nout, length(speech));
oSNR = 10 * log10( sum(xout.^2) / sum(nout.^2));


