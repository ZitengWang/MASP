% Geometric Source Separation
% reference: Parra, Lucas C., and Christopher V. Alvino.
% "Geometric source separation: Merging convolutive source separation 
%  with geometric beamforming." IEEE TASLP 10.6 (2002): 352-362.
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
saveDir = 'GeneratedData2\';
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
cfg.micCoordinate = [ 0.0463, 0.0, 0.0;
         0.02315,0.040096976, 0.0;
         -0.02315, 0.040096976, 0.0;
         -0.0463, 0.0,  0.0;
         -0.02315,-0.040096976, 0.0;
         0.02315,-0.040096976, 0.0]; % microphone array coordinates
    
cfg.az = 180;
cfg.el = 0;
cfg.dist = 2;

cfg.SNR = 30;
if cfg.SNR ~= inf               % inf means no noise
    cfg.noiseType = 'white';    % choice {'white' 'diffuse' 'recorded'}
    if strcmp(cfg.noiseType, 'recorded')
        % check first the noise is longer than speech!
        cfg.noiseFile = '';     
    end
end

cfg.SIR = 0;              
if cfg.SIR ~= inf
    % check the interference is longer than speech!
    cfg.interfFile = [speechDir '\mrgg0_si1829.wav'];         
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

Y = stft_multi_2(y, Nfft);
[Nframe, Nbin, Nch] = size(Y);

%%% assume the source number and source DOAs are known
Nsrc = 2;

% for speech target
steerVec = zeros(Nch,Nbin);
for bin=1:Nbin
    steerVec(:,bin) = exp(-2*1i*pi*(bin-1)/Nfft*fs*TDOA);
end
% for interference
steerVecITF = zeros(Nch,Nbin);
for bin=1:Nbin
    steerVecITF(:,bin) = exp(-2*1i*pi*(bin-1)/Nfft*fs*TDOAinterf);
end

%%% try to follow the ODAS codes
%%% The results are not as good as expected. Maybe I get it wrong.
lambda = 0.5;
mu = 0.01;
Wold = cell(Nbin,1);
Xest = zeros(Nframe,Nbin,Nsrc);
XestDS = zeros(Nframe,Nbin,Nsrc);
for frm=1:Nframe
    for bin=1:Nbin
        % If the source was previously active, load the corresponding demixing terms
        % Otherwise instantiate with terms of a weighted delay and sum beamformer
        if frm == 1
            Wnew(1,:) = 1/Nch * steerVec(:,bin)';
            Wnew(2,:) = 1/Nch * steerVecITF(:,bin)';
        else %%% may depend on the source activity
            Wnew = Wold{bin};
        end
        
        Yfrmbin = squeeze(Y(frm,bin,:));
        % Compute X = Wnew * Y
        X = Wnew * Yfrmbin;
        % Compute Rxx ?
        Rxx = X * X';
        % Compute E = Rxx - diag(Rxx)
        E = Rxx - diag(diag(Rxx));
        % Compute dJ = 4 * alpha * E * Wnew * Ryy = 4 * alpha * E * X * Y' ?
        alpha = (norm(Yfrmbin)^2)^2;
        dJ = 4 * E * X * Yfrmbin' / alpha;
        % Compute dJ2 = 2 * (Wnew * D - I) * D'
        D = [steerVec(:,bin) steerVecITF(:,bin)];
        dJ2 = 2 * (Wnew * D - eye(Nsrc)) * D';
        
        % Compute Wnew = (1 - lambda * mu) * Wnew - mu * (dJ + dJ2) ?
        Wnew = (1 - lambda * mu) * Wnew - mu * (dJ + dJ2);
      
        Xest(frm,bin,:) = Wnew * Yfrmbin;
        XestDS(frm,bin,1) = 1/Nch * steerVec(:,bin)' * Yfrmbin;
        XestDS(frm,bin,2) = 1/Nch * steerVecITF(:,bin)' * Yfrmbin;
        
        % Copy to demixing            
        Wold{bin} = Wnew;
    end
end
% output
x1est = istft_multi_2(Xest(:,:,1), length(speech));        
x2est = istft_multi_2(Xest(:,:,2), length(speech));  
audiowrite([saveDir prefix 'ODAS_src1' postfix '.wav'], x1est, fs);
audiowrite([saveDir prefix 'ODAS_src2' postfix '.wav'], x2est, fs);

x1est = istft_multi_2(XestDS(:,:,1), length(speech));        
x2est = istft_multi_2(XestDS(:,:,2), length(speech));  
audiowrite([saveDir prefix 'DS_src1' postfix '.wav'], x1est, fs);
audiowrite([saveDir prefix 'DS_src2' postfix '.wav'], x2est, fs);


%%% try to follow the paper
%%% The results are even worse. Needs further checking
% smoothFactor = 0.95;
% RyyOld = cell(Nbin,1);
% Wold = cell(Nbin,1);
% Xest = zeros(Nframe,Nbin,Nsrc);
% for frm=1:Nframe
%     for bin=1:Nbin
%         Yfrmbin = squeeze(Y(frm,bin,:));
%         
%         % If the source was previously active, load the corresponding demixing terms
%         % Otherwise instantiate with terms of a weighted delay and sum beamformer
%         if frm == 1
%             Wnew(1,:) = 1/Nch * steerVec(:,bin)';
%             Wnew(2,:) = 1/Nch * steerVecITF(:,bin)';
%             Ryy = Yfrmbin * Yfrmbin';
%         else %%% may depend on the source activity
%             Wnew = Wold{bin};
%             Ryy = smoothFactor * RyyOld{bin} + (1 - smoothFactor) * (Yfrmbin * Yfrmbin');
%         end
%         RyyOld{bin} = Ryy;
%         
%         % Compute Rxx = Wnew * Ryy * Wnew';
%         Rxx = Wnew * Ryy * Wnew';
%         % Compute E = Rxx - diag(Rxx)
%         E = Rxx - diag(diag(Rxx));
%         % Compute dJ = 4 * alpha * E * Wnew * Ryy
%         alpha = trace(Ryy * Ryy');
%         dJ = 4 * E * Wnew * Ryy / alpha;
%         % Compute dJ2 = 2 * (Wnew * D - I) * D'
%         D = [steerVec(:,bin) steerVecITF(:,bin)];
%         dJ2 = 2 * (Wnew * D - eye(Nsrc)) * D';
%         
%         % Compute Wnew = Wnew - mu * (dJ + lambda0 * lambda * dJ2)
%         lambda0 = 1/Nfft;
%         lambda = 1 / cond(D);
%         Wnew = Wnew - mu * (dJ + lambda0 * lambda * dJ2);
%         
%         Xest(frm,bin,:) = Wnew * Yfrmbin;
%         
%         % Copy to demixing            
%         Wold{bin} = Wnew;
%     end
% end
% % output
% x1est = istft_multi_2(Xest(:,:,1), length(speech));        
% x2est = istft_multi_2(Xest(:,:,2), length(speech));  
% audiowrite([saveDir prefix 'GSS_src1' postfix '.wav'], x1est, fs);
% audiowrite([saveDir prefix 'GSS_src2' postfix '.wav'], x2est, fs);

