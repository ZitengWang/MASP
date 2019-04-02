%% noise setup
% returns:  noise, if existed  (:,Nch) the background noise
%           interf, if existed (:,Nch) the interference
%           y = x + noise + interf, (:,Nch)  the observed signal

SNR = cfg.SNR;
SIR = cfg.SIR;

source = x;     %%% Generally, the reverberant speech is the target source.

%% simulate noise
if SNR ~= inf
    % noise collected from 'white' 'diffuse' 'recorded'
    noiseType = cfg.noiseType;
    
    if strcmp(noiseType, 'white')
        % % white noise
        noise = wgn(size(source,1), size(source,2), -10);
    end
    
    if strcmp(noiseType, 'diffuse')
        % % diffuse noise
        addpath('..\Simulation\INF-Generator')
        params.c = 343;
        params.fs = fs;
        params.N = 256;
        noise = sinf_3D(micCoordinate', size(source,1), params)';
    end
    
    if strcmp(noiseType, 'recorded')
        % % recorded multichannel noise !
        %rng(1)
        noiseFile = cfg.noiseFile;
        noise = audioread(noiseFile);
        fst = floor(rand * (size(noise,1) - size(source,1) - 10));
        noise = noise(fst+1:fst+size(source,1),:);
    end
    
    % scale
    noise = sqrt( sum(sum(source.^2))/sum(sum(noise.^2))/(10^(SNR/10)) )*noise;
    postfix = [postfix '_' noiseType '_SNR' num2str(SNR)];
end

%% simulate interference
if SIR ~= inf
    % % directional interference
    interfFile = cfg.interfFile;
    interfSource = audioread(interfFile);
    fst = floor(rand * (size(interfSource,1) - size(source,1) - 10));
    interfSource = interfSource(fst+1:fst+size(source,1)); 
    
    azITF = cfg.azITF;
    elITF = cfg.elITF;
    distITF = cfg.distITF;
    interfPose = micCenter + distITF * [cos(elITF/180*pi)*cos(azITF/180*pi) cos(elITF/180*pi)*sin(azITF/180*pi) sin(elITF/180*pi)];
    rirITFSimu = rir_generator(c, fs, micPose, interfPose, room, T60, beta);
    interf = zeros(length(speech), Nch);
    for ch=1:Nch
        interf(:,ch) = rirScale * fftfilt(rirITFSimu(ch,:), interfSource);
    end
    
    TDOAinterf = sqrt(sum((bsxfun(@minus, interfPose, micPose)).^2, 2)) / c;
    % scale
    interf = sqrt( sum(sum(source.^2)) / sum(sum(interf.^2))/(10^(SIR/10)) ) * interf;
    postfix = [postfix '_interf_SIR' num2str(SIR)];
end


%% collect data
y = x;

if SNR ~= inf
    % add noise
    y = y + noise; 
    audiowrite([saveDir 'noise' postfix '.wav'], noise ,fs);
    audiowrite([saveDir 'reverb' postfix '.wav'], y ,fs);
end
if SIR ~= inf
    % add interference
    y = y + interf;
    audiowrite([saveDir 'interf' postfix '.wav'], interf ,fs);
    audiowrite([saveDir 'reverb' postfix '.wav'], y ,fs);
end


