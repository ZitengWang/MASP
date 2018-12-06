%% room setup
% returns:  x,  (:,Nch) the scaled reverberant speech
%           TDOA, (Nch, 1)  the time difference of arrival

c = 343;                        % Sound velocity (m/s)
fs = cfg.fs;                    % Sample frequency (samples/s)
room = cfg.room;                % Room dimensions [x y z] (m)
T60 = cfg.T60;                  % Reverberation time (s)
beta = max(T60 * fs, 300);
SAVE_RVB_DETAILS = 0;           % save the early/late reverberant files

postfix = [postfix '_RT' num2str(T60)];

Nch = cfg.Nch;
micCenter = cfg.micCenter;          % array center (m)
micCoordinate= cfg.micCoordinate;   % microphone array coordinates
micPose = micCoordinate + repmat(micCenter, Nch, 1);

% place the source in the room
az = cfg.az;
el = cfg.el;
dist = cfg.dist;
sourcePose = dist * [cos(el/180*pi)*cos(az/180*pi) cos(el/180*pi)*sin(az/180*pi) sin(el/180*pi)] + micCenter;

if flatStart
    rirSimu = rir_generator(c, fs, micPose, sourcePose, room, T60, beta);
    save([saveDir 'rirSimu_RT' num2str(T60) '_tmp.mat'],'rirSimu');
else
    load([saveDir 'rirSimu_RT' num2str(T60) '_tmp.mat'],'rirSimu');
end

% observe
x = zeros(length(speech), Nch);
rirScale = 20;
for ch=1:Nch
    x(:,ch) = rirScale*fftfilt(rirSimu(ch,:), speech);
end

% early and late reverberation
TDOA = sqrt(sum((bsxfun(@minus, sourcePose, micPose)).^2, 2))/c;

if SAVE_RVB_DETAILS
    xd = zeros(length(speech), Nch);    % the direct sound
    xr = zeros(length(speech), Nch);    % the reverberant signal without the direct part 
    xe = zeros(length(speech), Nch);    % the early reverberation
    xl = zeros(length(speech), Nch);    % the late reverberation
    idxd = floor(max(TDOA)*fs) + 5;
    idxe = floor(max(TDOA)*fs) + 0.05*fs;  % threshold for early reverberation: 50 ms
    for ch=1:Nch
        xd(:,ch) = rirScale*fftfilt(rirSimu(ch,1:idxd), speech);
        xr(:,ch) = rirScale*fftfilt(rirSimu(ch,idxd+1:end), speech);
        xe(:,ch) = rirScale*fftfilt(rirSimu(ch,1:idxe), speech);
        xl(:,ch) = rirScale*fftfilt(rirSimu(ch,idxe+1:end), speech);
    end
    audiowrite([saveDir 'reverb' postfix '_direct.wav'], xd ,fs);
    audiowrite([saveDir 'reverb' postfix '_rvb.wav'], xr ,fs);
    audiowrite([saveDir 'reverb' postfix '_early.wav'], xe ,fs);
    audiowrite([saveDir 'reverb' postfix '_late.wav'], xl ,fs);
end

audiowrite([saveDir 'reverb' postfix '.wav'], x ,fs);
