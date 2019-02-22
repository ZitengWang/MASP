c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [3 3 1.5; 3 4 1.5; 3 5 1.5; 3 6 1.5;];              % Receiver position [x y z] (m)
s = [3 1 2];              % Source position [x y z] (m)
L = [6 8 5];                % Room dimensions [x y z] (m)
beta = 0.4;                 % Reflections Coefficients
n = beta*fs;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                  % Reflection order
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

h = rir_generator(c, fs, r, s, L, beta, n, mtype, order, dim, orientation, hp_filter);

sig = audioread('E:\Data\TIMIT\\test_wav_fmt\\fcmh1_sx233.wav');
for ch=1:size(h,1)
    out(:,ch) = 20*fftfilt(h(ch,:),sig);
end
