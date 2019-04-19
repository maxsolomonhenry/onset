%   Testing various methods of onset detection
%   Max Henry
%   MUMT 501 Digital Audio Signal Processing

%   This script reviews a series of signal processing techniques for onset
%   detection in music signals, as outlined in Bello et al. (2005).  Like
%   the paper, it is divided into sections various stages of processing in
%   onset detection:
%
%       (1) This first section will load the four sample sounds, to be
%       picked and processed later on.  Please run this section before
%       running any other part of the script.  The example signal to be
%       processed later on is determined on line 49.
%
%       (2) Preprocessing -- an example of decomposition into multiple
%       bands from Sheirer (1998) is demonstrated.
%
%       (3) Signal reduction in two parts:
%
%       (a) Reduction based on signal features -- two envelope follower
%       formulas are demonstrated.
%
%       (b) Reduction based on spectral features -- operations on the
%       STFT: (i) HFC weighting by Masri (1996), (ii) spectral flux as
%       used by Duxbury (2002).
%
%       (4) Thresholding and peak-picking, following the procedures as
%       outlined in the review paper.  A final waveform is displayed, with
%       asterisks indicating locations of detected note-onsets.
%
%       (5) More cowbell! A simple audio engine that replays the sample
%       signal with an added cowbell track, coinciding with the detected
%       note onsets.

clear;

fs = 44100;
nyq = fs/2;
T = 1/fs;

%   load test signals
x1 = audioread('piano.wav');
x2 = audioread('guitar.wav');
x3 = audioread('violin.wav');
x4 = audioread('dabka.wav');

%   ---------- SET EXAMPLE ANALYSIS FILE HERE ----------

EXAMPLE = x4;
t = (0:length(EXAMPLE)-1)'/fs;      % time index in seconds

%   ----------------------------------------------------

%   Transient demo-plots
figure();
subplot(1,2,1);
plot((0:length(x1)-1)/fs, x1);
title('Note event');
xlabel('Time (s)');
xlim([1.5 1.9]);
ylabel('Amplitude');

subplot(1,2,2);
plot((0:length(x1)-1)/fs, x1);
title('Transient');
xlabel('Time (s)');
xlim([1.53 1.56]);
ylabel('Amplitude');

%   Unprocessed sound plots
figure();
subplot(2, 2, 1);
plot((0:length(x1)-1)/fs, x1);
title('piano.wav');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 2, 2);
plot((0:length(x2)-1)/fs, x2);
title('guitar.wav');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 2, 3);
plot((0:length(x3)-1)/fs, x3);
title('violin.wav');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 2, 4);
plot((0:length(x4)-1)/fs, x4);
title('dabka.wav');
xlabel('Time (s)');
ylabel('Amplitude');

black = [0 0 0];                % set colour for plot overlays later on

%%  (2) Pre-processing Examples

%   -------------------- Multiple bands -- Sheirer (1998) --------------------

prex = EXAMPLE;

[z, p, k] = ellip(6, 3, 40, [200 400]/nyq, 'bandpass');             % use zero-pole-gain form to avoid instability
bpsos = zp2sos(z,p,k);                                              % convert to sos value for sosfilt function
bpprex = sosfilt(bpsos, prex);                                             

%   bandpass plots
figure();
freqz(bpsos);

figure();
plot(t, bpprex);
title('Bandpassed signal - Sheirer (1998)');
xlabel('Time (s)');
ylabel('Amplitude');

%%  (3a) Reduction Examples (Temporal)

%   -------------------- "Envelope Follower," equation 1. --------------------

envx = EXAMPLE;                         % set example input

N = 501;                                % set averaging window size
hN = floor(N/2);
window = hann(N);

%   Faster method for envelope follower calculation:

temp = filter(hann(N), 1, abs(envx));   % calculate moving average of amplitude
temp1 = temp/max(temp);                 % normalize
follower = circshift(temp1, -hN);       % shift to compensate for phase delay



%   OLD Explicit method for envelope follower calculation: necessitates
%   padding and a slow for loop.
%
% envxpad = padarray(envx, hN, 'both');   % pad signal to allow for windowing at edges
% envN = length(envxpad);                 % set length N to padded signal
% 
% followerpad = zeros(envN, 1);           % initialize follower
% 
% for n = (hN+1:envN - hN)
%     windowed = abs(envxpad(n - hN:n + hN)).* window;    % set local window    
%     followerpad(n) = sum(windowed)/N;                   % sum window contents
% end
% 
% follower = followerpad(hN+1:end-hN);                % unpad
% follower = follower/max(follower);                  % normalize


%  Comparison plots for "envelope follower" eq 1

figure();
plot(t, envx);
hold on;
plot(t, follower, 'Color', black);
hold off;
title(['Rectifying envelope follower, N = ' num2str(N)]);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original signal', 'Rectified signal (normalized)');

%   -------------------- "Energy Envelope Follower," equation 2. --------------------

%   Faster method for Energy envelope calculation:

temp2 = filter(hann(N), 1, (envx.^2));                  % calculate moving average of signal energy
temp3 = temp2/max(temp2);                               % normalize
energyfollower = circshift(temp3, -hN);                % shift to compensate for phase delay

%   OLD Explicit method for Energy follower calculation: necessitates padding
%   and for loop.
%
% energyfollowerpad = zeros(envN, 1);                     % initialize follower
% 
% for n = (hN+1:envN - hN)
%     windowed = (envxpad(n - hN:n + hN)).^2.* window;    % square values and multiply by window
%     energyfollowerpad(n) = sum(windowed)/N;             % sum window contents
% end
% 
% energyfollower = energyfollowerpad(hN+1:end-hN);        % unpad
% energyfollower = energyfollower/max(energyfollower);    % normalize

%   Comparison plots for envelope follower eq 2
figure();
plot(t, envx);
hold on;
plot(t, energyfollower, 'Color', black);
hold off;
title(['Local energy envelope follower, N = ' num2str(N)]);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original signal', 'Local energy envelope (normalized)');

%%  (3b) Reduction Examples (Spectral)
                                 
specx = EXAMPLE;                                                    % set example input
specPreN = length(specx);                                           % get sample length (unpadded)

%   -------------------- Explicit STFT --------------------

%   WINDOW SIZE MUST BE ODD
N = 513;                                                            % set window size
window = hann(N);
h = floor(N*0.25);                                                  % hop size for 75% overlap
hN = floor(N/2);                                                    % lower half window size

specxpad = padarray(specx, hN, 'both');                             % pad signal to allow for windowing at edges
specN = length(specxpad);                                           % set N value to padded length

nframes = length(specx(1:h:end));                                   % calculate approx number of frames needed

%   Initializing values
STFT = zeros(nframes, N);
centerSTFT = zeros(nframes, N);
buffer = zeros(1, N);

framepointer = 1;                                                   % initialize frame counter

for n = (hN+1:h:specN - hN)
    windowed = specxpad((n - hN):(n + hN)).* window;                % set centered local window
    
    buffer = fft(windowed);
    
    STFT(framepointer,:) = buffer;
    centerSTFT(framepointer,:) = fftshift(buffer);                  % center frames around 0
    
    framepointer = framepointer + 1;                % increment frame
end

hSTFT = centerSTFT(:, (hN+1:end));                                  % take positive half of spectrum

%   Positive frequency STFT Plot
% STFTplot = flipud(abs(hSTFT)');                     % clean up for presentation
% figure();
% imagesc(STFTplot);
% colorbar;
% xlabel('Frame number');
% set(gca, 'YTick', [0:0.1:1]*hN, 'YTickLabel', [1:-0.1:0]*fs/2);
% ylabel('Frequency (Hz)');
% title('Standard STFT');

%   ---------- STFT using spectrogram() -----------------

% N = 512;
% h = floor(N/4);
% [STFT, tSTFT] = spectrogram(specx', N, N-h);
% nframes = size(STFT, 2);                        % get number of frames

%   ---------- HFC weighting -- Masri (1996) ----------

preHFC = zeros(nframes, N);

kindex = 0;    % indexer to avoid non-positive indexing values

%   first half of HFC calculations
for framepointer = (1:nframes)
    for k = (-hN:hN)
        kindex = k + hN + 1;
        preHFC(framepointer,kindex) = abs(k) * abs(centerSTFT(framepointer,kindex)).^2;
    end
end

HFC = (1/N) * sum(preHFC, 2);                          % finish calculations on HFC


HFC = HFC/max(HFC);                                     % normalize
interpHFC = interp1(1:nframes, HFC, (1:specPreN)/h);    % interpolate up to sample rate
sigHFC = interpHFC(1:specPreN);                         % truncate to original signal size
sigHFC = sigHFC';

%   HFC Plot 
figure();
plot(t, specx); hold on;
plot(t, sigHFC, 'Color', black); hold off;
title('HFC - Masri (1996)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original signal','HFC function (normalized)');

%  ---------- Spectral flux: L2 norm on rectified difference -- Duxbury (2002) ----------

preSflux = zeros(nframes, N);
pre = 0;

for framepointer = (1:nframes)
        %   avoid negative indexing (pre = 0 at first frame only)
        if (framepointer > 1)
            pre = 1;
        end
        
    for q = (1:N - 1) 
        %   find the difference between the absolute values of a frame and its preceding frame.
        preSflux(framepointer, q) = abs(centerSTFT(framepointer, q)) - ...
            pre * abs(centerSTFT(framepointer - pre, q));
    end
end

preSflux(preSflux<0) = 0;                                   % set negative values to 0.

Sflux = sum(preSflux.^2, 2);                                % calculate spectral flux
Sflux = Sflux/max(Sflux);                                   % normalize
interpSflux = interp1(1:nframes, Sflux, (1:nframes*h)/h);   % interpolate up to sample rate
sigSflux = interpSflux(1:specPreN);                         % truncate to original signal size
sigSflux = sigSflux';

%   Plot spectral flux
figure();
plot(t, specx);
hold on;
plot(t, sigSflux, 'Color', black);
hold off;
xlabel('Time (s)');
ylabel('Amplitude');
title('Spectral flux - Duxbury (2002)');
legend('Original signal', 'Spectral flux (normalized)');

%% Peak Picking

peakx = EXAMPLE;

%   Select an example detection function.  Variables are:
%   
%   follower            --  envelope follower 1
%   energyfollower      --  envelope follower 2 (energy follower)
%   sigSflux            --  spectral flux
%   sigHFC              --  HFC weighted spectrum
%
%   Specify selected function below:   

detectx = sigSflux;
detectN = length(detectx);

%   You may now want to adjust:
%
%   cutoff              --  the detection function filter frequency
%   delta               --  fixed threshold amount
%   lambda              --  contribution of adaptive thresholding (0 - 1)
%   M                   --  adjust the size of the median window

%   Normalize by subtracting the mean and dividing by max absolute deviation

mean = nanmean(detectx);
detectx = (detectx - mean);
detectx = detectx/max(detectx);
detectx(isnan(detectx)) = 0;                            % set NaN values to 0

%   Plot normalized signal
figure();
plot(t, detectx);
xlabel('Time (s)');
ylabel('Level');
title('Normalized detection function');

%   ---- Post-Processing: filtering the detection signal -- Bello et al. (2005) ----

cutoff = 10;                                % filter cutoff frequency in Hz

Wn = cutoff*2/fs;                           % cutoff value in normalized frequency
windowsize = 5000;
b = fir1(windowsize, Wn, 'low',...
    blackman(windowsize+1));

shift = windowsize/2;

temp = filter(b, 1, detectx);
temp2 = circshift(temp, -shift);            % shift result to account for filter phase delay
filtdetectx = [temp2(1:end - shift); ...    % replace final values with signal mean
    repelem(-mean, shift)'];

%   Plot pre and post filter
figure();
plot(t, detectx); hold on; 
plot(t, filtdetectx, 'Color', black); hold off;
xlabel('Time (s)');
ylabel('Level');
title('Filtered detection function');
legend('Detection function ', 'Filtered detection function');

%   Downsample the detection function for quick calculations

downsamp = 225;                                         % downsample factor
Tdown = downsamp/fs;                                    % calculate downsampled period

downdetectx = resample(filtdetectx, 1, downsamp);       % downsample detection function
downN = length(downdetectx);
tdownsamp = (1:downN)*Tdown;                            % time vector for downsampled functions

%   ---- Adaptive Thesholding: Moving-Median Filter -- Bello et al. (2005)  ----

threshold = zeros(1, downN);                % initialize moving threshold

%  threshold parameters
delta = 0.005;                              % fixed threshold value
lambda = 0.9;                               % amount to scale the median window
M = 100;                                    % ~half median window length

downdetectxpad = padarray(downdetectx, M, 'both');  % pad signal to allow for median at edges

for q = (M+1 : downN)
    threshold(q - M) = delta + lambda * nanmedian( abs(downdetectxpad(q-M:q+M)) );
end

%   Calculate detection function above threshold 
threshfunction = downdetectx - threshold';

%   Find local maxima above 0
[pks, locs] = findpeaks(threshfunction,'MinPeakHeight', 0);
locstime = locs*Tdown;

%   Plot moving threshold vs detection function, and "thresholded" function
figure();
subplot(2, 1, 1);
plot(tdownsamp, downdetectx); hold on; 
plot(tdownsamp, threshold); hold off;
xlabel('Time (s)');
ylabel('Level');
title('Adaptive threshold -- moving median filter');
legend('Detection function ', 'Threshold');

subplot(2, 1, 2);
plot(tdownsamp, threshfunction); hold on;
scatter(locstime, pks, '*'); hold off;
xlabel('Time (s)');
ylabel('Level');
ylim([0 max(threshfunction)]);
title('Detection function above threshold with local maxima');

%   Plot detected peaks and original signal
figure();
plot(t, peakx); hold on;
scatter(locstime, zeros(1, length(locs)), '*'); hold off;
xlabel('Time (s)');
ylabel('Level');
title('Original signal and detected note onsets');

%% More Cowbell!  Sound Demonstation

cow = audioread('cowbell.wav');
cowN = length(cow);
outN = length(EXAMPLE);

uplocs = locs*downsamp;                             % upsample locations of note onsets

cowtrack = zeros(max(uplocs) + cowN, 1);

%   Every time there is a peak, place a cowbell
for q = (1:length(uplocs))
    pointer = uplocs(q);
    cowtrack(pointer:pointer + cowN - 1) = cow;
end

cowtrack = [cowtrack; zeros(outN, 1)];              % overpad to ensure cowbell is longer than example

output = 0.5*cowtrack(1:outN) + EXAMPLE;            % combine cowbell and sample track

soundsc(output, fs);                                % enjoy the fruits of your labour
