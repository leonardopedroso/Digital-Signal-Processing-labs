%% Lab 1 - Sampling and Aliasing - R4, R5, R6
%% 0. Init
clear;
%% 4
% Upload sound file
[x, Fs] = audioread('romanza_pe.wav');
Fs
%% play sampled sound
soundsc(x, Fs);
%% stop sound
clear sound
%% Choose window length
% Set window length
N= 4*[8 50 100 200 300 400];  
% Compute and display spectrogram of x 
figure();
for i = 1:6
    subplot(3,2,i);
    spectrogram(x(1:Fs*15), hann(N(i)), 3*N(i)/4, 4*N(i), Fs, 'yaxis');
end
%% Selected window length
N = 4*600;
figure();
spectrogram(x(1:Fs*15), hann(N), 3*N/4, 4*N, Fs, 'yaxis');
set(gca,'FontSize',35);
%% 5.
% Fs/5 sampling -> y
y = x(1:5:length(x));
Fs_y = Fs/5
%% play sampled sound
soundsc(y, Fs_y);
%% stop sound
clear sound
%% Plot sprectrogram
N = 4*600/5;
figure();
spectrogram(y(1:Fs_y*15), hann(N), 3*N/4, 4*N, Fs_y, 'yaxis');
set(gca,'FontSize',35);
%% 6
% 100th order Low-pass FIR filter
xf = filter(fir1(100, 0.2), 1, x);
% Cut-off freqeuncy (Hz)
fc = 0.2*pi/(2*pi*1/Fs)
N = 4*600;
figure();
spectrogram(xf(1:Fs*15), hann(N), 3*N/4, 4*N, Fs, 'yaxis');
set(gca,'FontSize',35);
%% Fs/5 sampling -> y
yf = xf(1:5:length(x));
Fs_yf = Fs/5;
%% play sampled sound
soundsc(yf, Fs_yf);
%% stop sound
clear sound
%%
N = 4*600/5;
figure();
spectrogram(yf(1:Fs_y*15), hann(N), 3*N/4, 4*N, Fs_y, 'yaxis');
set(gca,'FontSize',35);