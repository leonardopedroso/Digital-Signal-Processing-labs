%% Lab 3 R1 and R2
%% Init
clear;
%% R1.a) Original sound file
% Load sound file
[signal, fs] = audioread('fugee.wav');
% Remove the last sample just to have an even N for the DFT
signal = signal(1:end-1); 
soundsc(signal,fs); % Play the signal
%% Stop playing
clear sound;
%% R1.b) Time domain analysis of original signal
% Plot Nsamples samples startint at n = 50
ti = 50;
Nsamples = 10000; 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(ti:ti+Nsamples-1,signal(ti:ti+Nsamples-1),'LineWidth',3);
xlim([ti ti+Nsamples-1]);
saveas(gcf,'./figures/R1b.png');
%% R1.c) Frequency domain analysis of original signal
% Compute DFT
[absX,~,wX] = dft_custom(signal,length(signal));
% Plot magnitude spectrum
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$2\pi k/N \;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("Magnitude",'Interpreter','latex');
set(gca,'FontSize',35);
plot(wX,absX,'LineWidth',3);
set(gca,'Yscale','log');
xlim([min(wX) max(wX)]);
saveas(gcf,'./figures/R1c.png');
% Plot spectrogram
N = 80;
figure('units','normalized','outerposition',[0 0 1 1]);
spectrogram(signal(ti:ti+Nsamples-1), hann(N), 3*N/4, 4*N, fs, 'yaxis');
set(gca,'FontSize',35);
saveas(gcf,'./figures/R1c_spectrogram.png');
hold off;
%% R2.a) Butterworth low-pass filter
butterOrder = 10; % filter order
wco = pi/2; % cytoff frequency
fco =  wco/(2*pi*1/fs); % continuous cutoff frequency
Wn = 2*fco/fs; % Normalised cutoff frequency
[bLTI,aLTI] = butter(butterOrder,2*fco/fs,'low'); % Compute filter
[h,w] = freqz(bLTI,aLTI,1e3); % Get magnitude and phase of filter
% Plot Gain of the filter 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$\omega \;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("Gain",'Interpreter','latex');
set(gca,'FontSize',35);
plot(w,abs(h),'LineWidth',3);
plot([wco wco],[0 1.2],'--','LineWidth',3);
%set(gca,'Yscale','log');
xlim([min(wX) max(wX)]);
saveas(gcf,'./figures/R2a_gain.png');
% Plot phase of the filter
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$\omega \;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("Phase",'Interpreter','latex');
set(gca,'FontSize',35);
plot(w,angle(h),'LineWidth',3);
xlim([min(wX) max(wX)]);
saveas(gcf,'./figures/R2a_phase.png');
%% R2.b) Filtered signal
signalLTIfilter = filter(bLTI,aLTI,signal); % filter the signal
fco % cutoff frequency
%% R2.c) Analysis of the filtered signal
ti = 50;
Nsamples = 10000; 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(ti:ti+Nsamples-1,signal(ti:ti+Nsamples-1),'LineWidth',3);
plot(ti:ti+Nsamples-1,signalLTIfilter(ti:ti+Nsamples-1),'LineWidth',3);
xlim([ti ti+Nsamples-1]);
legend('Original','Filtered');
saveas(gcf,'./figures/R2c.png');
ti = 1100;
Nsamples = 50; 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(ti:ti+Nsamples-1,signal(ti:ti+Nsamples-1),'LineWidth',3);
plot(ti:ti+Nsamples-1,signalLTIfilter(ti:ti+Nsamples-1),'LineWidth',3);
xlim([ti ti+Nsamples-1]);
legend('Original','Filtered');
saveas(gcf,'./figures/R2c_smallAmp.png');
ti = 1100;
Nsamples = 500; 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(ti:ti+Nsamples-1,signal(ti:ti+Nsamples-1),'LineWidth',3);
plot(ti:ti+Nsamples-1,signalLTIfilter(ti:ti+Nsamples-1),'LineWidth',3);
xlim([ti ti+Nsamples-1]);
legend('Original','Filtered');
saveas(gcf,'./figures/R2c_zoomNoise.png');
ti = 2000;
Nsamples = 500; 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(ti:ti+Nsamples-1,signal(ti:ti+Nsamples-1),'LineWidth',3);
plot(ti:ti+Nsamples-1,signalLTIfilter(ti:ti+Nsamples-1),'LineWidth',3);
xlim([ti ti+Nsamples-1]);
legend('Original','Filtered');
saveas(gcf,'./figures/R2c_zoomNormal.png');
%% R2.d)
% Compute DFT of filtered signal
[absXLTIfilter,~,wXLTIfilter] = dft_custom(signalLTIfilter,length(signalLTIfilter));
% Plot magnitude sprectra
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$2\pi k/N \;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("Magnitude",'Interpreter','latex');
set(gca,'FontSize',35);
plot(wX,absX,'LineWidth',3);
plot(wXLTIfilter,absXLTIfilter,'LineWidth',3);
plot([wco wco],[1e-12 1e-2],'--','LineWidth',3);
set(gca,'Yscale','log');
xlim([min(wX) max(wX)]);
legend('Original','Filtered');
saveas(gcf,'./figures/R2d.png');
% Plot spectrogram
ti = 50;
Nsamples = 10000;
N = 80;
figure('units','normalized','outerposition',[0 0 1 1]);
spectrogram(signalLTIfilter(ti:ti+Nsamples-1), hann(N), 3*N/4, 4*N, fs, 'yaxis');
set(gca,'FontSize',35);
saveas(gcf,'./figures/R2d_spectrogram.png');
hold off;
%% R2.e) Listen to the filtered signal
soundsc(signalLTIfilter,fs);
%% Stop sound
clear sound;
%% R2.f) Improve Butterworth low-pass filter
butterOrder = 10;
wco = pi/4;
fco =  wco/(2*pi*1/fs);
Wn = 2*fco/fs; % Normalised cutoff frequency
[bLTI,aLTI] = butter(butterOrder,2*fco/fs,'low'); % Compute filters
[h,w] = freqz(bLTI,aLTI,1e3); 
% Plot gain
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$\omega \;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("Gain",'Interpreter','latex');
set(gca,'FontSize',35);
plot(w,abs(h),'LineWidth',3);
%set(gca,'Yscale','log');
xlim([min(wX) max(wX)]);
saveas(gcf,'./figures/R2f_filter.png');
signalLTIfilter = filter(bLTI,aLTI,signal);
fco
ti = 50;
Nsamples = 10000; 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(ti:ti+Nsamples-1,signal(ti:ti+Nsamples-1),'LineWidth',3);
plot(ti:ti+Nsamples-1,signalLTIfilter(ti:ti+Nsamples-1),'LineWidth',3);
xlim([ti ti+Nsamples-1]);
legend('Original','Filtered');
saveas(gcf,'./figures/R2f_zoomOut.png');
ti = 1100;
Nsamples = 50; 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(ti:ti+Nsamples-1,signal(ti:ti+Nsamples-1),'LineWidth',3);
plot(ti:ti+Nsamples-1,signalLTIfilter(ti:ti+Nsamples-1),'LineWidth',3);
xlim([ti ti+Nsamples-1]);
legend('Original','Filtered');
saveas(gcf,'./figures/R2f_smallAmp.png');
ti = 1100;
Nsamples = 500; 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(ti:ti+Nsamples-1,signal(ti:ti+Nsamples-1),'LineWidth',3);
plot(ti:ti+Nsamples-1,signalLTIfilter(ti:ti+Nsamples-1),'LineWidth',3);
xlim([ti ti+Nsamples-1]);
legend('Original','Filtered');
saveas(gcf,'./figures/R2f_zoomNoise.png');
ti = 2000;
Nsamples = 500; 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(ti:ti+Nsamples-1,signal(ti:ti+Nsamples-1),'LineWidth',3);
plot(ti:ti+Nsamples-1,signalLTIfilter(ti:ti+Nsamples-1),'LineWidth',3);
xlim([ti ti+Nsamples-1]);
legend('Original','Filtered');
saveas(gcf,'./figures/R2f_zoomNormal.png');
[absXLTIfilter,~,wXLTIfilter] = dft_custom(signalLTIfilter,length(signalLTIfilter));
% Plot |X(k)|
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$2\pi k/N \;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("Magnitude",'Interpreter','latex');
set(gca,'FontSize',35);
plot(wX,absX,'LineWidth',3);
plot(wXLTIfilter,absXLTIfilter,'LineWidth',3);
plot([wco wco],[1e-12 1e-2],'--','LineWidth',3);
set(gca,'Yscale','log');
xlim([min(wX) max(wX)]);
legend('Original','Filtered');
saveas(gcf,'./figures/R2f_spectra.png');
% Plot spectrogram
ti = 50;
Nsamples = 10000;
N = 80;
figure('units','normalized','outerposition',[0 0 1 1]);
spectrogram(signalLTIfilter(ti:ti+Nsamples-1), hann(N), 3*N/4, 4*N, fs, 'yaxis');
set(gca,'FontSize',35);
saveas(gcf,'./figures/R2f_spectrogram.png');
hold off;
%% Play filtered sound
soundsc(signalLTIfilter,fs);
%% Stop sound
clear sound;
