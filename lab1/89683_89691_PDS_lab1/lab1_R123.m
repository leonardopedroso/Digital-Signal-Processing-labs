%% Lab 1 - Sampling and Aliasing - R1, R2, R3

clear
close all

%% R1
k2 = 1000;  % s^(-3)
k1 = 0;     % s^(-2)
F0 = 0;     % Hz
phi0 = 0;   % dimensionless
fs = 8000;  % Hz
t0 = 0;     % s
tf = 2;     % s

% Creates the array of instants of time of the samples
nT = t0:1/fs:tf;

% Represents the continuous-time signal
xc = @(t) cos(2*pi*(k2*t.^3/3 + k1*t.^2/2 + F0*t + phi0));

% Performs the sampling
x = xc(nT);

%%
soundsc(x,fs)

%% R2 a)
% Sets various window lengths
Ns = 4*[8 20 50 80 100 120];

% Computes and displays the spectrogram of x 
figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(Ns)
    subplot(3,2,i);
    spectrogram(x, hann(Ns(i)), 3*Ns(i)/4, 4*Ns(i), fs, 'yaxis');
    title("N/4 = " + num2str(Ns(i)/4))
end

%% R2 b)
% Sets the chosen window length
N = 4*80;
figure('units','normalized','outerposition',[0 0 1 1]);
spectrogram(x,hann(N), 3*N/4, 4*N, fs, 'yaxis')
set(gca,'FontSize',20)
title("Spectrogram of x(n) with N="+num2str(N))

if false
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,'figures/spectrogram_x','-dpdf','-r0')
end

%% R3 a)
% Creates y
y = x(1:2:length(x));

%%
soundsc(y,fs/2)

%% R3 b)
% Creates Ny so that the window duration is maintained
Ny = 0.5*N;
figure('units','normalized','outerposition',[0 0 1 1]);
spectrogram(y,hann(Ny), 3*Ny/4, 4*Ny, fs/2, 'yaxis')
set(gca,'FontSize',20)
title("Spectrogram of y(n) with Ny="+num2str(Ny))

if false
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,'figures/spectrogram_y','-dpdf','-r0')
end

%% R3 c)
fis = [1e3, 2e3, 3e3, 4e3];
fss = [20e3, 8e3, 4e3];

z = @(f,t) cos(2*pi*f*t);

% For every frequency fi
for i = 1:length(fis)
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    % For every sampling frequency
    for j = 1:length(fss)
        % Creates a subplot
        subplot(1,3,j)
        nT = 0:1/fss(j):3/fis(i);
        % where plots the signal function
        fplot(@(t) z(fis(i),t), [0, 3/fis(i)])
        hold on
        % and the sampled signal
        scatter(nT,z(fis(i),nT))
        set(gca,'FontSize',24)
        title("Sampled at "+num2str(fss(j)) + "Hz")
        hold off
    end
    
    if false
        set(gcf,'Units','Inches');
        pathFigPos = get(gcf,'Position');
        set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
            'PaperSize',[pathFigPos(3), pathFigPos(4)])
        print(gcf,"figures/plot_z_"+num2str(fis(i))+"_Hz",'-dpdf','-r0')
    end
end

%% x at a very low sampling frequency
c = x(1:4:length(x));

%%
soundsc(c,fs/4)

%% R3 b)
% Creates Ny so that the window duration is maintained
Nc = 44;
figure('units','normalized','outerposition',[0 0 1 1]);
spectrogram(c,hann(Nc), 3*Nc/4, 4*Nc, fs/4, 'yaxis')
set(gca,'FontSize',20)
title("Spectrogram of c(n) with Nc="+num2str(Nc))

if false
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,'figures/spectrogram_c','-dpdf','-r0')
end