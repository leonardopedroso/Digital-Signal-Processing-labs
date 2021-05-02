% lab3.m
clear
close all

mkdir("./figures")
saveFigs = true;

[fugee,fs] = audioread("fugee.wav");
fugee = fugee(1:end-1); % so that the length is even

%%
soundsc(fugee,fs);

%%
clear sound

%% Median filter filtering
orderMedian = 3;
fugee_median = medfilt1(fugee,orderMedian);

%% filtered response zoom out
n = 1:1e4;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(n,fugee(n),'LineWidth',3)
plot(n,fugee_median(n),'LineWidth',3)
hold off
legend('Original','Median');

if saveFigs
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,"figures/comp_median_zoom_out",'-dpdf','-r0')
end

%% filtered response zoom in
n = 1100:1150;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(n,fugee(n),'LineWidth',3)
plot(n,fugee_median(n),'LineWidth',3)
hold off
legend('Original','Median');

if saveFigs
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,"figures/comp_median_zoom_in",'-dpdf','-r0')
end

%% filtered response zoom med
n = 2000:2500;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(n,fugee(n),'LineWidth',3)
plot(n,fugee_median(n),'LineWidth',3)
hold off
legend('Original','Median');

if saveFigs
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,"figures/comp_median_zoom_med",'-dpdf','-r0')
end

%% filtered response frequency
[absX,angX,f] = dft_custom(fugee,length(fugee));
[absX_median,angX_median,f_median] = ...
    dft_custom(fugee_median,length(fugee_median));

%%
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
set(gca,'FontSize',35);
plot(f,absX,'LineWidth',3)
plot(f_median,absX_median,'LineWidth',3)
set(gca,'Yscale','log');
grid on;
xlabel("$2\pi k/N \;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("Magnitude",'Interpreter','latex');
hold off
legend("Original","Median")

if saveFigs
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,"figures/dft_comp_median",'-dpdf','-r0')
end

%% Spectrogram
n = 50:50+1e4-1;
N = 80;
figure('units','normalized','outerposition',[0 0 1 1]);
spectrogram(fugee_median(n), hann(N), 3*N/4, 4*N, fs, 'yaxis');
set(gca,'FontSize',35);

if saveFigs
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,"figures/spectrogram_median",'-dpdf','-r0')
end

%%
soundsc(fugee_median,fs)

%%
clear sound

%% Median filter filtering
orderMedian = 7;
fugee_median7 = medfilt1(fugee,orderMedian);

%% filtered response zoom out
n = 1:1e4;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(n,fugee(n),'LineWidth',3)
plot(n,fugee_median7(n),'LineWidth',3)
hold off
legend('Original','Median');

if saveFigs
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,"figures/comp_median7_zoom_out",'-dpdf','-r0')
end

%% filtered response zoom in
n = 1100:1150;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(n,fugee(n),'LineWidth',3)
plot(n,fugee_median7(n),'LineWidth',3)
hold off
legend('Original','Median');

if saveFigs
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,"figures/comp_median7_zoom_in",'-dpdf','-r0')
end

%% filtered response zoom med
n = 2000:2500;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(n,fugee(n),'LineWidth',3)
plot(n,fugee_median7(n),'LineWidth',3)
hold off
legend('Original','Median');

if saveFigs
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,"figures/comp_median7_zoom_med",'-dpdf','-r0')
end

%% filtered response frequency
[absX,angX,f] = dft_custom(fugee,length(fugee));
[absX_median,angX_median,f_median] = ...
    dft_custom(fugee_median7,length(fugee_median7));

%%
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
set(gca,'FontSize',35);
plot(f,absX,'LineWidth',3)
plot(f_median,absX_median,'LineWidth',3)
set(gca,'Yscale','log');
grid on;
xlabel("$2\pi k/N \;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("Magnitude",'Interpreter','latex');
hold off
legend("Original","Median")

if saveFigs
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,"figures/dft_comp_median7",'-dpdf','-r0')
end

%% Spectrogram
n = 50:50+1e4-1;
N = 80;
figure('units','normalized','outerposition',[0 0 1 1]);
spectrogram(fugee_median7(n), hann(N), 3*N/4, 4*N, fs, 'yaxis');
set(gca,'FontSize',35);

if saveFigs
    set(gcf,'Units','Inches');
    pathFigPos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pathFigPos(3), pathFigPos(4)])
    print(gcf,"figures/spectrogram_median7",'-dpdf','-r0')
end

%%
soundsc(fugee_median7,fs)

%%
clear sound