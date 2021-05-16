%% Lab 02 - R1
%% Init
clear;

%% R1.a)
% Define constants
M = 512;
w_0 = 5.2*2*pi/M; % (rad)
% Create signal x
n = (0:M-1)';
x = 5*cos(w_0*n+1)+2*cos(2*w_0*n+2)+3*cos(5*w_0*n+3);

%% R1.b)
% Plot and save signal
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',35);
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)$",'Interpreter','latex');
plot(0:M-1,x,'LineWidth',3);
xlim([0,M-1]);
saveas(gcf,'R1b.png');

%% R1.c)
N = 512; % Set N
[absX,angX,wX] = dft_custom(x,N); % Perform Fast Fourier Tranform
% Plot |X(k)|
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$2\pi k/N \;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("$|X(k)|$",'Interpreter','latex');
set(gca,'FontSize',35);
plot(wX,absX,'LineWidth',3);
xlim([min(wX) max(wX)]);
saveas(gcf,'R1c_mag.png');
% Plot arg(X(k))
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$2\pi k/N\;[\mathrm{rad}]$",'Interpreter','latex');
ylabel("$\mathrm{arg}(X(k)) \;[\mathrm{rad}]$",'Interpreter','latex');
set(gca,'FontSize',35);
plot(wX,angX,'LineWidth',3);
xlim([min(wX) max(wX)]);
saveas(gcf,'R1c_arg.png');

%% R1.d)
% Find peaks of |X(k)|
[magpks,idxpks] = findpeaks(absX);
idxpks % indices of peaks + 1 (matlab indexes position as n as n+1 )
magpks % magnitude of peaks
wXpks = wX(idxpks) % normalized frequency of peaks
angpks = angX(idxpks) % phase of peaks

%% R1.e)
% Linear combination of the complex exponentials of each peak
xr = zeros(M,1);
for p = 1:3
    xr = xr + magpks(p)*exp(1i*angpks(p))*exp(1i*wXpks(p)*n);
end 
xr = real(xr);
% Plot linear combination of complex exponentials against original signal
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',35);
plot(0:M-1,x,'LineWidth',3);
plot(0:M-1,xr,'LineWidth',3);
xlim([0,M-1]);
legend({"$x(n)$","$x_r(n)$"},'Interpreter','latex');
saveas(gcf,'R1e.png');
% Compute fitness (sum of squared error)
SSE = sum((x-xr).^2)

%% R1.f)
N = 1024; % Set N
[absX,angX,wX] = dft_custom(x,N); % Perform Fast Fourier Tranform
% Find peaks of |X(k)|
[magpks_aux,idxpks_aux] = findpeaks(absX);
for p = 1:3
    [magpks(p),idx] = max(magpks_aux);
    magpks_aux(idx) = -inf;
    idxpks(p) = idxpks_aux(idx);
end

idxpks % indices of peaks
magpks % magnitude of peaks
wXpks = wX(idxpks) % normalized frequency of peaks
angpks = angX(idxpks) % phase of peaks

% Linear combination of the complex exponentials of each peak
xr = zeros(M,1);
for p = 1:3
    xr = xr + magpks(p)*exp(1i*angpks(p))*exp(1i*wXpks(p)*n);
end 
xr = real(xr);

% Plot linear combination of complex exponentials against original signal
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',35);
plot(0:M-1,x,'LineWidth',3);
plot(0:M-1,xr,'LineWidth',3);
xlim([0,M-1]);
legend({"$x(n)$","$x_r(n)$"},'Interpreter','latex');
saveas(gcf,'R1f.png');

% Compute fitness (sum of squared error)
SSE = sum((x-xr).^2)