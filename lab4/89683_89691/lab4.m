%% Lab 4 - Prediction
%% Init
clear
close all
clc
if ~isfolder("figures")
    mkdir("figures")
end
%% R1.b)
% Load datasert
load("energy_train.mat",'x_train');
N = 96;
% Lest squares fit
[a_est,E_res] = ARmodelR1a(x_train,N);
% Compute long term prediction
x_pred = zeros(length(x_train)-N,1);
for i=N+1:length(x_train)
    x_pred(i-N) = a_est*x_train(i-N);
end
% Compute train residuals
res = x_train(N+1:end)-x_pred;
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',25);
scatter(1:length(x_train),x_train,'LineWidth',2)
scatter(N+1:length(x_train),x_pred,'LineWidth',2)
legend("Train","Prediction")
hold off
saveas(gcf,'./figures/R1b_train.png');
% Compute residual
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)-\hat{x}(n)$",'Interpreter','latex');
set(gca,'FontSize',25);
scatter(N+1:length(x_train),res,'LineWidth',2)
hold off
saveas(gcf,'./figures/R1b_error.png');
%% R1.c)
a_est
E_res
% Root mean quared error
RMS = sqrt(E_res/(length(x_train)-N))
%% R1.e)
P = 6; % Horizon
% Lest squares short-term fit
[a_hat_r,E_err] = ARmodelR1d(res,P);
% Compute short term prediction
res_pred = zeros(length(res)-P,1);
for i = P+1:length(res)
    for p = 1:P
        res_pred(i-P) = res_pred(i-P) + a_hat_r(p)*res(i-p);
    end
end
% Compute shor-term prediction error 
err = res_pred - res(P+1:end);
% Compute predicted dtat with both models
x_pred2 = x_pred(P+1:end)+res_pred;

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',25);
scatter(N+1:length(x_train),res,'LineWidth',2)
scatter(N+P+1:length(x_train),res_pred,'LineWidth',2)
legend("Train","Prediction")
hold off
saveas(gcf,'./figures/R1e_residual.png');
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',25);
scatter(1:length(x_train),x_train,'LineWidth',2)
scatter(N+P+1:length(x_train),x_pred2,'LineWidth',2)
legend("Train","Prediction")
hold off
saveas(gcf,'./figures/R1e_pred.png');

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
ylabel("$\hat{r}(n)-r(n)$",'Interpreter','latex');
set(gca,'FontSize',25);
scatter(N+P+1:length(x_train),err,'LineWidth',2)
hold off
saveas(gcf,'./figures/R1e_error.png');
%% R1.f
a_hat_r
E_err
RMS = sqrt(E_err/(length(x_train)-N-P))
%% R2. Anomaly detection
%% R2.a)
threshold = 0.1; % Set anomaly detection to 0.1
% Anomalty detection is performed by function 
% [anomaly] = anomalyR2a(x,x_pred,thrshl)
%% R2.b)
% Convert variables of the model for the second part
a_hat_x = a_est; 
load("energy_test.mat",'x_test'); % Load test data
x_pred = zeros(length(x_test),1); % Allocate predicted data
for i = N+1:length(x_test) % Predict x(n)
    x_pred(i) = x_test(i-N)*a_hat_x;
end
r_test = x_test(N+1:end)-x_pred(N+1:end); % Allocate test residues
r_pred = zeros(length(r_test),1); % Allocate predicted residues
for i = P+1:length(r_pred) % Predict residues
    for p = 1:P
        r_pred(i) = r_pred(i) + a_hat_r(p)*r_test(i-p);
    end
end

% Plot 
if true
% Plot prdicted x(n)-r(n) (long term prediction)
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(1:length(x_test),x_test,'LineWidth',3);
scatter(N+1:length(x_test),x_pred(N+1:end),'LineWidth',3);
xlim([0 length(x_test)]);
legend("Test data","Long-term prediction");
saveas(gcf,'./figures/R2b_test.png');
% Plot predicted residue
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(N+1:length(x_test),r_test,'LineWidth',3);
scatter(N+P+1:length(x_test),r_pred(P+1:end),'LineWidth',3);
xlim([0 length(x_test)]);
legend("Residue test","Residue prediction");
saveas(gcf,'./figures/R2b_residual.png');
end
if false
% Plot residue error
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(1:length(x_test),x_test,'LineWidth',3);
scatter(N+P+1:length(x_test),x_pred(N+P+1:end)+r_pred(P+1:end),'LineWidth',3);
xlim([0 length(x_test)]);
legend("Test data","Data + residual prediction");
saveas(gcf,'./figures/R2b_data_residual.png');
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
ylabel("$\hat{r}(n)-r(n)$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(N+P+1:length(x_test),r_pred(P+1:end)-r_test(P+1:end),'LineWidth',3);
xlim([0 length(x_test)]);
saveas(gcf,'./figures/R2b_error.png');
end
% Dectect anomalous periods
anomaly = anomalyR2a(r_test(P+1:end),r_pred(P+1:end),threshold);
% Masks for all time series with anomalies
n_an = (P+1:length(r_test))';
n_an = n_an(anomaly==1);
x_test_an = x_test(N+P+1:end);
x_test_an = x_test_an(anomaly==1);
r_test_an = r_test(P+1:length(r_test));
r_test_an = r_test_an(anomaly==1);
x_test_res_an = x_test(N+P+1:end);
x_test_res_an = x_test_res_an(anomaly==1);
e_an = x_pred(N+P+1:end)+r_pred(P+1:end)-x_test(N+P+1:end);
e_an = e_an(anomaly==1);
if true
% Plot long term prediction with error 
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(1:length(x_test),x_test,'LineWidth',3);
scatter(N+1:length(x_test),x_pred(N+1:end),'LineWidth',3);
scatter(N+n_an,x_test_an,120,'LineWidth',4)
xlim([0 length(x_test)]);
legend("Test data","Long-term prediction","Anomaly");
saveas(gcf,'./figures/R2b_test_an.png');
% Plot residue prediction with error
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(N+1:length(x_test),r_test,'LineWidth',3);
scatter(N+P+1:length(x_test),r_pred(P+1:end),'LineWidth',3);
scatter(N+n_an,r_test_an,120,'LineWidth',4)
xlim([0 length(x_test)]);
legend("Residual test","Residue prediction","Anomaly");
saveas(gcf,'./figures/R2b_residual_an.png');
if false
% Plot predicted x + residue with identification of anomalies  
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(1:length(x_test),x_test,'LineWidth',3);
scatter(N+P+1:length(x_test),x_pred(N+P+1:end)+r_pred(P+1:end),'LineWidth',3);
scatter(N+n_an,x_test_res_an,120,'LineWidth',4)
xlim([0 length(x_test)]);
legend("Test data","Data + residual prediction","Anomaly");
saveas(gcf,'./figures/R2b_data_residual_an.png');
end
% Plot residue error to detect anomalies with threshold line
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
ylabel("$\hat{r}(n)-r(n)$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(N+P+1:length(x_test),r_pred(P+1:end)-r_test(P+1:end),'LineWidth',3);
scatter(N+n_an,e_an,120,'LineWidth',4)
plot([0 max(N+P+1:length(x_test))], [threshold threshold],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
plot([0 max(N+P+1:length(x_test))], -[threshold threshold],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
xlim([0 length(x_test)]);
legend("Error","Anomaly");
saveas(gcf,'./figures/R2b_error_an.png');
end

%% R2.c) Correct predictions affected by anomalies  - Method 1
% ---------- Step 1 ----------
anomaly_cor = anomaly; % Find data that is considered to be contaminated 
i = N+P+1;
% Identify the days of the first anomaly
while i<= length(x_test)
    if anomaly_cor(i-N-P)==1
        n = floor(i/N);
        % Flag data contamination in the whole that when an anomaliy is
        % foud
        for j = N*(n)+1:N*(n+1)
            if j > length(x_test)
                break;
            end
            anomaly_cor(j-P-N) = 1;
        end
        % The anomalies of the day a fter an anomaly is flagged are
        % ignored
        for j = N*(n+1)+1:N*(n+2)
            if j > length(x_test)
                break;
            end
            anomaly_cor(j-P-N) = 0;
        end
        i = N*(n+2);
    else
        i = i+1;
    end
end
% Create masks to point out the data that is considered to be contaminated
n_an_cor = (P+1:length(r_test))';
n_an_cor = n_an_cor(anomaly_cor==1);
e_an_cor = x_pred(N+P+1:end)+r_pred(P+1:end)-x_test(N+P+1:end);
e_an_cor = e_an_cor(anomaly_cor==1);
% Plot data contaminated
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
ylabel("$\hat{r}(n)-r(n)$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(N+P+1:length(x_test),r_pred(P+1:end)-r_test(P+1:end),'LineWidth',3);
scatter(N+n_an_cor,e_an_cor,120,'LineWidth',4);
plot([0 max(N+P+1:length(x_test))], [threshold threshold],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
plot([0 max(N+P+1:length(x_test))], -[threshold threshold],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
xlim([0 length(x_test)]);
legend("Error","Data contaminated");
saveas(gcf,'./figures/R2c_error_an.png');
% ---------- Step 2 ----------
% Correct anomalies
x_test_cor = x_test;
x_test_cor([zeros(N,1);anomaly_cor] == 1) = x_pred([zeros(N,1);anomaly_cor] == 1)+r_pred([zeros(N,1);anomaly_cor] == 1);
% Predict - update data iterations until convergence
x_pred_cor = zeros(length(x_test_cor),1);
while true % Iterate until convergence
    % Predict long term
    for i = N+1:length(x_test_cor)
        x_pred_cor(i) = x_test_cor(i-N)*a_hat_x;
    end
    r_test_cor = x_test_cor(N+1:end)-x_pred_cor(N+1:end);
    r_pred_cor = zeros(length(r_test_cor),1);
     % Predict residues
    for i = P+1:length(r_pred_cor)
        for p = 1:P
            r_pred_cor(i) = r_pred_cor(i) + a_hat_r(p)*r_test_cor(i-p);
        end
    end
    % Check convergence
    if norm((x_test_cor([zeros(N,1);anomaly_cor] == 1)-(x_pred_cor([zeros(N,1);anomaly_cor] == 1)+r_pred_cor([zeros(N,1);anomaly_cor] == 1))))/norm(x_test_cor([zeros(N,1);anomaly_cor] == 1)) < 1e-3
        break;
    end
    % Update data
    x_test_cor([zeros(N,1);anomaly_cor] == 1) = x_pred_cor([zeros(N,1);anomaly_cor] == 1)+r_pred_cor([zeros(N,1);anomaly_cor] == 1);
end
% Final corrected test dataset
x_cor = x_pred_cor(N+P+1:end)+r_pred_cor(P+1:end);
e_cor = x_cor-x_test(N+P+1:end);

% ---------- Step 3 ----------
% Detect anomalies without contamination 
anomaly = anomalyR2a(e_cor,zeros(length(e_cor),1),threshold);
% Create masks to point out the anomalies after contaminated data is
% ignored
n_an = (P+1:length(r_test))';
n_an = n_an(anomaly==1);
x_test_an = x_test(N+P+1:end);
x_test_an = x_test_an(anomaly==1);
e_an = x_pred_cor(N+P+1:end)+r_pred_cor(P+1:end)-x_test(N+P+1:end);
e_an = e_an(anomaly==1);

% Plot
if false
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(1:length(x_test),x_test,'LineWidth',3);
scatter(N+P+1:length(x_test),x_cor(1:end),'LineWidth',3);
scatter(N+n_an,x_test_an,120,'LineWidth',4)
xlim([0 length(x_test)]);
legend("Test data","Data + residual prediction","Anomalous");
saveas(gcf,'./figures/R2d_data_residual.png');
end

% Plot residue error and point out the flagged anomalies
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
ylabel("$\hat{r}(n)-r(n)$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(N+P+1:length(x_test),e_cor,'LineWidth',3);
scatter(N+n_an,e_an,120,'LineWidth',4)
plot([0 max(N+P+1:length(x_test))], [threshold threshold],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
plot([0 max(N+P+1:length(x_test))], -[threshold threshold],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
xlim([0 length(x_test)]);
legend("Error","Anomaly");
saveas(gcf,'./figures/R2c_error_an_after.png');


%% R2.c) Correct predictions affected by anomalies  - Method 2
% ---------- Step 1 ---------- 
anomaly_cor = anomaly; % Find data that is considered to be contaminated 
i = N+P+1;
x_test_cor = x_test; % Create corrected data set
% Identify the days of the first anomaly
while i<= length(x_test_cor)
    if anomaly_cor(i-N-P)==1
        n = floor(i/N);
        % Flag data contamination in the whole that when an anomaliy is
        % found
        for j = N*(n)+1:N*(n+1)
            if j > length(x_test)
                break;
            end
            anomaly_cor(j-P-N) = 1;
            % Copy day of previous window multiplied by the factor a
            x_test_cor(j) = x_test(j-N)*a_hat_x; 
        end
        % The anomalies of the day a fter an anomaly is flagged are
        % ignored
        for j = N*(n+1)+1:N*(n+2)
            if j > length(x_test)
                break;
            end
            anomaly_cor(j-P-N) = 0;
        end
        i = N*(n+2);
    else
        i = i+1;
    end
end
% Create masks to point out the data that is considered to be contaminated
n_an_cor = (P+1:length(r_test))';
n_an_cor = n_an_cor(anomaly_cor==1);
e_an_cor = x_pred(N+P+1:end)+r_pred(P+1:end)-x_test(N+P+1:end);
e_an_cor = e_an_cor(anomaly_cor==1);
% Plot data contaminated
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
ylabel("$\hat{r}(n)-r(n)$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(N+P+1:length(x_test),r_pred(P+1:end)-r_test(P+1:end),'LineWidth',3);
scatter(N+n_an_cor,e_an_cor,120,'LineWidth',4);
plot([0 max(N+P+1:length(x_test))], [threshold threshold],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
plot([0 max(N+P+1:length(x_test))], -[threshold threshold],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
xlim([0 length(x_test)]);
legend("Error","Data contaminated");
saveas(gcf,'./figures/R2c_error_an_m2.png');

% ---------- Step 2 ----------
% Predict
x_pred_cor = zeros(length(x_test_cor),1);
% Predict long term
for i = N+1:length(x_test_cor)
    x_pred_cor(i) = x_test_cor(i-N)*a_hat_x;
end
r_test_cor = x_test_cor(N+1:end)-x_pred_cor(N+1:end);
r_pred_cor = zeros(length(r_test_cor),1);
 % Predict residues
for i = P+1:length(r_pred_cor)
    for p = 1:P
        r_pred_cor(i) = r_pred_cor(i) + a_hat_r(p)*r_test_cor(i-p);
    end
end
% Final corrected test dataset
e_cor = r_test_cor(P+1:end)-r_pred_cor(P+1:end);

% ---------- Step 3 ----------
% Detect anomalies without contamination 
anomaly = anomalyR2a(e_cor,zeros(length(e_cor),1),threshold);
% Create masks to point out the anomalies after contaminated data is
% ignored
n_an = (P+1:length(r_test))';
n_an = n_an(anomaly==1);
x_test_an = x_test(N+P+1:end);
x_test_an = x_test_an(anomaly==1);
e_an = x_pred_cor(N+P+1:end)+r_pred_cor(P+1:end)-x_test(N+P+1:end);
e_an = e_an(anomaly==1);

% Plot residue error and point out the flagged anomalies
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
ylabel("$\hat{r}(n)-r(n)$",'Interpreter','latex');
set(gca,'FontSize',35);
scatter(N+P+1:length(x_test),e_cor,'LineWidth',3);
scatter(N+n_an,e_an,120,'LineWidth',4)
plot([0 max(N+P+1:length(x_test))], [threshold threshold],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
plot([0 max(N+P+1:length(x_test))], -[threshold threshold],'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3);
xlim([0 length(x_test)]);
legend("Error","Anomaly");
saveas(gcf,'./figures/R2c_error_an_after_m2.png');