% lab4_1.m
if ~isfolder("figures")
    mkdir("figures")
end

%%
clear
close all
clc

load energy_train.mat

N = 96;
[a_est,E_res] = ARmodelR1a(x_train,N);
x_pred = zeros(length(x_train)-N,1);
for i=N+1:length(x_train)
    x_pred(i-N) = a_est*x_train(i-N);
end

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
saveFigAsPDF(gcf,"figures/long_term_pred")

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
ylabel("$x(n)-\hat{x}(n)$",'Interpreter','latex');
set(gca,'FontSize',25);
scatter(N+1:length(x_train),res,'LineWidth',2)
hold off
saveFigAsPDF(gcf,"figures/residuals")

a_est
E_res

%%
P = 6;
[a_hat_r,E_err] = ARmodelR1d(res,P);

res_pred = zeros(length(res)-P,1);
for i = P+1:length(res)
    for p = 1:P
        res_pred(i-P) = res_pred(i-P) + a_hat_r(p)*res(i-p);
    end
end

err = res_pred - res(P+1:end);

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
saveFigAsPDF(gcf,"figures/short_term_pred_res")

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
set(gca,'FontSize',25);
scatter(1:length(x_train),x_train,'LineWidth',2)
scatter(N+P+1:length(x_train),x_pred2,'LineWidth',2)
legend("Train","Prediction")
hold off
saveFigAsPDF(gcf,"figures/short_term_pred_x")

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$n$",'Interpreter','latex');
ylabel("$\hat{r}(n)-r(n)$",'Interpreter','latex');
set(gca,'FontSize',25);
scatter(N+P+1:length(x_train),err,'LineWidth',2)
hold off
saveFigAsPDF(gcf,"figures/errors")

a_hat_r
E_err