% lab5.m
if ~exist('./figures','dir')
    mkdir('figures')
end

%% R1.b) Image visualization
clear
close all
clc

% Load
load("sar_image.mat",'I')
% Visualize
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$x$",'Interpreter','latex');
ylabel("$y$",'Interpreter','latex');
set(gca,'FontSize',25);
xlim([0 size(I,2)]);
ylim([0 size(I,1)]);
imagesc(I)
saveas(gcf,'./figures/R1b.png');

%% R1.c) Image cropping of water and ice sections
% rectangles
rectIce = [1 80 85 403];
rectWater = [270 1 300 320];

I_ice = imcrop(I,rectIce);
x_ice = I_ice(:);
I_water = imcrop(I,rectWater);
x_water = I_water(:);

% Visualize image with superimposed rectangles
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$x$",'Interpreter','latex');
ylabel("$y$",'Interpreter','latex');
set(gca,'FontSize',25);
xlim([0 size(I,2)]);
ylim([0 size(I,1)]);
imagesc(I)
rectangle('Position',rectIce,'EdgeColor','r','LineWidth',3);
rectangle('Position',rectWater,'EdgeColor','r','LineWidth',3);
saveas(gcf,'./figures/R1c.png');

%% R1.c) Parameter estimation
close all

% normal distribution
[mu_water,var_water] = MLEnormal(x_water)
[mu_ice,var_ice] = MLEnormal(x_ice)

% exponential distribution
lambda_water = MLEexponential(x_water)
lambda_ice = MLEexponential(x_ice)

% rayleigh distribution
sigma_water = MLErayleigh(x_water)
sigma_ice = MLErayleigh(x_ice)

%% R1.d) Estimated distributions and histogram of data
scale_water = 0:1:max(x_water);
water_rayleigh = raylpdf(scale_water,sigma_water);
water_exp = exppdf(scale_water,lambda_water);
water_normal = normpdf(scale_water,mu_water,sqrt(var_water));

scale_ice = 0:1:max(x_ice);
ice_rayleigh = raylpdf(scale_ice,sigma_ice);
ice_exp = exppdf(scale_ice,lambda_ice);
ice_normal = normpdf(scale_ice,mu_ice,sqrt(var_ice));

% water
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',25);
xlabel("Intensity",'Interpreter','latex');
h_water = histogram(x_water,'Normalization','pdf');
plot(scale_water,water_rayleigh)
plot(scale_water,water_exp)
plot(scale_water,water_normal)
legend('Data','Rayleigh','Exponential','Normal')
saveas(gcf,'./figures/R1d_water.png');

water_h_values = h_water.Values';
water_h_edges = h_water.BinEdges';
water_h_ray_values = raylpdf(water_h_edges(1:end-1),sigma_water);
water_h_exp_values = exppdf(water_h_edges(1:end-1),lambda_water);
water_h_norm_values = normpdf(water_h_edges(1:end-1),mu_water,sqrt(var_water));
water_ray_rmse = sqrt(sumsqr(water_h_ray_values-water_h_values)/length(water_h_values));
water_exp_rmse = sqrt(sumsqr(water_h_exp_values-water_h_values)/length(water_h_values));
water_norm_rmse = sqrt(sumsqr(water_h_norm_values-water_h_values)/length(water_h_values));
fprintf("rmse water\nray:%g exp:%g norm: %g\n",water_ray_rmse,water_exp_rmse,water_norm_rmse);

% ice
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
set(gca,'FontSize',25);
xlabel("Intensity",'Interpreter','latex');
h_ice = histogram(x_ice,'Normalization','pdf');
plot(scale_ice,ice_rayleigh)
plot(scale_ice,ice_exp)
plot(scale_ice,ice_normal)
legend('Data','Rayleigh','Exponential','Normal')
saveas(gcf,'./figures/R1d_ice.png');

ice_h_values = h_ice.Values';
ice_h_edges = h_ice.BinEdges';
ice_h_ray_values = raylpdf(ice_h_edges(1:end-1),sigma_ice);
ice_h_exp_values = exppdf(ice_h_edges(1:end-1),lambda_ice);
ice_h_norm_values = normpdf(ice_h_edges(1:end-1),mu_ice,sqrt(var_ice));
ice_ray_rmse = sqrt(sumsqr(ice_h_ray_values-ice_h_values)/length(ice_h_values));
ice_exp_rmse = sqrt(sumsqr(ice_h_exp_values-ice_h_values)/length(ice_h_values));
ice_norm_rmse = sqrt(sumsqr(ice_h_norm_values-ice_h_values)/length(ice_h_values));
fprintf("rmse ice\nray:%g exp:%g norm: %g\n",ice_ray_rmse,ice_exp_rmse,ice_norm_rmse);

%% R2.a) Initial pixel-wise segmentation
close all

isWater = classification(I,sigma_water,sigma_ice);

% visualization
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$x$",'Interpreter','latex');
ylabel("$y$",'Interpreter','latex');
set(gca,'FontSize',25);
imagesc(I);
imcontour(isWater,1,'r');
saveas(gcf,'./figures/R2a.png');

%% R2.b) Segmentation using a 5x5 region around each pixel
% filtering
patch = 5;
kernel = ones(patch)/patch^2;
I_filt = conv2(I,kernel,'same');

isWater_filt = classification(I_filt,sigma_water,sigma_ice);

% visualization
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
grid on;
xlabel("$x$",'Interpreter','latex');
ylabel("$y$",'Interpreter','latex');
set(gca,'FontSize',25);
imagesc(I);
imcontour(isWater_filt,1,'r');
saveas(gcf,'./figures/R2b.png');

%% R2.c) Correct answers in the training regions
close all

% water
% in the rectangle of water gets the obtained classification
isWater_water = imcrop(isWater,rectWater);
% gets a rate of correctness of the estimator in a region know to be of
% water
water_rate = sum(sum(isWater_water))/...
    (size(isWater_water,1)*size(isWater_water,2));
fprintf("In the water training region, the accuracy was %.4g\n",water_rate)

% ice
isWater_ice = imcrop(isWater,rectIce);
ice_rate = 1 - sum(sum(isWater_ice))/...
    (size(isWater_ice,1)*size(isWater_ice,2));
fprintf("In the ice training region, the accuracy was %.4g\n",ice_rate)

% water after filtering
isWater_water_filt = imcrop(isWater_filt,rectWater);
water_rate_filt = sum(sum(isWater_water_filt))/...
    (size(isWater_water_filt,1)*size(isWater_water_filt,2));
fprintf("In the water training region after filtering, the accuracy was %.4g\n",water_rate_filt)

% ice after filtering
isWater_ice_filt = imcrop(isWater_filt,rectIce);
ice_rate_filt = 1 - sum(sum(isWater_ice_filt))/...
    (size(isWater_ice_filt,1)*size(isWater_ice_filt,2));
fprintf("In the ice training region after filtering, the accuracy was %.4g\n",ice_rate_filt)