clear all; close all; clc;

ncdisp('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_1.nc');
% ncdisp('C:\Users\Ayden\NOAA Work\w_1.nc'); % Use this on Laptop

ncdisp('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\sst0n80.5e_dy.cdf');
% ncdisp('C:\Users\Ayden\NOAA Work\sst0n80.5e_dy.cdf'); % Laptop

w_1_depth = ncread('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_1.nc', 'depth');
% w_1_depth = ncread('C:\Users\Ayden\NOAA Work\w_1.nc', 'depth'); % Laptop

w_1_time = ncread('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_1.nc', 'time');
% w_1_time = ncread('C:\Users\Ayden\NOAA Work\w_1.nc', 'time'); % Laptop

w_1_time = datetime(2004,10,27,'Format','dd-MMM-yyyy HH:mm:ss') + days(w_1_time);
w_1 = ncread("C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_1.nc", '__xarray_dataarray_variable__');
% w_1 = ncread("C:\Users\Ayden\NOAA Work\w_1.nc", '__xarray_dataarray_variable__'); % Laptop

w_1 = w_1(:, 141:1812);
w_1_t = w_1_time(141:1812); % Indexed from 01-01-2009 to 30-07-2013 (same for tau_x below to match timescales)

w_1_m = movmean(w_1, 15, 2, 'omitnan');


% ncdisp('C:\Users\Ayden\NOAA Work\tau_x.nc')
tau_x = ncread("C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\tau_x.nc", 'taux');
% tau_x = ncread("C:\Users\Ayden\NOAA Work\tau_x.nc", 'taux'); % Laptop

tau_x_lon = ncread('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\tau_x.nc', 'longitude'); % index = 51 for 80.5E
% tau_x_lon = ncread('C:\Users\Ayden\NOAA Work\tau_x.nc', 'longitude'); % Laptop

tau_x_lat = ncread("C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\tau_x.nc", 'latitude'); % index = 31 for .5 N and 30 for .5 S
% tau_x_lat = ncread("C:\Users\Ayden\NOAA Work\tau_x.nc", 'latitude'); % Laptop

tau_x_time = ncread("C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\tau_x.nc", 'time');
% tau_x_time = ncread("C:\Users\Ayden\NOAA Work\tau_x.nc", 'time'); % Laptop

tau_x_time = datetime(1950,1,1,'Format','dd-MMM-yyyy HH:mm:ss') + days(tau_x_time);
tau_x_t  = tau_x_time(1:1672); % This is specific to our indexing from above
tau_x1 = tau_x(51, 30, 1:1672); % This is .5 S
tau_x2 = tau_x(51, 31, 1:1672); % This is .5 N
tau_x_1 = (tau_x1 + tau_x2)/2; % This is averaged tau_x between .5N and .5S

tau_x3 = tau_x(51, 29, 1:1672); % This is 1.5 S
tau_x4 = tau_x(51, 32, 1:1672); % This is 1.5 N
tau_x_2 = (tau_x3 + tau_x4)/2; % This is averaged tau_x between 1.5N and 1.5S

tau_x5 = tau_x(51, 28, 1:1645); % This is 2.5 S
tau_x6 = tau_x(51, 33, 1:1645); % This is 2.5 N
tau_x_3 = (tau_x5 + tau_x6)/2; % This is averaged tau_x between 2.5N and 2.5S

tau_x_m = movmean(tau_x_1, 15, 'omitnan');
tau_x_m = squeeze(tau_x_m); 
tau_x_m = tau_x_m';

tau_x_m_2 = movmean(tau_x_2, 15, 'omitnan');
tau_x_m_2 = squeeze(tau_x_m_2); 
tau_x_m_2 = tau_x_m_2';

tau_x_m_3 = movmean(tau_x_3, 15, 'omitnan');
tau_x_m_3 = squeeze(tau_x_m_3); 
tau_x_m_3 = tau_x_m_3';

m_t = mean(tau_x_m); % Going forward, maybe take mean of tau_x and w_1, not 15 day running means
m_w = mean(w_1_m, 2);

m_t_2 = mean(tau_x_m_2);

m_t_3 = mean(tau_x_m_3);


corr = zeros(41, 61);
w = zeros(41, 1672);
for i = 1:41
    w = w_1_m(i, :) - m_w(i);
    corr(i, :) = xcorr(w, (tau_x_m - m_t), 30, 'normalized');
end

% w_1_trial = w_1_m(3, :);
% corr = xcorr((tau_x_m - mean(tau_x_m)), (w_1_trial - mean(w_1_trial)), 30, 'normalized');
lags = (-30:30);
% plot(lags, squeeze(corr));

figure(1);
corr_plot = pcolor(lags, w_1_depth, corr);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with Zonal Wind Stress 15 RM (+ Lag, Wind Leads)');
% saveas(corr_plot,'Corr_W_Tau_15RM_Separ_1.png');


hold on
% sst_1 = ncread('C:\Users\Ayden\NOAA Work\sst0n80.5e_dy.cdf', 'T_25'); % Laptop
% sst_1 = ncread('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\sst0n80.5e_dy.cdf', 'T_25');

% sst_t = ncread('C:\Users\Ayden\NOAA Work\sst0n80.5e_dy.cdf', 'time'); % Laptop
% sst_t = ncread('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\sst0n80.5e_dy.cdf', 'time');

% sst_t = datetime(2009,1,1,'Format','dd-MMM-yyyy HH:mm:ss') + days(sst_t);
% sst_t = sst_t(1:1672);
% sst_1 = sst_1(:,:,1,1:1672);

% sst_1_m = movmean(sst_1, 15, 4, 'omitnan'); % 7/24: Lots of NaNs showing up around i = 700.. Not sure how to filter
% sst_1_m = squeeze(sst_1_m);
% sst_1_m = sst_1_m';

% m_sst = mean(sst_1_m);


% corr_2 = zeros(41, 61);
% w2 = zeros(41, 1672);
% for i = 1:41
%     w2 = w_1_m(i, :) - m_w(i);
%     corr_2(i, :) = xcorr(w2, (sst_1_m - m_sst), 30, 'normalized');
% end

% This is for 31 day running mean
w_1_m31 = movmean(w_1, 31, 2, 'omitnan');

tau_x_m31 = movmean(tau_x_1, 31, 'omitnan');
tau_x_m31 = squeeze(tau_x_m31); 
tau_x_m31 = tau_x_m31';

m_t31 = mean(tau_x_m31); % Going forward, maybe take mean of tau_x and w_1, not 15 day running means
m_w31 = mean(w_1_m31, 2);

corr31 = zeros(41, 61);
w31 = zeros(41, 1672);
for i = 1:41
    w31 = w_1_m31(i, :) - m_w31(i);
    corr31(i, :) = xcorr(w31, (tau_x_m31 - m_t31), 30, 'normalized');
end


figure(2);
corr31_plot = pcolor(lags, w_1_depth, corr31);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with Zonal Wind Stress 31 RM (+ Lag, Wind Leads)');
% saveas(corr31_plot,'Corr_W_Tau_31RM_Separ_1.png');

% ncdisp("w_2.nc");
w_2_time = ncread('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_2.nc', 'time');
% w_2_time = ncread('C:\Users\Ayden\NOAA Work\w_2.nc', 'time'); % Laptop

w_2_time = datetime(2004,10,27,'Format','dd-MMM-yyyy HH:mm:ss') + days(w_2_time); % This goes to 18-08-2014
w_2 = ncread("C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_2.nc", '__xarray_dataarray_variable__');
% w_2 = ncread("C:\Users\Ayden\NOAA Work\w_2.nc", '__xarray_dataarray_variable__'); % Laptop

w_2 = w_2(:, 139:1810);
w_2_t = w_2_time(139:1810); % Indexed from 01-01-2009 to 30-07-2013 (same for tau_x below to match timescales)

w_2_m = movmean(w_2, 15, 2, 'omitnan');

% ncdisp("w_3.nc");
w_3_time = ncread('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_3.nc', 'time');
% w_3_time = ncread('C:\Users\Ayden\NOAA Work\w_3.nc', 'time'); % Laptop

w_3_time = datetime(2004,10,27,'Format','dd-MMM-yyyy HH:mm:ss') + days(w_3_time); % This goes to 06-07-2013
w_3 = ncread("C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_3.nc", '__xarray_dataarray_variable__');
% w_3 = ncread("C:\Users\Ayden\NOAA Work\w_3.nc", '__xarray_dataarray_variable__'); % Laptop

w_3 = w_3(:, 139:1783);
w_3_t = w_3_time(139:1783); % Indexed from 01-01-2009 to 03-07-2013 (same for tau_x below to match timescales)

w_3_m = movmean(w_3, 15, 2, 'omitnan');


% This is for w_2:

m_w_2 = mean(w_2_m, 2);
corr_2 = zeros(41, 61);
w_1_5 = zeros(41, 1672);
for i = 1:41
    w_1_5 = w_2_m(i, :) - m_w_2(i);
    corr_2(i, :) = xcorr(w_1_5, (tau_x_m_2 - m_t_2), 30, 'normalized');
end

figure(3);
corr_2_plot = pcolor(lags, w_1_depth, corr_2);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_2_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with Zonal Wind Stress 15 RM (1.5N to 1.5S)');
% saveas(corr_2_plot,'Corr_W_Tau_15RM_Separ_2.png');

% This is for 31 day running mean w_2
w_2_m31 = movmean(w_2, 31, 2, 'omitnan');

tau_x_m31_2 = movmean(tau_x_2, 31, 'omitnan');
tau_x_m31_2 = squeeze(tau_x_m31_2); 
tau_x_m31_2 = tau_x_m31_2';

m_t31_2 = mean(tau_x_m31_2); % Going forward, maybe take mean of tau_x and w_1, not 15 day running means
m_w31_2 = mean(w_2_m31, 2);

corr31_2 = zeros(41, 61);
w31_2 = zeros(41, 1672);
for i = 1:41
    w31_2 = w_2_m31(i, :) - m_w31_2(i);
    corr31_2(i, :) = xcorr(w31_2, (tau_x_m31_2 - m_t31_2), 30, 'normalized');
end


figure(4);
corr31_2_plot = pcolor(lags, w_1_depth, corr31_2);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_2_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with Zonal Wind Stress 31 RM (1.5N to 1.5S)');
% saveas(corr31_2_plot,'Corr_W_Tau_31RM_Separ_2.png');

% This is for w_3:

m_w_3 = mean(w_3_m, 2);
corr_3 = zeros(41, 61);
w_2_5 = zeros(41, 1672);
for i = 1:41
    w_2_5 = w_3_m(i, :) - m_w_3(i);
    corr_3(i, :) = xcorr(w_2_5, (tau_x_m_3 - m_t_3), 30, 'normalized');
end

figure(5);
corr_3_plot = pcolor(lags, w_1_depth, corr_3);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_3_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with Zonal Wind Stress 15 RM (2.5N to 2.5S)');
% saveas(corr_3_plot,'Corr_W_Tau_15RM_Separ_3.png');

% This is for 31 day running mean w_3
w_3_m31 = movmean(w_3, 31, 2, 'omitnan');

tau_x_m31_3 = movmean(tau_x_3, 31, 'omitnan');
tau_x_m31_3 = squeeze(tau_x_m31_3); 
tau_x_m31_3 = tau_x_m31_3';

m_t31_3 = mean(tau_x_m31_3); % Going forward, maybe take mean of tau_x and w_1, not 15 day running means
m_w31_3 = mean(w_3_m31, 2);

corr31_3 = zeros(41, 61);
w31_3 = zeros(41, 1672);
for i = 1:41
    w31_3 = w_3_m31(i, :) - m_w31_3(i);
    corr31_3(i, :) = xcorr(w31_3, (tau_x_m31_3 - m_t31_3), 30, 'normalized');
end


figure(6);
corr31_3_plot = pcolor(lags, w_1_depth, corr31_3);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_3_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with Zonal Wind Stress 31 RM (2.5N to 2.5S)');
% saveas(corr31_3_plot,'Corr_W_Tau_31RM_Separ_3.png');

t = tiledlayout(2, 3);
title(t, 'Cross-Correlation of Vertical Velocities and Zonal Wind Stress');
nexttile;
corr_plot = pcolor(lags, w_1_depth, corr);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('15 RM (0.5N to 0.5S)');

nexttile;
corr_2_plot = pcolor(lags, w_1_depth, corr_2);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_2_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('15 RM (1.5N to 1.5S)');

nexttile;
corr_3_plot = pcolor(lags, w_1_depth, corr_3);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_3_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('15 RM (2.5N to 2.5S)');

nexttile;
corr31_plot = pcolor(lags, w_1_depth, corr31);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('31 RM (0.5N to 0.5S)');

nexttile;
corr31_2_plot = pcolor(lags, w_1_depth, corr31_2);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_2_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('31 RM (1.5N to 1.5S)');

nexttile;
corr31_3_plot = pcolor(lags, w_1_depth, corr31_3);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_3_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('31 RM (2.5N to 2.5S)');


% ncdisp("sst.day.mean.2009.nc");
% sst_t_1 = ncread("sst.day.mean.2009.nc", 'time');
% sst_t_2 = ncread("sst.day.mean.2010.nc", 'time');
% sst_t_3 = ncread("sst.day.mean.2011.nc", 'time');
% sst_t_4 = ncread("sst.day.mean.2012.nc", 'time');
% sst_t_5 = ncread("sst.day.mean.2013.nc", 'time');
% sst_t = [sst_t_1; sst_t_2; sst_t_3; sst_t_4; sst_t_5];
% sst_t = datetime(1800,1,1,'Format','dd-MMM-yyyy HH:mm:ss') + days(sst_t);
sst_lon = ncread("sst.day.mean.2009.nc", 'lon'); % average between 323 and 322 for 80.5E

sst_lat = ncread("sst.day.mean.2009.nc", 'lat'); % for 0.5N to 0.5S: avg(avg(362,363) and avg(358,359)) interval: 361-360
% For 1.5N to 1.5S: avg(avg(366,367) and avg(354,355)) interval: 365-356
% For 2.5N to 2.5S: avg(avg(370,371) and avg(350,351)) interval: 369-352

ncdisp("sst_new (1).nc")
sst = ncread("sst_new (1).nc", 'sst');
sst_latitude = ncread("sst_new (1).nc", 'lat');
sst_time = ncread("sst_new (1).nc", 'time');
sst_time = datetime(1800,1,1,'Format','dd-MMM-yyyy HH:mm:ss') + days(sst_time);
sst = (sst(1,:,:) + sst(2,:,:)) / 2;
sst_1_outer = (sst(:, 13, :) + sst(:, 14, :)) / 2;
sst_1_inner = (sst(:, 9, :) + sst(:, 10, :)) / 2;
sst_1 = (sst_1_outer + sst_1_inner + sst(:, 11, :) + sst(:, 12, :)) / 4;
sst_1 = sst_1(:,:,1:1672);

sst_1_m = movmean(sst_1, 15, 'omitnan');
sst_1_m = squeeze(sst_1_m);
sst_1_m = sst_1_m';

m_sst_1 = mean(sst_1_m);

corr_sst_1 = zeros(41, 61);
w_sst_1 = zeros(41, 1672);
for i = 1:41
    w_sst_1 = w_1_m(i, :) - m_w(i);
    corr_sst_1(i, :) = xcorr(w_sst_1, (sst_1_m - m_sst_1), 30, 'normalized');
end

figure(7);
corr_sst_1_plot = pcolor(lags, w_1_depth, corr_sst_1);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_sst_1_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with SST 15 RM (0.5N to 0.5S)');
saveas(corr_sst_1_plot,'Corr_W_SST_15RM_Separ_1.png');

sst31_1_m = movmean(sst_1, 31, 'omitnan');
sst31_1_m = squeeze(sst31_1_m);
sst31_1_m = sst31_1_m';

m31_sst_1 = mean(sst31_1_m);

corr31_sst_1 = zeros(41, 61);
w31_sst_1 = zeros(41, 1672);
for i = 1:41
    w31_sst_1 = w_1_m31(i, :) - m_w31(i);
    corr31_sst_1(i, :) = xcorr(w31_sst_1, (sst31_1_m - m31_sst_1), 30, 'normalized');
end


figure(8);
corr31_sst_1_plot = pcolor(lags, w_1_depth, corr31_sst_1);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_sst_1_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with SST 31 RM (0.5N to 0.5S)');
saveas(corr31_sst_1_plot,'Corr_W_SST_31RM_Separ_1.png');

% For SST and W2:

sst_2_outer = (sst(:, 17, :) + sst(:, 18, :)) / 2;
sst_2_inner = (sst(:, 5, :) + sst(:, 6, :)) / 2;

sst_2 = (sst_2_outer + sst_2_inner + sst(:, 7, :) + sst(:, 8, :) + sst(:, 9, :) + sst(:, 10, :) + sst(:, 11, :) + sst(:, 12, :) + sst(:, 13, :) + sst(:, 14, :) + sst(:, 15, :) + sst(:, 16, :)) / 12;

sst_2 = sst_2(:,:,1:1672);

sst_2_m = movmean(sst_2, 15, 'omitnan');
sst_2_m = squeeze(sst_2_m);
sst_2_m = sst_2_m';

m_sst_2 = mean(sst_2_m);

corr_sst_2 = zeros(41, 61);
w_sst_2 = zeros(41, 1672);
for i = 1:41
    w_sst_2 = w_2_m(i, :) - m_w_2(i);
    corr_sst_2(i, :) = xcorr(w_sst_2, (sst_2_m - m_sst_2), 30, 'normalized');
end

figure(9);
corr_sst_2_plot = pcolor(lags, w_1_depth, corr_sst_2);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_sst_2_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with SST 15 RM (1.5N to 1.5S)');
saveas(corr_sst_2_plot,'Corr_W_SST_15RM_Separ_2.png');

sst31_2_m = movmean(sst_2, 31, 'omitnan');
sst31_2_m = squeeze(sst31_2_m);
sst31_2_m = sst31_2_m';

m31_sst_2 = mean(sst31_2_m);

corr31_sst_2 = zeros(41, 61);
w31_sst_2 = zeros(41, 1672);
for i = 1:41
    w31_sst_2 = w_2_m31(i, :) - m_w31_2(i);
    corr31_sst_2(i, :) = xcorr(w31_sst_2, (sst31_2_m - m31_sst_2), 30, 'normalized');
end


figure(8);
corr31_sst_2_plot = pcolor(lags, w_1_depth, corr31_sst_2);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_sst_2_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with SST 31 RM (1.5N to 1.5S)');
saveas(corr31_sst_2_plot,'Corr_W_SST_31RM_Separ_2.png');

% For SST and W3:

sst_3_inner = (sst(:, 1, :) + sst(:, 2, :)) / 2;
sst_3_outer = (sst(:, 21, :) + sst(:, 22, :)) / 2;
sst_3 = (sst_3_outer + sst_3_inner + sst(:, 3, :) + sst(:, 4, :) + sst(:, 5, :) + sst(:, 6, :) + sst(:, 7, :) + sst(:, 8, :) + sst(:, 9, :) + sst(:, 10, :) ...
    + sst(:, 11, :) + sst(:, 12, :) + sst(:, 13, :) + sst(:, 14, :) + sst(:, 15, :) + sst(:, 16, :) + sst(:, 17, :) + sst(:, 18, :) + sst(:, 19, :) + sst(:, 20, :)) / 20;

sst_3 = sst_3(:, :, 1:1645);

sst_3_m = movmean(sst_3, 15, 'omitnan');
sst_3_m = squeeze(sst_3_m);
sst_3_m = sst_3_m';

m_sst_3 = mean(sst_3_m);

corr_sst_3 = zeros(41, 61);
w_sst_3 = zeros(41, 1672);
for i = 1:41
    w_sst_3 = w_3_m(i, :) - m_w_3(i);
    corr_sst_3(i, :) = xcorr(w_sst_3, (sst_3_m - m_sst_3), 30, 'normalized');
end

figure(9);
corr_sst_3_plot = pcolor(lags, w_1_depth, corr_sst_3);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_sst_3_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with SST 15 RM (2.5N to 2.5S)');
saveas(corr_sst_3_plot,'Corr_W_SST_15RM_Separ_3.png');

sst31_3_m = movmean(sst_3, 31, 'omitnan');
sst31_3_m = squeeze(sst31_3_m);
sst31_3_m = sst31_3_m';

m31_sst_3 = mean(sst31_3_m);

corr31_sst_3 = zeros(41, 61);
w31_sst_3 = zeros(41, 1672);
for i = 1:41
    w31_sst_3 = w_3_m31(i, :) - m_w31_3(i);
    corr31_sst_3(i, :) = xcorr(w31_sst_3, (sst31_3_m - m31_sst_3), 30, 'normalized');
end


figure(10);
corr31_sst_3_plot = pcolor(lags, w_1_depth, corr31_sst_3);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_sst_3_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with SST 31 RM (2.5N to 2.5S)');
saveas(corr31_sst_3_plot,'Corr_W_SST_31RM_Separ_3.png');

t_2 = tiledlayout(2, 3);
title(t_2, 'Cross-Correlation of Vertical Velocities and SST');
nexttile;
corr_sst_1_plot = pcolor(lags, w_1_depth, corr_sst_1);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_sst_1_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('15 RM (0.5N to 0.5S)');

nexttile;
corr_sst_2_plot = pcolor(lags, w_1_depth, corr_sst_2);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_sst_2_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('15 RM (1.5N to 1.5S)');

nexttile;
corr_sst_3_plot = pcolor(lags, w_1_depth, corr_sst_3);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_sst_3_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('15 RM (2.5N to 2.5S)');

nexttile;
corr31_sst_1_plot = pcolor(lags, w_1_depth, corr31_sst_1);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_sst_1_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('31 RM (0.5N to 0.5S)');

nexttile;
corr31_sst_2_plot = pcolor(lags, w_1_depth, corr31_sst_2);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_sst_2_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('31 RM (1.5N to 1.5S)');

nexttile;
corr31_sst_3_plot = pcolor(lags, w_1_depth, corr31_sst_3);
set(gca, 'Ydir', 'reverse');
colorbar;
corr31_sst_3_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('31 RM (2.5N to 2.5S)');