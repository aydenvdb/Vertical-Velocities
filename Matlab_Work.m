clear all; close all; clc;

ncdisp('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_1.nc');

ncdisp('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\sst0n80.5e_dy.cdf');

w_1_depth = ncread('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_1.nc', 'depth');
w_1_time = ncread('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_1.nc', 'time');
w_1_time = datetime(2004,10,27,'Format','dd-MMM-yyyy HH:mm:ss') + days(w_1_time);
w_1 = ncread("C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\w_1.nc", '__xarray_dataarray_variable__');
w_1 = w_1(:, 141:1812);
w_1_t = w_1_time(141:1812); % Indexed from 01-01-2009 to 30-07-2013 (same for tau_x below to match timescales)

w_1_m = movmean(w_1, 15, 2, 'omitnan');


ncdisp('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\tau_x.nc')
tau_x = ncread("C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\tau_x.nc", 'taux');
tau_x_lon = ncread('C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\tau_x.nc', 'longitude'); % index = 51 for 80.5E
tau_x_lat = ncread("C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\tau_x.nc", 'latitude'); % index = 31 for .5 N and 30 for .5 S
tau_x_time = ncread("C:\Users\Ayden van den Berg\Documents\MATLAB\NOAA Work\Data\tau_x.nc", 'time');
tau_x_time = datetime(1950,1,1,'Format','dd-MMM-yyyy HH:mm:ss') + days(tau_x_time);
tau_x_t  = tau_x_time(1:1672); % This is specific to our indexing from above
tau_x1 = tau_x(51, 30, 1:1672); % This is .5 S
tau_x2 = tau_x(51, 31, 1:1672); % This is .5 N
tau_x = (tau_x1 + tau_x2)/2; % This is averaged tau_x between .5N and .5S

tau_x_m = movmean(tau_x, 15, 'omitnan');
tau_x_m = squeeze(tau_x_m);
tau_x_m = tau_x_m';

m_t = mean(tau_x_m);
m_w = mean(w_1_m, 2);



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
% hold on

corr_plot = pcolor(lags, w_1_depth, corr);
set(gca, 'Ydir', 'reverse');
colorbar;
corr_plot.EdgeAlpha = .1;
xline(0, '--', 'r');
xlabel('Lag (days)');
ylabel('Depth (m)');
title('Cross-Correlation of W with Zonal Wind Stress 15 RM (+ Lag, Wind Leads)');
colorbar.Ticks = [-0.4, -0.3, -0.2, -0.1, 0 0.4];

