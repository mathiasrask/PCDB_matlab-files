clc
close all
clear all

%Generate Data

N=10;
std = 1;
mu = 200;

R = randn(N,1);

% Measure
mu_tmp = mean(R);
std_tmp = var(R)^(1/2);

% Normalise and denormalise
R = (R - mu_tmp) / std_tmp;
R = (R * std) + mu;


% Plot data

fit = fitdist(R,'Normal');

x_values = 190:.1:210;
y_values = pdf(fit,x_values);

x_hist = 190:1:210;
y_hist = hist(R,x_hist)/(N*2);

hold on
plot(x_values,y_values);
bar(x_hist,y_hist);