clc
clear all
close all
load('TData')

mu0 = TData(:,4);
bias = TData(:,5);
mu = mu0+bias;
stddiv = TData(:,6);
numberofMeasurements = TData(:,7);

figure()
hold on
plot(bias, stddiv, 'b.');


%% Calculate tolerance based on CPK 1.33

% Cpk = min(UL-mu, mu- UL)/(3 sigma) 
% mu = mu0 + b
% Cpk (3 sigma) = min(UL - (mu0 + b), (mu0 + b) - LL)
% Cpk (3 sigma) = min(UL - (mu0 + |b|), (mu0 + |b|) - LL)
% Cpk (3 sigma) = UL - mu0 - |b|
% UL - mu0 = symTol
% Cpk (3 sigma) = symTol - |b|
% symTol = b + 3 * sigma * Cpk


Cpk = 1.33;
symTol =  abs(bias) + Cpk * 3 * stddiv;
TData(:,8) = symTol;
