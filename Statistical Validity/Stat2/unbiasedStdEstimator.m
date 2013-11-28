clc
close all
clear

% test of unbiased estimator of standard deviation

mu = 0
sigma = 1



tic
for n=1:100000
    norminv(rand(),mu, sigma);
end
toc