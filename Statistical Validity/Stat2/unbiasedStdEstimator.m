clc
close all
clear

% test of unbiased estimator of standard deviation

mu = 0;
sigma = 1;
sampleSizeList = [3, 4, 5, 10];
sampleSizeListLength = length(sampleSizeList);
N = 10000;

std_samplesize = zeros(sampleSizeListLength, 1);
mean_samplesize = zeros(sampleSizeListLength, 1);

tic
for j = 1:sampleSizeListLength
    sampleSize = sampleSizeList(j);
    
    biasCorrectionFactor_c4 = sqrt(2/(sampleSize - 1)) * gamma(sampleSize/2)/gamma((sampleSize-1)/2);
    
    std_star_size = zeros(N,1);
    mean_star_size = zeros(N,1);
    
    parfor n=1:N
        samples = zeros(sampleSize,1);
        for i=1:sampleSize
            samples(i) = norminv(rand(),mu, sigma);
        end
        std_star_size(n) = std(samples)/biasCorrectionFactor_c4;
        mean_star_size(n) = mean(samples);
    end
    std_samplesize(j) = mean(std_star_size);
    mean_samplesize(j) = mean(mean_star_size);
    
end
toc


std_samplesize
mean_samplesize


