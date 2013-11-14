clc
clear all
close all

%%% statistical simulation of generalized product capability %%%

% 1 simulate effect of converting from IT dist to PCSL & dim and back

mu = 10;
sigma = 1;
pd = makedist('Normal', 'mu', mu, 'sigma', sigma);

plotInitialDist = 0
if (plotInitialDist)
    x = linspace(6,14,100);
    xpdf = pdf(pd, x);

    figure()
    plot(x,xpdf)
end

% create N number from the normal dist
N = 1000;
rng('default'); % for reproducibility

ITvals = icdf(pd, rand(N,1))

%histfit(ITvals)

dimension = 100;
dimensions = ones(N,1)*dimension;
tolerenceWidths = Utilities.dimAndITGradeToTol(dimensions, ITvals);

histfit(tolerenceWidths)

ITvals_2 = Utilities.dimAndTolToITGrade(dimensions, tolerenceWidths);

[ITvals, dimensions, tolerenceWidths, ITvals-ITvals_2]

%r = normrnd(10,1,10,1)

%% 

clear all
close all
clc

c_pk = 5/3; % 1.6667
target = 100; % mm  (m)
ITG = 10;
C_a = 0.1;

% C_a = 1 - |mu - m| / d
% C_a * d = d - |mu - m|
% d - C_a * d = mu - m   // release abs
% mu = d - C_a * d + m

tolerenceWidth = Utilities.dimAndITGradeToTol(target, ITG);
d = tolerenceWidth/2; % half specification width

mu = d - C_a * d + target;

sigma = (d - abs(mu - target) ) /(3*c_pk);

sample_size = 3; 
sample_sets = 10;

pd = makedist('Normal', 'mu', mu, 'sigma', sigma);

%rng('default');

runs = 1000;
h = waitbar(0,'Running montecarlo simulation')

for i= 1:runs
    for j = 1:sample_sets     
        samples(:,j) = icdf(pd, rand(sample_size,1));

        sample_std(j) = std(samples(:,j));
        sample_mean(j) = mean(samples(:,j));
        sample_meanshift(j) = target - sample_mean(j);
        sample_PCSL(j) = Utilities.meanshiftAndStdAndCpkToPCSL(sample_meanshift(j),sample_std(j),c_pk);

        ITvals_3(j) = Utilities.dimAndTolToITGrade(target,sample_PCSL(j)*2);
    end

    %[sample_std' , sample_meanshift', sample_tol', ITvals_3']

    sets_mean(i) = mean(ITvals_3);
    
    waitbar(i/runs)
end
close(h)

monte_std = std(sets_mean);
monte_mean = mean(sets_mean)
pd = makedist('Normal', 'mu', monte_mean, 'sigma', monte_std);

lowerConfidenceLimit = icdf(pd, 0.025)
upperConfidenceLimit = icdf(pd, 0.975)
confidenreceIntevalWidth = upperConfidenceLimit-lowerConfidenceLimit

histfit(sets_mean)

