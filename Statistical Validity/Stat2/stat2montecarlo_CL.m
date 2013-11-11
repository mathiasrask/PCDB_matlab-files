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

sample_size = 12; 
sample_sets = 20;

ITG_mean = 10;
ITG_std = 1;

pd = makedist('Normal', 'mu', ITG_mean, 'sigma', ITG_std);

ITG = icdf(pd, rand(sample_sets,1));;
C_a = 0.2;

% C_a = 1 - |mu - m| / d
% C_a * d = d - |mu - m|
% d - C_a * d = mu - m   // release abs
% mu = d - C_a * d + m

tolerenceWidth = Utilities.dimAndITGradeToTol(target, ITG);
d = tolerenceWidth/2; % half specification width

mu = (1-C_a) * d + target;

sigma = (d - abs(mu - target) ) /(3*c_pk);



pd = makedist('Normal', 'mu', mu, 'sigma', sigma);

%rng('default');

runs = 100;
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
    sets_std(i) = std(ITvals_3);
    
    pd = makedist('Normal', 'mu', sets_mean, 'sigma', sets_std);
    
    
    CL10(i) = icdf(pd, 0.1);
    CL30(i) = icdf(pd, 0.3);
    CL50(i) = icdf(pd, 0.5);
    CL70(i) = icdf(pd, 0.7);
    CL90(i) = icdf(pd, 0.9);
   
    waitbar(i/runs)
end
close(h)


monte_std = std(CL90);
monte_mean = mean(CL90)
pd = makedist('Normal', 'mu', monte_mean, 'sigma', monte_std);

lowerConfidenceLimit = icdf(pd, 0.025)
upperConfidenceLimit = icdf(pd, 0.975)
confidenreceIntevalWidth = upperConfidenceLimit-lowerConfidenceLimit

histfit(CL90)

