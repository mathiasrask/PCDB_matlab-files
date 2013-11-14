clc
clear all
close all

%%% statistical simulation of generalized product capability %%%

% 1 simulate effect of converting from IT dist to PCSL & dim and back
if (0)
mu = 10;
sigma = 1;
pd = makedist('Normal', 'mu', mu, 'sigma', sigma);

plotInitialDist = 1
if (plotInitialDist)
    x = linspace(6,14,100);
    xpdf = pdf(pd, x);

    figure()
    plot(x,xpdf)
end

% create N number from the normal dist
N = 1000;
%rng('default'); % for reproducibility

ITvals = icdf(pd, rand(N,1))

%histfit(ITvals)

dimension = 100;
dimensions = ones(N,1)*dimension;
tolerenceWidths = Utilities.dimAndITGradeToTol(dimensions, ITvals);

figure()
histfit(tolerenceWidths)

ITvals_2 = Utilities.dimAndTolToITGrade(dimensions, tolerenceWidths);

[ITvals, dimensions, tolerenceWidths, ITvals-ITvals_2]
end 
%r = normrnd(10,1,10,1)

%% 
if (1)

c_pk = 5/3; % 1.6667
target = 100; % mm  (m)

sample_size = 10; 
sample_sets = 20;

ITG_mean = 10;
ITG_std = 1;
C_a = 0.2;

ITGradesPD = makedist('Normal', 'mu', ITG_mean, 'sigma', ITG_std);


runs = 40;
h = waitbar(0,'Running montecarlo simulation')

probability =0.025:0.025:0.975;

sets_mean = zeros( sample_sets, 1);
sets_std = zeros( sample_sets, 1);
ITG_samplesets = zeros(sample_sets, 1);
ITGatProperbility = zeros(sample_sets,length(probability));

for i= 1:runs

    % Calculate the ITgrades for each sample set
    ITG = icdf(ITGradesPD, rand(sample_sets,1));
    
    % C_a = 1 - |mu - m| / d
    % C_a * d = d - |mu - m|
    % d - C_a * d = mu - m   // release abs
    % mu = d - C_a * d + m

    tolerenceWidth = Utilities.dimAndITGradeToTol(target, ITG);
    d = tolerenceWidth/2; % half specification width
    mu = (1-C_a) * d + target;
    sigma = (d - abs(mu - target) ) /(3*c_pk);
    
    
    for j = 1:sample_sets  
        
        sampleSetPD = makedist('Normal', 'mu', mu(j), 'sigma', sigma(j));
        
        samples = icdf(sampleSetPD, rand(sample_size,1));

        sample_std = std(samples);
        sample_mean = mean(samples);
        sample_meanshift = target - sample_mean;
        sample_PCSL = Utilities.meanshiftAndStdAndCpkToPCSL(sample_meanshift,sample_std,c_pk);

        ITG_samplesets(j) = Utilities.dimAndTolToITGrade(target,sample_PCSL*2);
    end

    %[sample_std' , sample_meanshift', sample_tol', ITvals_3']

    sets_mean(i) = mean(ITG_samplesets);
    sets_std(i) = std(ITG_samplesets);
    
    setsdist = makedist('Normal', 'mu', sets_mean(i), 'sigma', sets_std(i));
    
    
    for k = 1:length(probability)
        ITGatProperbility(i,k) = icdf(setsdist, probability(k));
    end
   
    waitbar(i/runs)
end
close(h)

for i = 1:length(probability)
    [low(i),mid(i),upp(i),wid(i)] = Utilities.ConfidenceLimit(ITGatProperbility(:,i));
end

figure ()
hold on
plot(low,probability,'r')
plot(mid,probability,'b')
plot(upp,probability,'r')

figure()
plot(probability, wid)


end



%%%


