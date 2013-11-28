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

figure()


end 
%r = normrnd(10,1,10,1)

%%  Calculate 
if (1)

c_pk = 5/3; % 1.6667
target = 100; % mm  (m)

sample_size = 100; 
sample_sets = 20;

biasCorrectionFactor_c4 = sqrt(2/(sample_size - 1)) * gamma(sample_size/2)/gamma((sample_size-1)/2)


ITG_mean = 10;
ITG_std = 1;
C_a = 0.6;

%ITGradesPD = makedist('Normal', 'mu', ITG_mean, 'sigma', ITG_std);


runs = 1000;
h = waitbar(0,'Running montecarlo simulation');

probability =0.025:0.025:0.975;

sets_mean = zeros( sample_sets, 1);
sets_std = zeros( sample_sets, 1);
ITG_samplesets = zeros(sample_sets, 1);
ITGatProperbility = zeros(sample_sets,length(probability));

for i= 1:runs                
    ITG = norminv(rand(sample_sets,1),ITG_mean,ITG_std);
    tolerenceWidth = Utilities.dimAndITGradeToTol(target, ITG);
    d = tolerenceWidth/2; % half specification width
    mu = (1-C_a) * d + target;
    sigma = (d - abs(mu - target) ) /(3*c_pk);

    for j = 1:sample_sets  
        
        samples = norminv(rand(sample_size,1),mu(j),sigma(j));
        sample_std = std(samples)/biasCorrectionFactor_c4;
        sample_mean = mean(samples);
        sample_meanshift = target - sample_mean;
        sample_PCSL = Utilities.meanshiftAndStdAndCpkToPCSL(sample_meanshift,sample_std,c_pk);

        ITG_samplesets(j) = Utilities.dimAndTolToITGrade(target,sample_PCSL*2);
    end

    sets_mean(i) = mean(ITG_samplesets);
    sets_std(i) = std(ITG_samplesets);

    for k = 1:length(probability)
        ITGatProperbility(i,k) = norminv(probability(k),sets_mean(i),sets_std(i));
    end
    
    if (rem(i,runs/20)==0)
        waitbar(i/runs)
    end
end
close(h)

for i = 1:length(probability)
    [low(i),mid(i),upp(i),wid(i)] = Utilities.ConfidenceLimit(ITGatProperbility(:,i));
end

% probability(20) = 0.5

mean = mean(sets_mean)

figure ()
hold on
plot(low,probability,'r')
plot(mid,probability,'b')
plot(upp,probability,'r')

x = linspace(6,14,100);

%xpdf = normpdf(x, ITG_mean, ITG_std);
xcdf = normcdf(x, ITG_mean, ITG_std);
[xc, wsi] = Utilities.wilson(xcdf*sample_sets,sample_sets,0.05);
%[xc, wsi] = binofit(xcdf*sample_sets,sample_sets,0.05);

plot(x,xc, 'b-.')
plot(x, wsi(:,1), 'r-.')
plot(x, wsi(:,2), 'r-.')


figure()
plot(probability, wid)


end


