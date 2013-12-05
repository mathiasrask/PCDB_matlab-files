clc
clear all
close all

addpath('../lib/export_fig');

if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 8
end

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
runs = 10000;
h = waitbar(0,'Running montecarlo simulation');
    
    
c_pk = 5/3; % 1.6667
target = 100; % mm  (m)

sample_size = 10; 
sample_sets = 20;

biasCorrectionFactor_c4 = sqrt(2/(sample_size - 1)) * gamma(sample_size/2)/gamma((sample_size-1)/2);
biasCorrectionFactor_c4sampleSets = sqrt(2/(sample_sets - 1)) * gamma(sample_sets/2)/gamma((sample_sets-1)/2);

%biasCorrectionFactor_c4 = 1;

ITG_mean = 10;
ITG_std = 1;
C_a = 1;

%ITGradesPD = makedist('Normal', 'mu', ITG_mean, 'sigma', ITG_std);




probability =0.025:0.025:0.975;

sets_mean = zeros( sample_sets, 1);
sets_std = zeros( sample_sets, 1);
ITGatProperbility = zeros(runs,length(probability));

parfor i= 1:runs       
    
    ITG_samplesets = zeros(sample_sets, 1);
    
    
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
    sets_std(i) = std(ITG_samplesets)/biasCorrectionFactor_c4sampleSets;
    
%     if (rem(i,runs/20)==0)
%         waitbar(i/runs)
%     end
end

for i = 1:runs
   
    for k = 1:length(probability)
        ITGatProperbility(i,k) = norminv(probability(k),sets_mean(i),sets_std(i));
    end
    
end

close(h)

for i = 1:length(probability)
    [low(i),mid(i),upp(i),wid(i)] = Utilities.ConfidenceLimit(ITGatProperbility(:,i));
end

% probability(20) = 0.5

ITG_mean_cal= mean(sets_mean)
ITG_std_cal = mean(sets_std)



%% PLOT
% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 9)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 9)


f1 = figure ()
hold on


x = linspace(7,13,100);

xcdf = normcdf(x, ITG_mean, ITG_std);
[xc, wsi] = Utilities.wilson(xcdf*sample_sets,sample_sets,0.05);
%[xc, wsi] = binofit(xcdf*sample_sets,sample_sets,0.05);

plot(x,xcdf, 'k-.')
plot(mid,probability,'b')
plot(low,probability,'r')
gp1 = plot(x, wsi(:,1), '--')
plot(upp,probability,'r')
gp2 = plot(x, wsi(:,2), '--')

set(gp1,'Color',[0 .5 0]);
set(gp2,'Color',[0 .5 0]);


legend('input', 'MC', 'MC conf.', 'Wil. conf.', 'location', 'SouthEast')
xlabel('Process capability driven tolerance (IT Grade)')
ylabel('Probability')
xlim([7 13])



set(gcf, 'Color', 'w');
set(f1, 'units', 'centimeters', 'pos', [0 0 8 8])


export_fig('confidenceIntervals.pdf')


figure()
plot(probability, wid)


end


