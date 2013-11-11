clc
clear all
close all

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 9)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 9)


addpath('lib/export_fig');

v = 3:30;

sigma = 1;
S = sigma;
S2 = S^2;

alpha = 0.05;

confL = sqrt((v-1).*S2./chi2inv(1-alpha/2, v-1));
confH = sqrt((v-1).*S2./chi2inv(alpha/2, v-1));
sigmaseries = v*0+sigma;

f1 = figure();
hold on
plot(v,sigmaseries, 'r')
plot(v,confL, 'b-.')
plot(v,confH, 'b-.')
xlabel('n - number of samples')
ylabel('Standard deviation (Std)')
%title('Sample standard deviation')
legend('sample std.', '95 % confidence interval')
set(gcf, 'Color', 'w');
set(f1, 'units', 'centimeters', 'pos', [0 0 8 8])

export_fig('stats_std_confidence.pdf')



bias = 0;
biasseries = v*bias;
pd = makedist('Normal');
z_025 = icdf(pd,.025);
z_975 = icdf(pd,.974);
confL = z_025*sqrt(1./v);
confH = z_975*sqrt(1./v);

f2 = figure();
hold on
plot(v,biasseries, 'r')
plot(v,confL, 'b-.')
plot(v,confH, 'b-.')
xlabel('n - number of samples')
ylabel('Bias (std)')
legend('sample bias', '95 % confidence interval')
set(gcf, 'Color', 'w');
set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

export_fig('stats_bias_confidence.pdf')

%bConfL = 




