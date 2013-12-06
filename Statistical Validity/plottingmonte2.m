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

% IT grade
for run=1:0
    load 'Stat2/Size3_sample10.mat'

    f2=figure()
    hold on
    plot(x,xc, 'K-.')
    plot(mid,probability,'b')
    plot(low,probability,'r')
    plot(upp,probability,'r')

    legend('Actual', 'Estimated','Confidence','location','northwest')
    xlabel('IT grade')
    ylabel('Probability')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])
end


% Standard deviation - confidence interval
for i=1:1
    load('Stat2/CLW90_std.mat')
    
    CLW90_std_sym = CLW90_norm1./2;
    
%     x = linspace(min(ITG_std), max(ITG_std), 100);
%     
%     p = polyfit(ITG_std,CLW90_std_sym,2)
%     
%     y = p(1)*x.^2 + p(2)*x + p(3); 

    f= fit(ITG_std',CLW90_std_sym','exp2')
    
    f2=figure()
    hold on
    plot(f,ITG_std,CLW90_std_sym)

    %legend('Actual', 'Estimated','Confidence','location','northwest')
    xlabel('Standard deviation of IT Grade input distribution')
    ylabel('Normalized symetric confidence interval')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])
    hLeg = legend('example')
    set(hLeg,'visible','off')
    
    export_fig('CLW90_std.pdf')
end
