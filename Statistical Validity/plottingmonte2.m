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

export_fig('CLW90_lines.pdf')
