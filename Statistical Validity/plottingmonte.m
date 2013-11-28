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


load 'Stat2/CLW90.mat'

[x,y] = size(CLW90)

f2=figure()
hold on
box off
plot(sample_sets,CLW90)
legend('3', '7','11','15','19','23')
xlabel('Measurement sets')
ylabel('Cofidence interval limit width')
set(gcf, 'Color', 'w');
set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])
hold off
export_fig('CLW90_lines.pdf')

f3 = figure()
surf(sample_size,sample_sets,CLW90)
xlabel('Measurement sets')
ylabel('Sample size')
zlabel('Cofidence interval limit width')
set(gcf, 'Color', 'w');
set(f3, 'units', 'centimeters', 'pos', [0 0 8 8])

export_fig('CLW90_surf.pdf')