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

goodvector = [1,2,4,5]

if 1
f2=figure()
hold on
box off
plot(sample_sets,CLW90(:,goodvector)./2)
legend('3', '5','12','20')
xlabel('Measurement sets')
ylabel('Symetric cofidence interval (IT Grade)')
set(gcf, 'Color', 'w');
set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])
hold off
export_fig('CLW90_lines.pdf')
end

if  0
f3 = figure()
surf(sample_size(goodvector),sample_sets,CLW90(:,goodvector))
xlabel('Measurement sets')
ylabel('Sample size')
zlabel('Cofidence interval limit width')
set(gcf, 'Color', 'w');
set(f3, 'units', 'centimeters', 'pos', [0 0 8 8])
export_fig('CLW90_surf.pdf')
end

if 0
    fi1 = fit(sample_sets',CLW90(:,1),'exp2')
    fi2 = fit(sample_sets',CLW90(:,2),'exp2')
    fi3 = fit(sample_sets',CLW90(:,4),'exp2')
    fi4 = fit(sample_sets',CLW90(:,5),'exp2')
    f4 = figure ()
    hold on
    plot(fi1,sample_sets',CLW90(:,1))
    plot(fi2,sample_sets',CLW90(:,2))
    plot(fi3,sample_sets',CLW90(:,4))
    plot(fi4,sample_sets',CLW90(:,:))

end