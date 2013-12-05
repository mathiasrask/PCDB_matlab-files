clc
clear all
close all

%parpool('awersome',4)

sample_sets = [1:10,12:2:50];
sample_size = [12];

count_max = length(sample_sets) * length(sample_size);
count = 0;

h = waitbar(0,'Running montecarlo simulation');

CLW90 = zeros(length(sample_sets),length(sample_size));

for i = 1:length(sample_sets)
    for j = 1:length(sample_size)
        tic
        count = count + 1;
        CLW90(i,j) = Utilities.montecarloCLwidth(1000,sample_size(j),sample_sets(i),0.9);
        t(count) = toc;
        time = max(t) * (count_max-count);
        waitbar(count/count_max,h,sprintf(strcat('Running montecarlo simulation. Time left:' , datestr(time/24/3600, 'HH:MM:SS'))))
    end
end
close(h)

%figure ()
%surf( sample_size, sample_sets, CLW90)

figure()
hold on

    %plot(sample_sets,CLW90(:,:))
    %plot(sample_sets,CLW90(:))
hold off
%%
% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 9)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 9)

addpath('lib/export_fig');

f = fit(sample_sets',CLW90,'exp2')

f2 = figure ()
plot(f,sample_sets,CLW90)
xlabel('No. Measurement sets')
ylabel('Confidence width')
set(gcf, 'Color', 'w');
set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

export_fig('regresion_sample_sets.pdf')