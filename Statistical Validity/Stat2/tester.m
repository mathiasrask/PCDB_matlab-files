clc
clear all
close all

%parpool('awersome',4)

sample_sets = [20];
sample_size = 12;
ITG_std = [0.1, 0.15,  0.2, 0.3, 0.5, 1, 2, 3];

var1 = sample_sets;
var2 = ITG_std;

count_max = length(var1) * length(var2);
count = 0;

h = waitbar(0,'Running montecarlo simulation');

CLW90 = zeros(length(var1),length(var2));

for i = 1:length(var1)
    for j = 1:length(var2)
        tic
        count = count + 1;
        CLW90(i,j) = Utilities.montecarloCLwidthSTD(1000,sample_size,sample_sets(i),0.9, ITG_std(j));
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

    plot(sample_sets,CLW90(:,:))
    %plot(sample_sets,CLW90(:))
hold off
%%
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', 'Times New Roman')
% set(0,'DefaultAxesFontSize', 9)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', 'Times New Roman')
% set(0,'DefaultTextFontSize', 9)
% 
% addpath('lib/export_fig');
% 
% f = fit(sample_sets',CLW90,'exp2')
% 
% f2 = figure ()
% plot(f,sample_sets,CLW90)
% xlabel('No. Measurement sets')
% ylabel('Confidence width')
% set(gcf, 'Color', 'w');
% set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])
% 
% export_fig('regresion_sample_sets.pdf')

%% dicide section

for i=1:length(ITG_std)
   CLW90_norm1(:,i) = CLW90(:,i)./ITG_std(i);
end
%CLW90_norm2 = CLW90_norm1(1,:)./CLW90_norm1(2,:);