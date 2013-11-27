clc
clear all
close all

%parpool('awersome',4)

sample_sets = [2,3,4,5,6,7,8,9,10:5:100];
sample_size = 12;

count_max = length(sample_sets) * length(sample_size);
count = 0;

h = waitbar(0,'Running montecarlo simulation');

CLW90 = zeros(length(sample_sets),length(sample_size));

for i = 1:length(sample_sets)
    for j = 1:length(sample_size)
        tic
        count = count + 1;
        CLW90(i,j) = Utilities.montecarloCLwidth(100,sample_size(j),sample_sets(i),0.9);
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
    plot(sample_sets,CLW90(:))
hold off
%%
figure ()
f = fit(sample_sets',CLW90,'exp2')
plot(f,sample_sets,CLW90)