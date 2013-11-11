
clc
clear all
close all



% Data in
file = fopen('data/deep_groove_ball_bearing.csv');
D = textscan(file, '%f %f %s %s %d');

Dim = D{1};
Var = D{2}/1000;

%bearing = importdata('data/deep_groove_ball_bearing.csv');

%Dim = bearing(:,1);
%Var = abs(bearing(:,2)/1000);


%% Calculation IT grade
[ITG] = ITGrade(Dim,Var)';

%% sorting Base on size

F150 = [NaN, NaN];
F500 = [NaN, NaN];
F3000 = [NaN, NaN];

j =0;
k =0;
l =0;

for i = 1:length(Dim)
    if Dim(i) < 150
        j = j+1;
        F150(j,1) = ITG(i); 
    elseif Dim(i) < 500
        k= k+1;
        F500(k,1) = ITG(i);
    else 
        l =l+1;
        F3000(l,1) = ITG(i);
    end
end

F150(:,2) = linspace(0,1,length(F150));
F500(:,2) = linspace(0,1,length(F500));
F3000(:,2) = linspace(0,1,length(F3000));

F150 = sort(F150,1);
F500 = sort(F500,1);
F3000 = sort(F3000,1);

%% Sorting based on type

type = D{4};

j=0;
k=0;
Fnorm = [NaN, NaN];
Fexp = [NaN, NaN];
for i=1:length(Dim)
    if strcmp(type(i), 'Norm')
        j = j+1;
        Fnorm(j,1) = ITG(i);
    else
        k=k+1;
        Fexp(k,1) = ITG(i);
    end
end

Fnorm(:,2) = linspace(0,1,length(Fnorm));
Fexp(:,2) = linspace(0,1,length(Fexp));

Fnorm = sort(Fnorm,1);
Fexp = sort(Fexp,1);


%% plotting

figure(1)
hold on
plot(F150(:,1),F150(:,2)*100,'r.-')
plot(F500(:,1),F500(:,2)*100,'k.-')
plot(F3000(:,1),F3000(:,2)*100,'b.-')
legend('0-150 mm', '150-500 mm', '500-3000 mm', 'location','southeast')
title('Accumulated IT grade distribution - SKF Deep grove ball bearings')
xlabel('IT grade')
ylabel('Procent')
hold off

% figureHandle = gcf;
%# make all text in the figure to size 14 and bold
% set(findall(figureHandle,'type','text'),'fontSize',18,'fontWeight','bold')
% set(gca,'FontSize',18,'fontWeight','bold')

figure(3)
hold on
plot(Fexp(:,1),Fexp(:,2)*100,'r.-')
plot(Fnorm(:,1),Fnorm(:,2)*100,'k.-')
legend('exp', 'norm', 'location','southeast')
title('Accumulated IT grade distribution - SKF Deep grove ball bearings')
xlabel('IT grade')
ylabel('Procent')
hold off


% 
% figure (2)
% hold on
% xout = linspace(4,15,23);
% n= hist(ITG,xout);
% bar(xout,n)
% title('Total distribution of samples')
% hold off
% 
