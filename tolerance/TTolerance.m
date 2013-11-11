
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


% Data in
load('data/TData');

Dim = TData(:,4);
Bias = TData(:,5);
Std = TData(:,6);
CPK = 1.66;
SymTol = (abs(Bias) + abs(Std).*CPK*3);
numberOfSamples = TData(:,7);
Machine = TData(:,9);

%% Calculation IT grade
ITG = ITGrade(Dim,SymTol)';
SYMTOL = sort(ITGradeINV(ITG,10));

%% sorting Base on size


AllData = [NaN, NaN];
j =0;

for i = 1:length(Dim)
    if Dim(i) < 150
        j = j+1;
        AllData(j,1) = ITG(i); 
    end
end

AllData(:,2) = linspace(0,1,length(AllData));
AllData = sort(AllData,1);


%% Sorting based on number of measurements

j=0;
k=0;
l=0;
F1031 = [NaN, NaN];
F1032 = [NaN, NaN];
F1145 = [NaN, NaN];


for i=1:length(Dim)
    if Machine(i) == 1031
        j = j+1;
        F1031(j,1) = ITG(i);
    elseif Machine(i) == 1032
        k=k+1;
        F1032(k,1) = ITG(i);
    else 
        l=l+1;
        F1145(l,1) = ITG(i);
    end
end

F1031(:,2) = linspace(0,1,length(F1031));
F1032(:,2) = linspace(0,1,length(F1032));
F1145(:,2) = linspace(0,1,length(F1145));

F1031 = sort(F1031,1);
F1032 = sort(F1032,1);
F1145 = sort(F1145,1);

%% Fitting Normal distribution


[x_lst , y_cdf] = plotCdf(AllData(:,1));

[x_1031 , y_1031] = plotCdf(F1031(:,1));
[x_1032 , y_1032] = plotCdf(F1032(:,1));
[x_1145 , y_1145] = plotCdf(F1145(:,1));


%% plotting

f1 = figure(1)
hold on
plot(AllData(:,1),AllData(:,2),'r.')
plot(x_lst,y_cdf)

% legend('500-3000 mm', 'location','southeast')
xlabel('IT grade')
%title('Accumulated IT grade distribution - Thornton actual data')
ylabel('Probability')
hold off
set(gcf, 'Color', 'w');

set(f1, 'units', 'centimeters', 'pos', [0 0 8 8])
export_fig('Acum_freqIT.pdf')

f2=figure(2)
hold on
plot(SYMTOL,AllData(:,2)*100,'r.-')

% legend('500-3000 mm', 'location','southeast')
% title('Accumulated IT grade distribution CPK 1.33')
xlabel('Sym Tolerance (mm) based on nominel 10 mm')
ylabel('Procent')
hold off
set(gcf, 'Color', 'w');

% figureHandle = gcf;
%# make all text in the figure to size 14 and bold
% set(findall(figureHandle,'type','text'),'fontSize',18,'fontWeight','bold')
% set(gca,'FontSize',18,'fontWeight','bold')

set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])
export_fig('Acum_freq.pdf')

f3 = figure(3)
hold on
plot(F1031(:,1),F1031(:,2),'r.')
plot(F1032(:,1),F1032(:,2),'k.')
%plot(F1145(:,1),F1145(:,2),'b.')
plot(x_1031,y_1031,'r-')
plot(x_1032,y_1032,'k-')
%plot(x_1145,y_1145,'b-')
legend('1031', '1032', 'location','southeast')
%title('Accumulated IT grade distribution - Machine number,  Thornton actual data - CPK 1.33')
xlabel('IT grade')
ylabel('Probability')
hold off
set(gcf, 'Color', 'w');

set(f3, 'units', 'centimeters', 'pos', [0 0 8 8])
export_fig('Acum_freqF3.pdf')
