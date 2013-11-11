clc 
close all
clear all

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 9)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 9)


addpath('lib/export_fig');


D = linspace(1,200,100)
ITG = [13]

for i=1:length(ITG)
 for j=1:length(D)
    
    if D(j) > 500
    T(j,i)=10^(0.2*(ITG(i)-1))*(0.004*D(j)+2.1);
    else 
    T(j,i)=10^(0.2*(ITG(i)-1))*(0.45*D(j)^(1/3)+D(j)*0.001);
    end
    
 end
end

x = [1,3,6,10,15,22,30,40,53,70,90,115,150,200]
y = [0.13, 0.15,0.17,0.20,0.22,0.25,0.28,0.32,0.37,0.44,0.50,0.60,0.75,0.95]


x2 = [1,3,6,10,15,22,30,40,53,70,90,120,160,200]
y2 = [0.12,0.14,0.16,0.18,0.2,0.22,0.26,0.3,0.34,0.4,0.48,0.58,0.7,0.86]

x3 = [0.5,3,6,30,120,400]
y3 = [0.1,0.1,0.2,0.3,0.5,0.8]

x4 = [3,6,10,18,30,50,80,120,180,250]
y4 = [0.14,0.18,0.22,0.27,0.33,0.39,0.46,0.54,0.63,0.72]

f2 = figure()
hold on
plot(D,T/1000)
plot(x,y,'r')
plot(x2,y2,'g')
plot(x3,y3,'k')
%plot(x4,y4,'m')
axis([0 200 0 1.2])
legend('ANSI B4-2', 'NFT58000','DS812','DIN7168','ISO286-2','location','southeast')
xlabel('Dimension')
ylabel('Tolerance')
set(gcf, 'Color', 'w');
set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

export_fig('Tolerance_standards.pdf')
