clc 
close all
clear all


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


figure()
hold on
plot(D,T/1000)
plot(x,y,'r')
plot(x2,y2,'g')
legend('IT 13', 'NFT','DS812','location','southeast')
xlabel('Dimension')
ylabel('Tolerance')