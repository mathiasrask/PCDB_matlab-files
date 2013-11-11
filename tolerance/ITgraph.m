clc
close all
clear all


ITG =linspace(5,9,5);

D = linspace(1,750);

T=zeros(length(D),length(ITG));


for i=1:length(ITG)
 for j=1:length(D)
    
    if D(j) > 500
    T(j,i)=10^(0.2*(ITG(i)-1))*(0.004*D(j)+2.1);
    else 
    T(j,i)=10^(0.2*(ITG(i)-1))*(0.45*D(j)^(1/3)+D(j)*0.001);
    end
    
 end
end


hold on
plot(D,T(:,:),'LineWidth',2)
legend('IT 5','IT 6', 'IT 7','IT 8','IT 9','location','northwest')

