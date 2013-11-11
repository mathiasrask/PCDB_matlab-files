
clc
clear all
close all


%% Generate Data

% Genrate 600 data series with smaples size between 1-30 samples with a
% variance correspondig to a ITgrade between 4-16.


    Proces = {'machining', 'casting', 'thermoshaping'};
    Material = {'Steel', 'Plastic', 'Aluminium', 'composite'};
    Type = {'A', 'R', 'L', 'P'};
    
    
    NumPoint = 600;

for i =1:NumPoint
    N = randi(30,1);
    mu = randi(800,1);
    
    
    ITG = (randi(120,1)+40)/10;
    if mu > 500
        T=10^(0.2*(ITG-1))*(0.004*mu+2.1);
    else 
        T=10^(0.2*(ITG-1))*(0.45*mu^(1/3)+mu*0.001);
    end
    
    std = T/3/1000; % 3 Sigmas
    
    Tag = [Proces(randi(3,1)) ; Material(randi(4,1))];
    
    Sample.(['A' num2str(i)]).Mean = mu;
    Sample.(['A' num2str(i)]).Raw = DataGeneratorFunction(N,std,mu);
    Sample.(['A' num2str(i)]).Tags = Tag;
    Sample.(['A' num2str(i)]).ID = i;
    Sample.(['A' num2str(i)]).Type = Type(randi(4,1));
   
end


%% Analyse Database

%sorting tags
TagID.str = {};
for i= 1:NumPoint
   
    tags = Sample.(['A' num2str(i)]).Tags;
    for j = 1:length(tags)
        ind=find(ismember(TagID.str,tags{j}));
        
        if ind > 0
            TagID.(tags{j}) = [TagID.(tags{j})  i ];
        else 
            TagID.(tags{j}) = i;
            TagID.str = [TagID.str;  tags{j}];
        end
    end
end


%calculating Variance and IT grade

for i = 1:NumPoint
    Dim = Sample.(['A' num2str(i)]).Mean;
    Var = var(Sample.(['A' num2str(i)]).Raw)^(1/2);
    ITG(i) = ITGrade(Dim,Var)';
end


% sort on size

sortSize = 1;
sortTags = 0;

j=1;
k=1;
l=1;


if sortSize
    for i = 1:NumPoint
        if Sample.(['A' num2str(i)]).Mean < 150
            y_150(j,1)=ITG(i);
            j = j+1;
        elseif Sample.(['A' num2str(i)]).Mean < 500
            y_500(k,1) = ITG(i);
            k= k+1;
        else
            y_3000(l,1) = ITG(i);
            l =l+1;
        end
    end
        
    y_150(:,2) = linspace(0,1,j-1);
    y_500(:,2) = linspace(0,1,k-1);
    y_3000(:,2) = linspace(0,1,l-1);

    y_150 = sort(y_150,1);
    y_500 = sort(y_500,1);
    y_3000 = sort(y_3000,1);
        
    
    hold on 
    plot(y_150(:,1),y_150(:,2),'r')
    plot(y_500(:,1),y_500(:,2),'k')
    plot(y_3000(:,1),y_3000(:,2),'b')
    legend('0-150mm','150-500 mm','500- mm','location','southeast')
end





