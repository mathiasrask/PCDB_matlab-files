function [ITG] = ITGrade (Dim,Var)

% This function recives two evenly long vecotrs of the nominal dimension
% in milimeters and tolerance in  micrometers


for i = 1:length(Var(:,1))
    if Dim(i) > 500
    ITG(i) =5*log10(2*Var(i)*1000/(0.0004*Dim(i)+2.1))+1;
    else
    ITG(i) =5*log10(2*Var(i)*1000/(0.45*Dim(i)^(1/3)+Dim(i)*0.001))+1;
    end
end

