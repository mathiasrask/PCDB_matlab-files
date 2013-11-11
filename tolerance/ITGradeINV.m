function [TOL] = ITGradeINV(ITG,Dim)

% This function recives two evenly long vecotrs of the nominal dimension
% in milimeters and tolerance in  micrometers

TOL = zeros(length(ITG),1);

for i = 1:length(ITG(:,1))
    if Dim >= 500
    TOL(i) =10^(0.2*(ITG(i)-1))*(0.004*Dim+2.1)/1000;
    else
    TOL(i) =10^(0.2*(ITG(i)-1))*(0.45*Dim^(1/3)+0.001*Dim)/1000;
    end
end

