classdef Utilities
    %Utilities Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
    end
    
    methods (Static)
       
        function [ITG] = dimAndTolToITGrade (dimension,tolerencewidth)

            % This function recives two evenly long vecotrs of the nominal dimension
            % in milimeters and tolerance in  micrometers
            
            vectorLength = length(tolerencewidth);
            ITG = zeros(vectorLength,1);
            
            for i = 1:vectorLength
                if dimension(i) > 500
                    ITG(i) =5*log10(tolerencewidth(i)*1000/(0.0004*dimension(i)+2.1))+1;
                else
                    ITG(i) =5*log10(tolerencewidth(i)*1000/(0.45*dimension(i)^(1/3)+dimension(i)*0.001))+1;
                end
            end
        end
        
        function [tolerancewidth] = dimAndITGradeToTol(dimension, ITgrade)

            % This function recives two evenly long vecotrs of the nominal dimension
            % in milimeters and tolerance in  micrometers
            
            
            vectorLength = length(ITgrade);
            tolerancewidth = zeros(vectorLength,1);
            
            for i = 1:vectorLength
                if dimension(i) >= 500
                tolerancewidth(i) =10^(0.2*(ITgrade(i)-1))*(0.004*dimension(i)+2.1)/1000;
                else
                tolerancewidth(i) =10^(0.2*(ITgrade(i)-1))*(0.45*dimension(i)^(1/3)+0.001*dimension(i))/1000;
                end
            end

        end
        
        function [PCSL] = meanshiftAndStdAndCpkToPCSL(meanshift, std, cpk)
            PCSL = 3*cpk.*std + abs(meanshift);
        end
    end
end

