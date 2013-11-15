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
            if length(dimension) < length(ITgrade)
                dimension = dimension * ones(length(ITgrade));
            elseif length(dimension) > length(ITgrade)
                ITgrade = ITgrade * ones(length(dimension));
            end
                
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
        
        function [lowerConfidenceLimit,midConfidenceLimit,upperConfidenceLimit,confidenreceIntevalWidth] = ConfidenceLimit (list_values)
            monte_std = std(list_values);
            monte_mean = mean(list_values);
            pd = makedist('Normal', 'mu', monte_mean, 'sigma', monte_std);

            lowerConfidenceLimit = icdf(pd, 0.025);
            midConfidenceLimit = icdf(pd, 0.5);
            upperConfidenceLimit = icdf(pd, 0.975);
            confidenreceIntevalWidth = upperConfidenceLimit-lowerConfidenceLimit;
        end
        
        function [xc, wsi] = wilson(x,n,alpha)
            %WILSON Centers and Wilson Score Intervals for binomial data.
            %   XC = WILSON(X,N) Returns the center of the Wilson Score Interval for the
            %   binomial distribution. X and N are scalars containing the number of
            %   successes and the number of trials, respectively.  If X and N are vectors,
            %   WILSON returns a vector of estimates whose I-th element is the parameter
            %   estimate for X(I) and N(I).  A scalar value for X or N is expanded to the
            %   same size as the other input.
            %
            %   [XC, WSI] = WILSON(X,N,ALPHA) gives the centers and 100(1-ALPHA) 
            %   percent Score Intervals given the data. Each row of WSI contains
            %   the lower and upper bounds for the corresponding element of XC.
            %   By default, the optional parameter ALPHA = 0.05 corresponding to 95%
            %   score interval.
            %
            %   See also BINOFIT. 

            %   Mike Sheppard
            %   MIT Lincoln Laboratory
            %   michael.sheppard@ll.mit.edu
            %   Original: 2-Jan-2011


            %The error catching is the same as BINOFIT. 
            %Included here, in slightly modified form, so the Statistics Toolbox is 
            %not required for the Wilson Score Interval function. Only the error
            %catching part of the code is used.
            %--------------
            if nargin < 3 
                alpha = 0.05;
            end
            % Initialize params to zero.
            [row, col] = size(x);
            if min(row,col) ~= 1
               error('WILSON:VectorRequired','First argument must be a vector.');
            end
            [r1,c1] = size(n);
            if ~isscalar(n)
               if row ~= r1 || col ~= c1
                  error('WILSON:InputSizeMismatch',...
                        'The first two inputs must match in size.');
               end
            end
            if ~isfloat(x)
               x = double(x);
            end
            if any(n<0) || any(n~=round(n)) || any(isinf(n)) || any(x>n)
                error('WILSON:InvalidN',...
                      'All N values must be non-negative integers at least as large as X.')
            end
            if any(x<0)
                error('WILSON:InvalidX','X must not be negative.')
            end
            %--------------



            %Compute Wilson Score Intervals
            phat = x./n;
            z=sqrt(2).*erfcinv(alpha);
            den=1+(z^2./n);
            xc=(phat+(z^2)./(2*n))./den;
            halfwidth=(z*sqrt((phat.*(1-phat)./n)+(z^2./(4*(n.^2)))))./den;
            wsi=[xc(:) xc(:)]+[-halfwidth(:) halfwidth(:)];


            end
    end
end

