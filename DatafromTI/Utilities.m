classdef Utilities
    %Utilities Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
    end
    
    methods (Static)
       %%
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
        %%
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
        %%
        function [PCSL] = meanshiftAndStdAndCpkToPCSL(meanshift, err_std, cpk)
            vectorlength = length(err_std);
            
            for i = 1:vectorlength
                PCSL(i) = 3*cpk*err_std(i) + abs(meanshift(i));
            end
        end
        %%
        function [lowerConfidenceLimit,midConfidenceLimit,upperConfidenceLimit,confidenreceIntevalWidth] = ConfidenceLimit (list_values)
            monte_std = std(list_values);
            monte_mean = mean(list_values);
            pd = makedist('Normal', 'mu', monte_mean, 'sigma', monte_std);

            lowerConfidenceLimit = icdf(pd, 0.025);
            midConfidenceLimit = icdf(pd, 0.5);
            upperConfidenceLimit = icdf(pd, 0.975);
            confidenreceIntevalWidth = upperConfidenceLimit-lowerConfidenceLimit;
        end
        %%
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
        %%
        function [cofidentWidth] = montecarloCLwidth (runs, sample_size, sample_sets, probability)
            
            % Generate data
            c_pk = 5/3; % 1.6667
            target = 100; % mm  (m)
            
            ITG_mean = 10;
            ITG_std = 1;
            C_a = 0.9;

            %ITGradesPD = makedist('Normal', 'mu', ITG_mean, 'sigma', ITG_std);
            
            % initialise
            sets_mean = zeros( sample_sets, 1);
            sets_std = zeros( sample_sets, 1);
            ITG_samplesets = zeros(sample_sets, 1);
            ITGatProperbility = zeros(sample_sets,length(probability));
            
            for i= 1:runs
               
                %ITG = icdf(ITGradesPD, rand(sample_sets,1));
                ITG = norminv(rand(sample_sets,1),ITG_mean,ITG_std)
                tolerenceWidth = Utilities.dimAndITGradeToTol(target, ITG);
                d = tolerenceWidth/2; % half specification width
                mu = (1-C_a) * d + target;
                sigma = (d - abs(mu - target) ) /(3*c_pk);

                for j = 1:sample_sets  

                    %sampleSetPD = makedist('Normal', 'mu', mu(j), 'sigma', sigma(j));
                    %samples = icdf(sampleSetPD, rand(sample_size,1));
                    samples = norminv(rand(sample_size,1),mu(j),sigma(j));
                    sample_std = std(samples);
                    sample_mean = mean(samples);
                    sample_meanshift = target - sample_mean;
                    sample_PCSL = Utilities.meanshiftAndStdAndCpkToPCSL(sample_meanshift,sample_std,c_pk);

                    ITG_samplesets(j) = Utilities.dimAndTolToITGrade(target,sample_PCSL*2);
                end

                sets_mean(i) = mean(ITG_samplesets);
                sets_std(i) = std(ITG_samplesets);

                %setsdist = makedist('Normal', 'mu', sets_mean(i), 'sigma', sets_std(i));

                for k = 1:length(probability)
                    %ITGatProperbility(i,k) = icdf(setsdist, probability(k));
                    ITGatProperbility(i,k) = norminv(probability(k),sets_mean(i),sets_std(i));
                end
            end 
            for i = 1:length(probability)
                [~,~,~,wid(i)] = Utilities.ConfidenceLimit(ITGatProperbility(:,i));
            end
            
            cofidentWidth = wid;
        end 
        
        %%
        function [x,y] = dataToAccumFrequncyPlot (lst_std, lst_mean)
          
            y = 0.001:0.01:0.999; %linspace(lst_min,lst_max,100);
            
            x = icdf('norm',y,lst_mean,lst_std);
            
        end
        %%
        function [s_std , s_mean, s_num, s_index, s_x, s_y] = listSorting (data, sort_lst)
            
            % Take a list of data and a sorting index vector
            % returns a standard deviation and mean of data of each group
            
            vectorlength = length(sort_lst);
            
            nandata = find(isnan(data)==1);
            
            sort_lst(nandata) = NaN;
            
            %
            indicies =[];
            for i= 1:vectorlength
                if isnan(sort_lst(i))
                elseif ismember(sort_lst(i),indicies)
                else
                    indicies = [indicies sort_lst(i)];
                end
            end
            
            indiciesLength = length(indicies);
            
             
            s = NaN(indiciesLength,vectorlength);
            s_y = NaN(indiciesLength,vectorlength);
           
            
            for i = 1:indiciesLength
                index = find(sort_lst == indicies(i));
                s_num(i) = length(index);
                s_ylin = linspace(0,1,s_num(i));
                for j=1:length(index)
                    s(i,j) = data(index(j));
                    s_y(i,j) = s_ylin(j);
                end
                s_x(i,:) = sort(s(i,:));
            end
           
           
            
            
            
            for i=1:indiciesLength
            s_std(i) = nanstd(s(i,:));
            s_mean(i) = nanmean(s(i,:));
            end
            
            s_index = indicies;
            
            
        end
    end
end

