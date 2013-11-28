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

load 'DATA.mat'


% Assembling data
count = 0;
for i =1:3
    for j=1:35
        count = count +1;
        Data(count,:) = rotor(j,(1:12)+(i*12-12));
    end
end

for i =1:11
    count = count +1;
    Data(count,:) = [axel(i,1:8) , NaN, NaN, NaN, NaN];
end

for i =1:11
    count = count +1;
    Data(count,:) = [axel(i,9:19) , NaN];
end
% Caculating ITG

List_Target = [rotor_spec(:,1);rotor_spec(:,1);rotor_spec(:,1);axel_spec(:,1);axel_spec(:,1)];
List_SpecUp = [rotor_spec(:,2);rotor_spec(:,2);rotor_spec(:,2);axel_spec(:,2);axel_spec(:,2)];

cpk = 1.66;
for i=1:length(Data(:,1))
    List_Std(i) = nanstd(Data(i,:));
    List_Mean(i) = nanmean(Data(i,:));
    List_MeanShift(i) = List_Mean(i) - List_Target(i); 
end

List_PCSL = Utilities.meanshiftAndStdAndCpkToPCSL(abs(List_MeanShift), List_Std, cpk);
List_ITG = Utilities.dimAndTolToITGrade(abs(List_Target) ,List_PCSL*2);
List_SpecITG = Utilities.dimAndTolToITGrade(abs(List_Target) ,List_SpecUp);

List_Ca = 1 - abs(List_MeanShift)./List_PCSL;

% removing ugly data

uglyvector = [7,8,9,10,42,43,44,45,77,78,79,80,114,125];

List_ITG(uglyvector) = NaN;
List_PCSL(uglyvector) = NaN;

% Sorting for Dimension
[num, center] = hist(List_Target,10);

for i = 1:length(center)
    low_lim(i) = center(i) - (center(2) - center(1))/2;
    hig_lim(i) = center(i) + (center(2) - center(1))/2;
end

List_SortDim = zeros(0,length(List_Target));


for j = 1:length(List_Target)
    for i = 1:length(center)
        if low_lim(i) < List_Target(j) 
            List_SortDim(j) = i;
        end
    end         
end

% Sorting Material

List_material = [ones(1,105),2*ones(1,22)];

% Sorting for Molds

rotor_mold = [1,1,1,1, 1,1, NaN,NaN,NaN,NaN, 1,1, 23, 2,3, 1,1, 1, 2,3, 13,12, 1, 1,1,1,1, 2,3, 23,23, 2,3, 23,23];

axel_mold = [45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45];

List_Mold = [rotor_mold, rotor_mold, rotor_mold, axel_mold, axel_mold];

% sorting inside, outside and translation
% 1 = inside
% 2 = outside
% 3 = translational surfaces

rotor_inside = [1,1,1,1, 1,1, NaN,NaN,NaN,NaN, 1,1, 2, 2,2, 1,1, 1, 1,1, 3,3, 1, 3,3,3,3, 2,2, 2,2, 2,2, 2,2 ];

axel_inside = [2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2];

List_inside = [rotor_inside, rotor_inside, rotor_inside, axel_inside, axel_inside];

% sorting of radius and diameters

rotor_rad = [0,0,0,0, 0,0, NaN,NaN,NaN,NaN, 0,0, 0, 2,2, 0,0, 0, 2,2, 0,0, 2, 0,0,0,0, 0,0, 0,0, 0,0, 0,0 ];

axel_rad = [0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2];

List_rad = [rotor_rad, rotor_rad, rotor_rad, axel_rad, axel_rad];

% sorting for color

List_color = [ones(35,1);2*ones(35,1);3*ones(35,1);NaN(22,1)];


% COMPARISON OF SORTING LISTS
if 0

    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_Ca, List_inside);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    figure()
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
end



% PLOT OF DIMENSION VS. TOLERANCE
if 1
    
    % ISO line
    dimension = [2:1:80];
    ITgrade = 10:1:14;
    for i = 1:length(ITgrade)
    tol(i,:) = Utilities.dimAndITGradeToTol(dimension, ITgrade(i));
    end
   
    figure ()
    hold on 
    plot(abs(List_Target),abs(List_PCSL*2),'.')
    plot(dimension,tol,'r')
    
end

if 1
    figure ()
    normplot(List_ITG)
    
end
