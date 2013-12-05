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

for i =1:25
    count = count +1;
    Data(count,:) = [rotor_kina(i,1:5),NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end
for i =1:41
    count = count +1;
    Data(count,:) = [axel_kina(i,1:5),NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end
for i =1:14
    count = count +1;
    Data(count,:) = [helper_kina(i,1:5),NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end
for i =1:10
    count = count +1;
    Data(count,:) = [insert_kina(i,1:5),NaN,NaN,NaN,NaN,NaN,NaN,NaN];
end


List_Target = [rotor_spec(:,1);rotor_spec(:,1);rotor_spec(:,1);axel_spec(:,1);axel_spec(:,1);rotor_kina_spec(:,1);axel_kina_spec(:,1);helper_kina_spec(:,1);insert_kina_spec(:,1)];
List_SpecUp = [rotor_spec(:,2);rotor_spec(:,2);rotor_spec(:,2);axel_spec(:,2);axel_spec(:,2);rotor_kina_spec(:,2);axel_kina_spec(:,2);helper_kina_spec(:,2);insert_kina_spec(:,2)];

cpk = 1.66;

% Data sets sort

rotor_k_set = [1:25];
rotor_set = [1,1,1,1, 2,2, 3,3,3,3, 4,4, 5, 6,6, 7,7, 8, 9,9, 10,10, 11, 12,12,12,12, 13,13,14,14,15,15,16,16]+25;

axel_k_set = [1:9,9:22,22:39]+25+16;
axel_set = [1:11]+25+16+39;

helper_set = [1:14]+25+16+39+11;
insert_set = [1:10]+25+16+39+11+14;

List_sets = [rotor_set,rotor_set,rotor_set,axel_set,axel_set,rotor_k_set,axel_k_set,helper_set,insert_set];

% Aggregation of data

sets_max = max(List_sets);

Agg_control = 1;

if Agg_control
    for i=1:sets_max
        set_index = find(List_sets==i);
        
        set_target(i) = List_Target(set_index(1));
        set_specup(i) = List_SpecUp(set_index(1));
        
        lst_data = [];
        for j=1:length(set_index)
            lst_data = [lst_data, Data(set_index(j),:)];
        end

        s_size = length(find(isnan(lst_data)==0));
        biasCorrectionFactor_c4 = sqrt(2/(s_size - 1)) * gamma(s_size/2)/gamma((s_size-1)/2);
        List_Std(i) = nanstd(lst_data) / biasCorrectionFactor_c4;
        List_Mean(i) = nanmean(lst_data);
        List_MeanShift(i) = List_Mean(i) - set_target(i);
    end    
    List_Target = set_target;
    List_SpecUp = set_specup';
    
    uglyvector = [28,89,25,97,98,99,100,101,102,103,104,105,110,];
    
else
    for i=1:length(Data(:,1))
        s_size = length(find(isnan(Data(i,:))==0));
        biasCorrectionFactor_c4 = sqrt(2/(s_size - 1)) * gamma(s_size/2)/gamma((s_size-1)/2);
        List_Std(i) = nanstd(Data(i,:)) / biasCorrectionFactor_c4;
        List_Mean(i) = nanmean(Data(i,:));
        List_MeanShift(i) = List_Mean(i) - List_Target(i); 
    end
    uglyvector = [7,8,9,10,42,43,44,45,77,78,79,80,114,125];
end

% removing ugly data

List_Target(uglyvector) = NaN;
List_MeanShift(uglyvector) = NaN;
List_Std(uglyvector) = NaN;

% Caculating ITG
List_PCSL = Utilities.meanshiftAndStdAndCpkToPCSL(abs(List_MeanShift), List_Std, cpk);
List_ITG = Utilities.dimAndTolToITGrade(abs(List_Target) ,List_PCSL*2);
List_SpecITG = Utilities.dimAndTolToITGrade(abs(List_Target) ,List_SpecUp*2);

List_Ca = 1 - abs(List_MeanShift)./List_PCSL;
List_Ca_spec = 1 - abs(List_MeanShift)./List_SpecUp';
List_Ca_spec2 = List_MeanShift./List_SpecUp';

List_Cp_spec = List_SpecUp'./ (3 * List_Std);
List_Cp = List_PCSL./ (3 * List_Std);

% Sorting for Dimension
[num, center] = hist(List_Target,3);

for i = 1:length(center)
    low_lim(i) = center(i) - (center(2) - center(1))/2;
    hig_lim(i) = center(i) + (center(2) - center(1))/2;
end

List_SortDim = zeros(0,length(List_Target));

low_lim = [0,5,10];

for j = 1:length(List_Target)
    for i = 1:length(low_lim)
        if low_lim(i) < List_Target(j) 
            List_SortDim(j) = i;
        end
    end         
end

% Sorting Material

List_material = [ones(1,105),2*ones(1,22),ones(1,25),2*ones(1,(41+14+10))];



% Sorting for Molds

rotor_mold = [1,1,1,1, 1,1, NaN,NaN,NaN,NaN, 1,1, 2, 1,1, 1,1, 1, 1,1, 2,2, 1, 1,1,1,1, 1,1, 2,2, 1,1,2,2];
axel_mold = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];

rotor_kina_mold = [1,1,1,1,1,1,1,1,2,1,1,2,1,1,1,1,1,1,1,2,1,1,2,1,1];
axel_kina_mold = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,1,1,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1];
helper_kina_mold = [2,1,1,1,1,1,1,1,1,1,1,1,1,1];
insert_kina_mold = [1,1,1,2,1,1,1,2,1,1];

List_Mold = [rotor_mold, rotor_mold, rotor_mold, axel_mold, axel_mold,rotor_kina_mold,axel_kina_mold,helper_kina_mold,insert_kina_mold];

% sorting inside, outside and translation
% 1 = inside
% 2 = outside
% 3 = translational surfaces

rotor_inside = [1,1,1,1, 1,1, NaN,NaN,NaN,NaN, 1,1, 2, 2,2, 1,1, 1, 1,1, NaN,NaN, 1, NaN,NaN,NaN,NaN, 2,2, 2,2, 2,2, 2,2 ];
axel_inside = [2, NaN, 2, 2, 2, 2, 2, 2, 2, 2, 2];

rotor_k_inside = [1,1,1,NaN,2,NaN,1,1,2,2,1,NaN,NaN,1,2,NaN,NaN,NaN,2,2,2,2,2,1,NaN];
axel_k_inside = [1,1,NaN,NaN,2,NaN,NaN,NaN,NaN,NaN,NaN,2,NaN,NaN,NaN,2,2,2,2,2,2,2,2,2,NaN,2,2,2,1,NaN,NaN,2,2,NaN,NaN,2,2,NaN,NaN,NaN,NaN];
helper_inside = [1,NaN,NaN,2,NaN,NaN,2,2,NaN,NaN,2,NaN,NaN,NaN];
insert_inside = [1,1,1,2,NaN,1,2,2,2,2];

List_inside = [rotor_inside, rotor_inside, rotor_inside, axel_inside, axel_inside, rotor_k_inside, axel_k_inside, helper_inside, insert_inside];


% sort kina

List_kina = [zeros(1,127),ones(1,25+41+14+10)];

if Agg_control
    for i=1:sets_max
        set_index = find(List_sets==i);
        
        lst_inside(i) = List_inside(set_index(1));
        lst_kina(i) = List_kina(set_index(1));
        lst_mold(i) = List_Mold(set_index(1));
        lst_material(i) = List_material(set_index(1));
        
    end
    List_inside = lst_inside;
    List_kina = lst_kina;
    List_Mold = lst_mold;
    List_material = lst_material;
end


% sorting of radius and diameters

rotor_rad = [0,0,0,0, 0,0, NaN,NaN,NaN,NaN, 0,0, 0, 2,2, 0,0, 0, 2,2, 0,0, 2, 0,0,0,0, 0,0, 0,0, 0,0, 0,0 ];

axel_rad = [0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2];

List_rad = [rotor_rad, rotor_rad, rotor_rad, axel_rad, axel_rad];

% sorting for color

List_color = [ones(35,1);2*ones(35,1);3*ones(35,1);NaN(22,1)];

%Total Plot

List_ALL = [ones(25,1); ...
            ones(16,1); ...
            ones(39,1); ...
            ones(11,1); ...
            ones(14,1); ...
            ones(10,1)];

% COMPARISON OF SORTING LISTS

% Plot Ca .. sigma/d
for run = 1:0
    sig_d = List_Std' ./ List_SpecUp;
    
    cpm = 1.66;
    x_cpm = linspace(0,1,40);
    for i=1:40
    y_cpm(i) = (1/3)*sqrt(-9*x_cpm(i)^2*cpm^2+18*x_cpm(i)*cpm^2-9*cpm^2+1)/cpm;
    end

    f2 = figure();
    hold on
    plot(List_Ca_spec,sig_d,'.')
    plot([0,1],[0,1/(3*1.66)],'r')
    plot(x_cpm,y_cpm,'g')
    
    
    xlabel('C_{a} closeness to target')
    ylabel('\sigma / d')
    legend('Data','C_{pk}','C_{pm}')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('Ca_sigma_d.pdf')  
end

% Plot Ca .. Cp
for run = 1:0

    
    
    f2 = figure();
    hold on
    plot(List_Ca_spec,List_Cp_spec,'.')
    
    xlabel('C_{a} closeness to target')
    ylabel('C_{p} Precision index')
    legend('Data')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('Ca_cp.pdf')    
end

% Plot ITG .. total
for run = 1
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_ITG, List_ALL);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf,'r');
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    xlim([7,16])
    
    xlabel('Process capability driven tolerance (IT Grade)')
    ylabel('Probability')
    %legend('fitted curve','Data','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('ITG_total.pdf')    
end

% Plot Cp .. total
for run = 1:0
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_Cp, List_ALL);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    %plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    
    xlabel('C_{p} Precision index')
    ylabel('Probability')
    %legend('fitted curve','Data','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('Cp_total.pdf')    
end

% Plot ITG .. Sort Dimension
for run = 1:0
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_ITG, List_SortDim);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    
    xlabel('Precision (IT Grade)')
    ylabel('Probability')
    legend('0-5 mm','5-10 mm','10- mm','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('ITG_Dimsort.pdf')    
end

% Plot ITG .. material
for run = 1:0
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_ITG, List_material);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    
    xlabel('Precision (IT Grade)')
    ylabel('Probability')
    legend('ABS+PC','POM','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('ITG_material.pdf')    
end

% Plot cb .. material
for run = 1:0
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_Ca_spec2, List_material);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    
    xlabel('C_{b} meanshift / specified tolerance width')
    ylabel('Probability')
    legend('ABS+PC','POM','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('Cb_material.pdf')    
end

% Plot ITG .. Mold lines
for run = 1:0
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_ITG, List_Mold);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    
    xlabel('Precision (IT Grade)')
    ylabel('Probability')
    legend('in mold','across moldline','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('ITG_moldline.pdf')    
end

% Plot Cb .. Mold lines
for run = 1:0
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_Ca_spec2, List_Mold);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    
    xlabel('C_{b} meanshift / specified tolerance width')
    ylabel('Probability')
    legend('in mold','across moldline','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('Cb_moldline.pdf')    
end

% Plot Ca .. Mold lines
for run = 1:0
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_Ca_spec, List_Mold);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    
    xlabel('C_{a} closeness to target')
    ylabel('Probability')
    legend('in mold','across moldline','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('Ca_moldline.pdf')    
end

% Plot Cb .. Hole and Shaft
for run = 1:0
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_Ca_spec2, List_inside);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    
    xlabel('C_{b} meanshift / specified tolerance width')
    ylabel('Probability')
    legend('Hole','Shaft','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('Cb_holeshaft.pdf')    
end

% Plot ITG .. Hole and Shaft
for run = 1:0
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_ITG, List_inside);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    
    xlabel('Precision (IT Grade)')
    ylabel('Probability')
    legend('Hole','Shaft','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('ITG_holeshaft.pdf')    
end

% Plot ITG .. Kina Measurements vs. DTI
for run = 1:0

    
    
    [s_std, s_mean, s_num, s_index, s_x, s_y] = Utilities.listSorting(List_ITG, List_kina);

    for i=1:length(s_std)

        s_min = s_mean(i) - s_std(i)*3;
        s_max = s_mean(i) + s_std(i)*3;
        x(:,i) = linspace(s_min,s_max,100);

        PD = makedist('Normal', 'mu', s_mean(i), 'sigma', s_std(i));
        xcdf(:,i) = cdf(PD,x(:,i));

        factor = 1;

        [xc(:,i), wsi(:,:,i)] = Utilities.wilson(xcdf(:,i)*floor(s_num(i)/factor),floor(s_num(i)/factor),0.05);
    end
    f2 = figure();
    hold on
    %plot(x,xc, 'b-.')
    plot(x,xcdf);
    %plot(x, wsi(:,1), '-.')
    %plot(x, wsi(:,2), '-.')
    plot(s_x',s_y','.')
    
    xlabel('Precision (IT grade)')
    ylabel('Probability')
    legend('External control','Production','location','southeast')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('ITG_china.pdf')   
end

% PLOT OF DIMENSION VS. TOLERANCE
for run = 1:0
    
    % ISO line
    dimension = [2:1:80];
    ITgrade = 10:1:14;
    for i = 1:length(ITgrade)
    tol(i,:) = Utilities.dimAndITGradeToTol(dimension, ITgrade(i));
    end
   
    f2 = figure ();
    hold on 
    plot(abs(List_Target),abs(List_PCSL*2),'.')
    plot(dimension,tol)
    xlabel('Dimension')
    ylabel('Tolerance')
    legend('Data','ITG10','ITG11','ITG12','ITG13','ITG14')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('DimTol.pdf')
end

% NORMAL PLOT
for run = 1:0
    f2 = figure ();
    normplot(List_ITG)'
    xlabel('IT Grade')
    title('')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])
    
    export_fig('Normal_plot.pdf')

end

% IT GRADE vs. Specified IT GRADE
for run = 1
    f2 = figure();
    hold on 
    box off
    
    plot(List_SpecITG,List_ITG, '.')
    plot([7.5,17],[7.5,17],'r')
    plot([7.5,17],[14.22, 14.22],'--','color',[0 0.5 0])
    axis([7.5 17 7.5 17])
    ylabel('Process capability driven tolerance (IT-grade)')
    xlabel('Specified Tolerance (IT-grade)')
    set(gcf, 'Color', 'w');
    set(f2, 'units', 'centimeters', 'pos', [0 0 8 8])

    export_fig('ITG_ITGSpec.pdf')
end
