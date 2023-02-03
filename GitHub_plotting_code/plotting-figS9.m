
clear;
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Variables/AOD_matfiles
load('Master_AOD_coarse_2007_2018.mat')
load('Master_AOD_fine_2007_2018.mat')
load('Master_AOD_total_2007_2018.mat')
load('Master_AOD_total_pixelcounts_2007_2018.mat')
load('Dimensions_AOD.mat')

%%
% Only select pixels with greater than or equal to X number of observations

pixel_threshold_2007_2018 = Master_AOD_total_pixelcounts_2007_2018 >= 2;
threshold_total_AOD_2007_2018 = Master_AOD_total_2007_2018;
threshold_coarse_AOD_2007_2018 = Master_AOD_coarse_2007_2018;
threshold_fine_AOD_2007_2018 = Master_AOD_fine_2007_2018; 

threshold_total_AOD_2007_2018(~pixel_threshold_2007_2018) = NaN;
threshold_coarse_AOD_2007_2018(~pixel_threshold_2007_2018) = NaN;
threshold_fine_AOD_2007_2018(~pixel_threshold_2007_2018) = NaN;

% Any negative values should also get thrown out (as NAN)

threshold_coarse_AOD_2007_2018(threshold_coarse_AOD_2007_2018<0) = NaN;
threshold_fine_AOD_2007_2018(threshold_fine_AOD_2007_2018<0) = NaN;
threshold_total_AOD_2007_2018(threshold_total_AOD_2007_2018<0) = NaN;



t1 = datetime(2007,01,01);
t2 = datetime(2018,12,31);
times_2007_2018 = t1:caldays(1):t2; 
times_2007_2018 = times_2007_2018';
clear t1 t2




[winter_x, winter_y] = find(times_2007_2018.Month >= 6 & times_2007_2018.Month <= 8);
% Can check if this worked by typing 'times(winter_x)' in command window
times_2007_2018(winter_x)
% 
AOD_total_winter = threshold_total_AOD_2007_2018; 
AOD_total_winter = AOD_total_winter(:,:,winter_x);
AOD_fine_winter = threshold_fine_AOD_2007_2018;
AOD_fine_winter = AOD_fine_winter(:,:,winter_x);
AOD_coarse_winter = threshold_coarse_AOD_2007_2018;
AOD_coarse_winter = AOD_coarse_winter(:,:,winter_x);

Pixel_Counts_winter = Master_AOD_total_pixelcounts_2007_2018(:,:,winter_x);

% For Spring: 
[spring_x, spring_y] = find(times_2007_2018.Month >= 9 & times_2007_2018.Month <= 11);
times_2007_2018(spring_x)
AOD_total_spring = threshold_total_AOD_2007_2018; 
AOD_total_spring = AOD_total_spring(:,:, spring_x); 
AOD_fine_spring = threshold_fine_AOD_2007_2018;
AOD_fine_spring = AOD_fine_spring(:,:,spring_x);
AOD_coarse_spring = threshold_coarse_AOD_2007_2018;
AOD_coarse_spring = AOD_coarse_spring(:,:,spring_x);

Pixel_Counts_spring = Master_AOD_total_pixelcounts_2007_2018(:,:,spring_x);

% For Summer: 

[summer_x, summer_y] = find(times_2007_2018.Month >= 1 & times_2007_2018.Month <=2 | times_2007_2018.Month == 12);

times_2007_2018(summer_x)
AOD_total_summer = threshold_total_AOD_2007_2018; 
AOD_total_summer = AOD_total_summer(:,:, summer_x); 
AOD_fine_summer = threshold_fine_AOD_2007_2018;
AOD_fine_summer = AOD_fine_summer(:,:,summer_x);
AOD_coarse_summer = threshold_coarse_AOD_2007_2018;
AOD_coarse_summer = AOD_coarse_summer(:,:,summer_x);


Pixel_Counts_summer = Master_AOD_total_pixelcounts_2007_2018(:,:,summer_x);

% For Fall: 

[fall_x, fall_y] = find(times_2007_2018.Month >= 3 & times_2007_2018.Month <= 5); 
times_2007_2018(fall_x) 
AOD_total_fall = threshold_total_AOD_2007_2018; 
AOD_total_fall = AOD_total_fall(:,:, fall_x); 
AOD_fine_fall = threshold_fine_AOD_2007_2018;
AOD_fine_fall = AOD_fine_fall(:,:,fall_x);
AOD_coarse_fall = threshold_coarse_AOD_2007_2018;
AOD_coarse_fall = AOD_coarse_fall(:,:,fall_x);

Pixel_Counts_fall = Master_AOD_total_pixelcounts_2007_2018(:,:, fall_x);

% 
AOD_total_winter_mean = mean(AOD_total_winter,3 ,'omitnan');
AOD_fine_winter_mean = mean(AOD_fine_winter, 3, 'omitnan');
AOD_coarse_winter_mean = mean(AOD_coarse_winter, 3, 'omitnan');
Pixel_Counts_winter_mean = mean(Pixel_Counts_winter, 3, 'omitnan');

AOD_total_spring_mean = mean(AOD_total_spring, 3 , 'omitnan'); 
AOD_fine_spring_mean = mean(AOD_fine_spring, 3, 'omitnan');
AOD_coarse_spring_mean = mean(AOD_coarse_spring, 3, 'omitnan');
Pixel_Counts_spring_mean = mean(Pixel_Counts_spring, 3, 'omitnan');

AOD_total_summer_mean = mean(AOD_total_summer, 3, 'omitnan'); 
AOD_fine_summer_mean = mean(AOD_fine_summer, 3, 'omitnan');
AOD_coarse_summer_mean = mean(AOD_coarse_summer, 3, 'omitnan');
Pixel_Counts_summer_mean = mean(Pixel_Counts_summer, 3, 'omitnan');

AOD_total_fall_mean = mean(AOD_total_fall, 3, 'omitnan'); 
AOD_fine_fall_mean = mean(AOD_fine_fall, 3, 'omitnan');
AOD_coarse_fall_mean = mean(AOD_coarse_fall, 3, 'omitnan');
Pixel_Counts_fall_mean = mean(Pixel_Counts_fall, 3, 'omitnan');

%%

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis

%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.02], [0.01 0.05], [0.1 0.1]);

if ~make_it_tight,  clear subplot;  end

fig = figure;clf;


ax(1) = subplot(2,4,1);
plot_BSea_figure_subplot_allvars(AOD_total_winter_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.15]); 
colormap(ax(1),cmocean('tempo'))

title({'Winter'}, 'FontSize',13);






ax(2) = subplot(2,4,2);
plot_BSea_figure_subplot_allvars(AOD_total_spring_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.15]); 
colormap(ax(2),cmocean('tempo'))

title({'Spring'}, 'FontSize',13);


ax(3) = subplot(2,4,3);
plot_BSea_figure_subplot_allvars(AOD_total_summer_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.15]); 
colormap(ax(3),cmocean('tempo'))

title({'Summer'}, 'FontSize',13);


ax(4) = subplot(2,4,4);
plot_BSea_figure_subplot_allvars(AOD_total_fall_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.15]); 
colormap(ax(4),cmocean('tempo'))
title({'Fall'}, 'FontSize',13);

    
ax(5) = subplot(2,4,5);
plot_BSea_figure_subplot_allvars(AOD_fine_winter_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.15]); 
colormap(ax(5),cmocean('tempo'))




ax(6) = subplot(2,4,6);
plot_BSea_figure_subplot_allvars(AOD_fine_spring_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.15]); 
colormap(ax(6),cmocean('tempo'))

ax(7) = subplot(2,4,7);
plot_BSea_figure_subplot_allvars(AOD_fine_summer_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.15]); 
colormap(ax(7),cmocean('tempo'))

ax(8) = subplot(2,4,8);
plot_BSea_figure_subplot_allvars(AOD_fine_fall_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.15]); 
colormap(ax(8),cmocean('tempo'))

hp4 = get(subplot(2,4,4),'Position');
h = colorbar('Position', [0.93 0.2  0.015  0.6]);
h.FontWeight = 'bold';
h.FontSize = 12;
hy = ylabel(h, sprintf('AOD'), 'FontSize', 12);
% hy.FontSize = 18;

%%
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures
set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'AODT_AODF_climatologies','-dpng','-r300');       %  *// 300 dpi

