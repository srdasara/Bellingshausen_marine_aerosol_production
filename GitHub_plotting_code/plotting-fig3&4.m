
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_vars

% Scatterplot Code:
clear;
load('Southernmost_complete_TT_daily.mat')
load('Corner_complete_TT_daily.mat')
load('Central_complete_TT_daily.mat')

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures/Diurnal_winds

%% Also checking winter vs summer as in Dror et al 2018

% Corner spatial region, daily, subsetting
Winter_TR_Corner_daily_index = month(Corner_complete_TT_daily.Total_Profile_Time_New_Surface) >=6 &...
    month(Corner_complete_TT_daily.Total_Profile_Time_New_Surface) <=8;

Spring_TR_Corner_daily_index = month(Corner_complete_TT_daily.Total_Profile_Time_New_Surface) >=9 &...
    month(Corner_complete_TT_daily.Total_Profile_Time_New_Surface) <=11;

Summer_TR_Corner_daily_index = month(Corner_complete_TT_daily.Total_Profile_Time_New_Surface) >=1 &...
    month(Corner_complete_TT_daily.Total_Profile_Time_New_Surface) <=2 | ...
    month(Corner_complete_TT_daily.Total_Profile_Time_New_Surface) == 12;

Fall_TR_Corner_daily_index = month(Corner_complete_TT_daily.Total_Profile_Time_New_Surface) >=3 &...
    month(Corner_complete_TT_daily.Total_Profile_Time_New_Surface) <=5;

Winter_TT_Corner_BSea_daily = Corner_complete_TT_daily(Winter_TR_Corner_daily_index,:);
Spring_TT_Corner_BSea_daily = Corner_complete_TT_daily(Spring_TR_Corner_daily_index,:);
Summer_TT_Corner_BSea_daily = Corner_complete_TT_daily(Summer_TR_Corner_daily_index,:);
Fall_TT_Corner_BSea_daily  = Corner_complete_TT_daily(Fall_TR_Corner_daily_index,:);


% Central spatial region, daily, subsetting
Winter_TR_Central_daily_index = month(Central_complete_TT_daily.Total_Profile_Time_New_Surface) >=6 &...
    month(Central_complete_TT_daily.Total_Profile_Time_New_Surface) <=8;

Spring_TR_Central_daily_index = month(Central_complete_TT_daily.Total_Profile_Time_New_Surface) >=9 &...
    month(Central_complete_TT_daily.Total_Profile_Time_New_Surface) <=11;

Summer_TR_Central_daily_index = month(Central_complete_TT_daily.Total_Profile_Time_New_Surface) >=1 &...
    month(Central_complete_TT_daily.Total_Profile_Time_New_Surface) <=2 | ...
    month(Central_complete_TT_daily.Total_Profile_Time_New_Surface) == 12;

Fall_TR_Central_daily_index = month(Central_complete_TT_daily.Total_Profile_Time_New_Surface) >=3 &...
    month(Central_complete_TT_daily.Total_Profile_Time_New_Surface) <=5;

Winter_TT_Central_BSea_daily = Central_complete_TT_daily(Winter_TR_Central_daily_index,:);
Spring_TT_Central_BSea_daily = Central_complete_TT_daily(Spring_TR_Central_daily_index,:);
Summer_TT_Central_BSea_daily = Central_complete_TT_daily(Summer_TR_Central_daily_index,:);
Fall_TT_Central_BSea_daily  = Central_complete_TT_daily(Fall_TR_Central_daily_index,:);


% Southernmost spatial region, daily, subsetting
Winter_TR_Southernmost_daily_index = month(Southernmost_complete_TT_daily.Total_Profile_Time_New_Surface) >=6 &...
    month(Southernmost_complete_TT_daily.Total_Profile_Time_New_Surface) <=8;

Spring_TR_Southernmost_daily_index = month(Southernmost_complete_TT_daily.Total_Profile_Time_New_Surface) >=9 &...
    month(Southernmost_complete_TT_daily.Total_Profile_Time_New_Surface) <=11;

Summer_TR_Southernmost_daily_index = month(Southernmost_complete_TT_daily.Total_Profile_Time_New_Surface) >=1 &...
    month(Southernmost_complete_TT_daily.Total_Profile_Time_New_Surface) <=2 | ...
    month(Southernmost_complete_TT_daily.Total_Profile_Time_New_Surface) == 12;

Fall_TR_Southernmost_daily_index = month(Southernmost_complete_TT_daily.Total_Profile_Time_New_Surface) >=3 &...
    month(Southernmost_complete_TT_daily.Total_Profile_Time_New_Surface) <=5;

Winter_TT_Southernmost_BSea_daily = Southernmost_complete_TT_daily(Winter_TR_Southernmost_daily_index,:);
Spring_TT_Southernmost_BSea_daily = Southernmost_complete_TT_daily(Spring_TR_Southernmost_daily_index,:);
Summer_TT_Southernmost_BSea_daily = Southernmost_complete_TT_daily(Summer_TR_Southernmost_daily_index,:);
Fall_TT_Southernmost_BSea_daily  = Southernmost_complete_TT_daily(Fall_TR_Southernmost_daily_index,:);

%% Scatter plot, All Regions, daily data:
% SUBPLOTTED

% DAILY:
% Coarse AOD & MAOD

%%
clear x y nanVals p yfit R
x_Corner_Wind_Daily_forAODc = Corner_complete_TT_daily.MASTER_Winds_daytime;
y_Corner_AODc_Daily = Corner_complete_TT_daily.Corner_AODc_openocean;

nanVals_Corner_AODc_Wind = ismissing(x_Corner_Wind_Daily_forAODc) | ismissing(y_Corner_AODc_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_Daily_forAODc = x_Corner_Wind_Daily_forAODc(~nanVals_Corner_AODc_Wind);
y_Corner_AODc_Daily = y_Corner_AODc_Daily(~nanVals_Corner_AODc_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_Daily_AODc_Wind = polyfit(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_Daily_AODc_Wind = polyval(p_Corner_Daily_AODc_Wind,x_Corner_Wind_Daily_forAODc);

[R_Corner_Daily_AODc_Wind, pval_Corner_Daily_AODc_Wind] = corrcoef(x_Corner_Wind_Daily_forAODc, y_Corner_AODc_Daily );


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_Wind_Daily_forAODc = Central_complete_TT_daily.MASTER_Winds_daytime;
y_Central_AODc_Daily = Central_complete_TT_daily.Central_AODc_openocean;

nanVals_Central_AODc_Wind = ismissing(x_Central_Wind_Daily_forAODc) | ismissing(y_Central_AODc_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_Daily_forAODc = x_Central_Wind_Daily_forAODc(~nanVals_Central_AODc_Wind);
y_Central_AODc_Daily = y_Central_AODc_Daily(~nanVals_Central_AODc_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_Daily_AODc_Wind = polyfit(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_Daily_AODc_Wind = polyval(p_Central_Daily_AODc_Wind,x_Central_Wind_Daily_forAODc);

[R_Central_Daily_AODc_Wind, pval_Central_Daily_AODc_Wind] = corrcoef(x_Central_Wind_Daily_forAODc, y_Central_AODc_Daily );


% Southernmost daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_Wind_Daily_forAODc = Southernmost_complete_TT_daily.MASTER_Winds_daytime;
y_Southernmost_AODc_Daily = Southernmost_complete_TT_daily.Southernmost_AODc_openocean;

nanVals_Southernmost_AODc_Wind = ismissing(x_Southernmost_Wind_Daily_forAODc) | ismissing(y_Southernmost_AODc_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_Daily_forAODc = x_Southernmost_Wind_Daily_forAODc(~nanVals_Southernmost_AODc_Wind);
y_Southernmost_AODc_Daily = y_Southernmost_AODc_Daily(~nanVals_Southernmost_AODc_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_Daily_AODc_Wind = polyfit(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_Daily_AODc_Wind = polyval(p_Southernmost_Daily_AODc_Wind,x_Southernmost_Wind_Daily_forAODc);

[R_Southernmost_Daily_AODc_Wind, pval_Southernmost_Daily_AODc_Wind] = corrcoef(x_Southernmost_Wind_Daily_forAODc, y_Southernmost_AODc_Daily );

%%


clear x y nanVals p yfit R
x_Corner_Wind_Daily_forMAOD = Corner_complete_TT_daily.MASTER_Winds_nighttime;
y_Corner_MAOD_Daily = Corner_complete_TT_daily.MAOD_nocleanair_openocean;

nanVals_Corner_MAOD_Wind = ismissing(x_Corner_Wind_Daily_forMAOD) | ismissing(y_Corner_MAOD_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_Daily_forMAOD = x_Corner_Wind_Daily_forMAOD(~nanVals_Corner_MAOD_Wind);
y_Corner_MAOD_Daily = y_Corner_MAOD_Daily(~nanVals_Corner_MAOD_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_Daily_MAOD_Wind = polyfit(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_Daily_MAOD_Wind = polyval(p_Corner_Daily_MAOD_Wind,x_Corner_Wind_Daily_forMAOD);

[R_Corner_Daily_MAOD_Wind, pval_Corner_Daily_MAOD_Wind] = corrcoef(x_Corner_Wind_Daily_forMAOD, y_Corner_MAOD_Daily );


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_Wind_Daily_forMAOD = Central_complete_TT_daily.MASTER_Winds_nighttime;
y_Central_MAOD_Daily = Central_complete_TT_daily.MAOD_nocleanair_openocean;

nanVals_Central_MAOD_Wind = ismissing(x_Central_Wind_Daily_forMAOD) | ismissing(y_Central_MAOD_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_Daily_forMAOD = x_Central_Wind_Daily_forMAOD(~nanVals_Central_MAOD_Wind);
y_Central_MAOD_Daily = y_Central_MAOD_Daily(~nanVals_Central_MAOD_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_Daily = polyfit(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_Daily = polyval(p_Central_Daily,x_Central_Wind_Daily_forMAOD);

[R_Central_Daily_MAOD_Wind, pval_Central_Daily_MAOD_Wind] = corrcoef(x_Central_Wind_Daily_forMAOD, y_Central_MAOD_Daily );


% Southernmost daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_Wind_Daily_forMAOD = Southernmost_complete_TT_daily.MASTER_Winds_nighttime;
y_Southernmost_MAOD_Daily = Southernmost_complete_TT_daily.MAOD_nocleanair_openocean;

nanVals_Southernmost_MAOD_Wind = ismissing(x_Southernmost_Wind_Daily_forMAOD) | ismissing(y_Southernmost_MAOD_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_Daily_forMAOD = x_Southernmost_Wind_Daily_forMAOD(~nanVals_Southernmost_MAOD_Wind);
y_Southernmost_MAOD_Daily = y_Southernmost_MAOD_Daily(~nanVals_Southernmost_MAOD_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_Daily_MAOD_Wind = polyfit(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_Daily_MAOD_Wind = polyval(p_Southernmost_Daily_MAOD_Wind,x_Southernmost_Wind_Daily_forMAOD);

[R_Southernmost_Daily_MAOD_Wind, pval_Southernmost_Daily_MAOD_Wind] = corrcoef(x_Southernmost_Wind_Daily_forMAOD, y_Southernmost_MAOD_Daily );


%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

%Corner Daily Scatterplot
subplot(2,3,1)
scatter(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_Daily_forMAOD,yfit_Corner_Daily_MAOD_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 30])
set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)

title('Region A') 

txt = ['y = ' num2str(round(p_Corner_Daily_MAOD_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_MAOD_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_MAOD_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_MAOD_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_Daily_forMAOD),3,'significant'))];
Region = ('MAOD');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10, 0.36, p_txt,'FontSize', 15)
text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region



% Corner Daily Scatterplot:
subplot(2,3,4)
scatter(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_Daily_forAODc,yfit_Corner_Daily_AODc_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 30])
set(gca,'FontSize', 15)


txt = ['y = ' num2str(round(p_Corner_Daily_AODc_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_AODc_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_AODc_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_AODc_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_Daily_forAODc),3,'significant'))];
Region = ('AOD_C');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10, 0.36, p_txt,'FontSize', 15)
text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% Central Daily Scatterplot
subplot(2,3,2)
scatter(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_Daily_forMAOD,yfit_Central_Daily,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 30])
set(gca,'xtick',[],'xticklabel',[], 'ytick',[], 'FontSize', 15)

title('Region B')

txt = ['y = ' num2str(round(p_Central_Daily(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_MAOD_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_MAOD_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_Daily_forMAOD),3,'significant'))];
Region = ('MAOD');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10, 0.36, p_txt,'FontSize', 15)

text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region

% Central Daily Scatterplot
subplot(2,3,5)
scatter(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_Daily_forAODc,yfit_Central_Daily_AODc_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 30])
set(gca, 'ytick',[], 'FontSize', 15)



txt = ['y = ' num2str(round(p_Central_Daily(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_AODc_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_AODc_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_Daily_forAODc),3,'significant'))];
Region = ('AOD_C');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10, 0.36, p_txt,'FontSize', 15)

text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% Southernmost Daily Scatterplot:
subplot(2,3,3)
scatter(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Southernmost_Wind_Daily_forMAOD,yfit_Southernmost_Daily_MAOD_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 30])
set(gca,'xtick',[],'xticklabel',[],'ytick',[], 'FontSize', 15)
 
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
title('Region C')
txt = ['y = ' num2str(round(p_Southernmost_Daily_MAOD_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_MAOD_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_MAOD_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_MAOD_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_Daily_forMAOD),3,'significant'))];
Region = ('MAOD');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10,0.36,p_txt, 'FontSize', 15)
text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
clear txt R_txt N_txt Region



% Southernmost Daily Scatterplot:
subplot(2,3,6)
scatter(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;

hold on
plot(x_Southernmost_Wind_Daily_forAODc,yfit_Southernmost_Daily_AODc_Wind,'-',...
    'LineWidth', 2)
ylim([0 0.5])
xlim([0 30])
set(gca,'ytick',[],'FontSize', 15)
 
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));

txt = ['y = ' num2str(round(p_Southernmost_Daily_AODc_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_AODc_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_AODc_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_AODc_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_Daily_forAODc),3,'significant'))];
Region = ('AOD_C');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10,0.36,p_txt, 'FontSize', 15)
text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
clear txt R_txt N_txt Region



%%

set(gcf,'PaperPositionMode','auto')
print(gcf,'DIURNALWINDS_OPENOCEAN_MAOD_CoarseAOD_allregions_daily_scatterplot_v2.png','-dpng','-r300');       %  *// 300 dpi


%%
% 
% % Total AOD & MAOD vs Daily Diurnal Winds:
% 
% %%
% clear x y nanVals p yfit R y_Corner_AODc_Daily y_Central_AODc_Daily y_Southernmost_AODc_Daily
% 
% x_Corner_Wind_Daily_forAODc = Corner_complete_TT_daily.MASTER_Winds_daytime;
% y_Corner_AODc_Daily = Corner_complete_TT_daily.Corner_AOD_total_openocean;
% 
% nanVals_Corner_AODc_Wind = ismissing(x_Corner_Wind_Daily_forAODc) | ismissing(y_Corner_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Corner_Wind_Daily_forAODc = x_Corner_Wind_Daily_forAODc(~nanVals_Corner_AODc_Wind);
% y_Corner_AODc_Daily = y_Corner_AODc_Daily(~nanVals_Corner_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Corner_Daily_AODc_Wind = polyfit(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Corner_Daily_AODc_Wind = polyval(p_Corner_Daily_AODc_Wind,x_Corner_Wind_Daily_forAODc);
% 
% [R_Corner_Daily_AODc_Wind, pval_Corner_Daily_AODc_Wind] = corrcoef(x_Corner_Wind_Daily_forAODc, y_Corner_AODc_Daily );
% 
% 
% % Central daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% x_Central_Wind_Daily_forAODc = Central_complete_TT_daily.MASTER_Winds_daytime;
% y_Central_AODc_Daily = Central_complete_TT_daily.Central_AOD_total_openocean;
% 
% nanVals_Central_AODc_Wind = ismissing(x_Central_Wind_Daily_forAODc) | ismissing(y_Central_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Central_Wind_Daily_forAODc = x_Central_Wind_Daily_forAODc(~nanVals_Central_AODc_Wind);
% y_Central_AODc_Daily = y_Central_AODc_Daily(~nanVals_Central_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Central_Daily_AODc_Wind = polyfit(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Central_Daily_AODc_Wind = polyval(p_Central_Daily_AODc_Wind,x_Central_Wind_Daily_forAODc);
% 
% [R_Central_Daily_AODc_Wind, pval_Central_Daily_AODc_Wind] = corrcoef(x_Central_Wind_Daily_forAODc, y_Central_AODc_Daily );
% 
% 
% % Southernmost daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% clear x y
% % load count.dat
% x_Southernmost_Wind_Daily_forAODc = Southernmost_complete_TT_daily.MASTER_Winds_daytime;
% y_Southernmost_AODc_Daily = Southernmost_complete_TT_daily.Southernmost_AOD_total_openocean;
% 
% nanVals_Southernmost_AODc_Wind = ismissing(x_Southernmost_Wind_Daily_forAODc) | ismissing(y_Southernmost_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Southernmost_Wind_Daily_forAODc = x_Southernmost_Wind_Daily_forAODc(~nanVals_Southernmost_AODc_Wind);
% y_Southernmost_AODc_Daily = y_Southernmost_AODc_Daily(~nanVals_Southernmost_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Southernmost_Daily_AODc_Wind = polyfit(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Southernmost_Daily_AODc_Wind = polyval(p_Southernmost_Daily_AODc_Wind,x_Southernmost_Wind_Daily_forAODc);
% 
% [R_Southernmost_Daily_AODc_Wind, pval_Southernmost_Daily_AODc_Wind] = corrcoef(x_Southernmost_Wind_Daily_forAODc, y_Southernmost_AODc_Daily );
% 
% %%
% 
% 
% clear x y nanVals p yfit R
% x_Corner_Wind_Daily_forMAOD = Corner_complete_TT_daily.MASTER_Winds_nighttime;
% y_Corner_MAOD_Daily = Corner_complete_TT_daily.MAOD_nocleanair_openocean;
% 
% nanVals_Corner_MAOD_Wind = ismissing(x_Corner_Wind_Daily_forMAOD) | ismissing(y_Corner_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Corner_Wind_Daily_forMAOD = x_Corner_Wind_Daily_forMAOD(~nanVals_Corner_MAOD_Wind);
% y_Corner_MAOD_Daily = y_Corner_MAOD_Daily(~nanVals_Corner_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Corner_Daily_MAOD_Wind = polyfit(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Corner_Daily_MAOD_Wind = polyval(p_Corner_Daily_MAOD_Wind,x_Corner_Wind_Daily_forMAOD);
% 
% [R_Corner_Daily_MAOD_Wind, pval_Corner_Daily_MAOD_Wind] = corrcoef(x_Corner_Wind_Daily_forMAOD, y_Corner_MAOD_Daily );
% 
% 
% % Central daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% x_Central_Wind_Daily_forMAOD = Central_complete_TT_daily.MASTER_Winds_nighttime;
% y_Central_MAOD_Daily = Central_complete_TT_daily.MAOD_nocleanair_openocean;
% 
% nanVals_Central_MAOD_Wind = ismissing(x_Central_Wind_Daily_forMAOD) | ismissing(y_Central_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Central_Wind_Daily_forMAOD = x_Central_Wind_Daily_forMAOD(~nanVals_Central_MAOD_Wind);
% y_Central_MAOD_Daily = y_Central_MAOD_Daily(~nanVals_Central_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Central_Daily = polyfit(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Central_Daily = polyval(p_Central_Daily,x_Central_Wind_Daily_forMAOD);
% 
% [R_Central_Daily_MAOD_Wind, pval_Central_Daily_MAOD_Wind] = corrcoef(x_Central_Wind_Daily_forMAOD, y_Central_MAOD_Daily );
% 
% 
% % Southernmost daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% clear x y
% % load count.dat
% x_Southernmost_Wind_Daily_forMAOD = Southernmost_complete_TT_daily.MASTER_Winds_nighttime;
% y_Southernmost_MAOD_Daily = Southernmost_complete_TT_daily.MAOD_nocleanair_openocean;
% 
% nanVals_Southernmost_MAOD_Wind = ismissing(x_Southernmost_Wind_Daily_forMAOD) | ismissing(y_Southernmost_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Southernmost_Wind_Daily_forMAOD = x_Southernmost_Wind_Daily_forMAOD(~nanVals_Southernmost_MAOD_Wind);
% y_Southernmost_MAOD_Daily = y_Southernmost_MAOD_Daily(~nanVals_Southernmost_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Southernmost_Daily_MAOD_Wind = polyfit(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Southernmost_Daily_MAOD_Wind = polyval(p_Southernmost_Daily_MAOD_Wind,x_Southernmost_Wind_Daily_forMAOD);
% 
% [R_Southernmost_Daily_MAOD_Wind, pval_Southernmost_Daily_MAOD_Wind] = corrcoef(x_Southernmost_Wind_Daily_forMAOD, y_Southernmost_MAOD_Daily );
% 
% 
% %%
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
% if ~make_it_tight,  clear subplot;  end
% 
% fig = figure; clf;
% 
% %Corner Daily Scatterplot
% subplot(2,3,1)
% scatter(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Corner_Wind_Daily_forMAOD,yfit_Corner_Daily_MAOD_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.5])
% xlim([0 30])
% set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)
% 
% title('Region A') 
% 
% txt = ['y = ' num2str(p_Corner_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Corner_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Corner_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Corner_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Corner_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10, 0.36, p_txt,'FontSize', 15)
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% % Corner Daily Scatterplot:
% subplot(2,3,4)
% scatter(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Corner_Wind_Daily_forAODc,yfit_Corner_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.5])
% xlim([0 30])
% set(gca,'FontSize', 15)
% 
% 
% txt = ['y = ' num2str(p_Corner_Daily_AODc_Wind(1)) 'x + ' num2str(p_Corner_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Corner_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Corner_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Corner_Wind_Daily_forAODc))];
% Region = ('AOD_T');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10, 0.36, p_txt,'FontSize', 15)
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% % Central Daily Scatterplot
% subplot(2,3,2)
% scatter(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Central_Wind_Daily_forMAOD,yfit_Central_Daily,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.5])
% xlim([0 30])
% set(gca,'xtick',[],'xticklabel',[], 'ytick',[], 'FontSize', 15)
% 
% title('Region B')
% 
% txt = ['y = ' num2str(p_Central_Daily(1)) 'x + ' num2str(p_Central_Daily(2))];
% R_txt = ['R = ' num2str(R_Central_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Central_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Central_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10, 0.36, p_txt,'FontSize', 15)
% 
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% % Central Daily Scatterplot
% subplot(2,3,5)
% scatter(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Central_Wind_Daily_forAODc,yfit_Central_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.5])
% xlim([0 30])
% set(gca, 'ytick',[], 'FontSize', 15)
% 
% 
% 
% txt = ['y = ' num2str(p_Central_Daily(1)) 'x + ' num2str(p_Central_Daily(2))];
% R_txt = ['R = ' num2str(R_Central_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Central_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Central_Wind_Daily_forAODc))];
% Region = ('AOD_T');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10, 0.36, p_txt,'FontSize', 15)
% 
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% % Southernmost Daily Scatterplot:
% subplot(2,3,3)
% scatter(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Southernmost_Wind_Daily_forMAOD,yfit_Southernmost_Daily_MAOD_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.5])
% xlim([0 30])
% set(gca,'xtick',[],'xticklabel',[],'ytick',[], 'FontSize', 15)
%  
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% title('Region C')
% txt = ['y = ' num2str(p_Southernmost_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Southernmost_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Southernmost_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Southernmost_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Southernmost_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10,0.36,p_txt, 'FontSize', 15)
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% clear txt R_txt N_txt Region
% 
% 
% 
% % Southernmost Daily Scatterplot:
% subplot(2,3,6)
% scatter(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% 
% hold on
% plot(x_Southernmost_Wind_Daily_forAODc,yfit_Southernmost_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% ylim([0 0.5])
% xlim([0 30])
% set(gca,'ytick',[],'FontSize', 15)
%  
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% 
% txt = ['y = ' num2str(p_Southernmost_Daily_AODc_Wind(1)) 'x + ' num2str(p_Southernmost_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Southernmost_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Southernmost_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Southernmost_Wind_Daily_forAODc))];
% Region = ('AOD_T');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10,0.36,p_txt, 'FontSize', 15)
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% clear txt R_txt N_txt Region
% 
% 
% 
% %%
% 
% set(gcf,'PaperPositionMode','auto')
% print(gcf,'DIURNALWINDS_OPENOCEAN_MAOD_TotalAOD_allregions_daily_scatterplot.png','-dpng','-r300');       %  *// 300 dpi
% 
% 
% 
% %%
% 
% % Total AOD & MAOD vs Daily Diurnal Winds:
% 
% % %%
% % clear x y nanVals p yfit R y_Corner_AODc_Daily y_Central_AODc_Daily y_Southernmost_AODc_Daily
% % 
% % x_Corner_Wind_Daily_forAODc = Corner_complete_TT_daily.MASTER_Winds_daytime;
% % y_Corner_AODc_Daily = Corner_complete_TT_daily.Corner_AOD_fine_openocean;
% % 
% % nanVals_Corner_AODc_Wind = ismissing(x_Corner_Wind_Daily_forAODc) | ismissing(y_Corner_AODc_Daily); % indices of values that are NaN in x or y
% % % resampling x and y to exclude any entries with NaN Values
% % x_Corner_Wind_Daily_forAODc = x_Corner_Wind_Daily_forAODc(~nanVals_Corner_AODc_Wind);
% % y_Corner_AODc_Daily = y_Corner_AODc_Daily(~nanVals_Corner_AODc_Wind);
% % 
% % % Use polyfit to compute a linear regression that predicts y from x:
% % p_Corner_Daily_AODc_Wind = polyfit(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,1);
% % 
% % % Call polyval to use p to predict y, calling the result yfit:
% % yfit_Corner_Daily_AODc_Wind = polyval(p_Corner_Daily_AODc_Wind,x_Corner_Wind_Daily_forAODc);
% % 
% % [R_Corner_Daily_AODc_Wind, pval_Corner_Daily_AODc_Wind] = corrcoef(x_Corner_Wind_Daily_forAODc, y_Corner_AODc_Daily );
% % 
% % 
% % % Central daily scatterplot
% % 
% % 
% % % Trying to plot MAOD against windspeed, let's start with corner region.
% % % I've already plotted out the histogram. :
% % 
% % x_Central_Wind_Daily_forAODc = Central_complete_TT_daily.MASTER_Winds_daytime;
% % y_Central_AODc_Daily = Central_complete_TT_daily.Central_AOD_fine_openocean;
% % 
% % nanVals_Central_AODc_Wind = ismissing(x_Central_Wind_Daily_forAODc) | ismissing(y_Central_AODc_Daily); % indices of values that are NaN in x or y
% % % resampling x and y to exclude any entries with NaN Values
% % x_Central_Wind_Daily_forAODc = x_Central_Wind_Daily_forAODc(~nanVals_Central_AODc_Wind);
% % y_Central_AODc_Daily = y_Central_AODc_Daily(~nanVals_Central_AODc_Wind);
% % 
% % % Use polyfit to compute a linear regression that predicts y from x:
% % p_Central_Daily_AODc_Wind = polyfit(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,1);
% % 
% % % Call polyval to use p to predict y, calling the result yfit:
% % yfit_Central_Daily_AODc_Wind = polyval(p_Central_Daily_AODc_Wind,x_Central_Wind_Daily_forAODc);
% % 
% % [R_Central_Daily_AODc_Wind, pval_Central_Daily_AODc_Wind] = corrcoef(x_Central_Wind_Daily_forAODc, y_Central_AODc_Daily );
% % 
% % 
% % % Southernmost daily scatterplot
% % 
% % 
% % % Trying to plot MAOD against windspeed, let's start with corner region.
% % % I've already plotted out the histogram. :
% % 
% % clear x y
% % % load count.dat
% % x_Southernmost_Wind_Daily_forAODc = Southernmost_complete_TT_daily.MASTER_Winds_daytime;
% % y_Southernmost_AODc_Daily = Southernmost_complete_TT_daily.Southernmost_AOD_fine_openocean;
% % 
% % nanVals_Southernmost_AODc_Wind = ismissing(x_Southernmost_Wind_Daily_forAODc) | ismissing(y_Southernmost_AODc_Daily); % indices of values that are NaN in x or y
% % % resampling x and y to exclude any entries with NaN Values
% % x_Southernmost_Wind_Daily_forAODc = x_Southernmost_Wind_Daily_forAODc(~nanVals_Southernmost_AODc_Wind);
% % y_Southernmost_AODc_Daily = y_Southernmost_AODc_Daily(~nanVals_Southernmost_AODc_Wind);
% % 
% % % Use polyfit to compute a linear regression that predicts y from x:
% % p_Southernmost_Daily_AODc_Wind = polyfit(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily,1);
% % 
% % % Call polyval to use p to predict y, calling the result yfit:
% % yfit_Southernmost_Daily_AODc_Wind = polyval(p_Southernmost_Daily_AODc_Wind,x_Southernmost_Wind_Daily_forAODc);
% % 
% % [R_Southernmost_Daily_AODc_Wind, pval_Southernmost_Daily_AODc_Wind] = corrcoef(x_Southernmost_Wind_Daily_forAODc, y_Southernmost_AODc_Daily );
% % 
% % %%
% % 
% % 
% clear x y nanVals p yfit R
% x_Corner_Wind_Daily_forMAOD = Corner_complete_TT_daily.MASTER_Winds_nighttime;
% y_Corner_MAOD_Daily = Corner_complete_TT_daily.MAOD_nocleanair_openocean;
% 
% nanVals_Corner_MAOD_Wind = ismissing(x_Corner_Wind_Daily_forMAOD) | ismissing(y_Corner_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Corner_Wind_Daily_forMAOD = x_Corner_Wind_Daily_forMAOD(~nanVals_Corner_MAOD_Wind);
% y_Corner_MAOD_Daily = y_Corner_MAOD_Daily(~nanVals_Corner_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Corner_Daily_MAOD_Wind = polyfit(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Corner_Daily_MAOD_Wind = polyval(p_Corner_Daily_MAOD_Wind,x_Corner_Wind_Daily_forMAOD);
% 
% [R_Corner_Daily_MAOD_Wind, pval_Corner_Daily_MAOD_Wind] = corrcoef(x_Corner_Wind_Daily_forMAOD, y_Corner_MAOD_Daily );
% 
% 
% % Central daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% x_Central_Wind_Daily_forMAOD = Central_complete_TT_daily.MASTER_Winds_nighttime;
% y_Central_MAOD_Daily = Central_complete_TT_daily.MAOD_nocleanair_openocean;
% 
% nanVals_Central_MAOD_Wind = ismissing(x_Central_Wind_Daily_forMAOD) | ismissing(y_Central_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Central_Wind_Daily_forMAOD = x_Central_Wind_Daily_forMAOD(~nanVals_Central_MAOD_Wind);
% y_Central_MAOD_Daily = y_Central_MAOD_Daily(~nanVals_Central_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Central_Daily = polyfit(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Central_Daily = polyval(p_Central_Daily,x_Central_Wind_Daily_forMAOD);
% 
% [R_Central_Daily_MAOD_Wind, pval_Central_Daily_MAOD_Wind] = corrcoef(x_Central_Wind_Daily_forMAOD, y_Central_MAOD_Daily );
% 
% 
% % Southernmost daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% clear x y
% % load count.dat
% x_Southernmost_Wind_Daily_forMAOD = Southernmost_complete_TT_daily.MASTER_Winds_nighttime;
% y_Southernmost_MAOD_Daily = Southernmost_complete_TT_daily.MAOD_nocleanair_openocean;
% 
% nanVals_Southernmost_MAOD_Wind = ismissing(x_Southernmost_Wind_Daily_forMAOD) | ismissing(y_Southernmost_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Southernmost_Wind_Daily_forMAOD = x_Southernmost_Wind_Daily_forMAOD(~nanVals_Southernmost_MAOD_Wind);
% y_Southernmost_MAOD_Daily = y_Southernmost_MAOD_Daily(~nanVals_Southernmost_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Southernmost_Daily_MAOD_Wind = polyfit(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Southernmost_Daily_MAOD_Wind = polyval(p_Southernmost_Daily_MAOD_Wind,x_Southernmost_Wind_Daily_forMAOD);
% 
% [R_Southernmost_Daily_MAOD_Wind, pval_Southernmost_Daily_MAOD_Wind] = corrcoef(x_Southernmost_Wind_Daily_forMAOD, y_Southernmost_MAOD_Daily );
% 
% 
% %%
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
% if ~make_it_tight,  clear subplot;  end
% 
% fig = figure; clf;
% 
% %Corner Daily Scatterplot
% subplot(2,3,1)
% scatter(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Corner_Wind_Daily_forMAOD,yfit_Corner_Daily_MAOD_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.5])
% xlim([0 30])
% set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)
% 
% title('Region A') 
% 
% txt = ['y = ' num2str(p_Corner_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Corner_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Corner_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Corner_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Corner_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10, 0.36, p_txt,'FontSize', 15)
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% % Corner Daily Scatterplot:
% subplot(2,3,4)
% scatter(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Corner_Wind_Daily_forAODc,yfit_Corner_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.5])
% xlim([0 30])
% set(gca,'FontSize', 15)
% 
% 
% txt = ['y = ' num2str(p_Corner_Daily_AODc_Wind(1)) 'x + ' num2str(p_Corner_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Corner_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Corner_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Corner_Wind_Daily_forAODc))];
% Region = ('AOD_f');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10, 0.36, p_txt,'FontSize', 15)
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% % Central Daily Scatterplot
% subplot(2,3,2)
% scatter(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Central_Wind_Daily_forMAOD,yfit_Central_Daily,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.5])
% xlim([0 30])
% set(gca,'xtick',[],'xticklabel',[], 'ytick',[], 'FontSize', 15)
% 
% title('Region B')
% 
% txt = ['y = ' num2str(p_Central_Daily(1)) 'x + ' num2str(p_Central_Daily(2))];
% R_txt = ['R = ' num2str(R_Central_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Central_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Central_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10, 0.36, p_txt,'FontSize', 15)
% 
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% % Central Daily Scatterplot
% subplot(2,3,5)
% scatter(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Central_Wind_Daily_forAODc,yfit_Central_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.5])
% xlim([0 30])
% set(gca, 'ytick',[], 'FontSize', 15)
% 
% 
% 
% txt = ['y = ' num2str(p_Central_Daily(1)) 'x + ' num2str(p_Central_Daily(2))];
% R_txt = ['R = ' num2str(R_Central_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Central_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Central_Wind_Daily_forAODc))];
% Region = ('AOD_f');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10, 0.36, p_txt,'FontSize', 15)
% 
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% % Southernmost Daily Scatterplot:
% subplot(2,3,3)
% scatter(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Southernmost_Wind_Daily_forMAOD,yfit_Southernmost_Daily_MAOD_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.5])
% xlim([0 30])
% set(gca,'xtick',[],'xticklabel',[],'ytick',[], 'FontSize', 15)
%  
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% title('Region C')
% txt = ['y = ' num2str(p_Southernmost_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Southernmost_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Southernmost_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Southernmost_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Southernmost_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10,0.36,p_txt, 'FontSize', 15)
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% clear txt R_txt N_txt Region
% 
% 
% 
% % Southernmost Daily Scatterplot:
% subplot(2,3,6)
% scatter(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% 
% hold on
% plot(x_Southernmost_Wind_Daily_forAODc,yfit_Southernmost_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% ylim([0 0.5])
% xlim([0 30])
% set(gca,'ytick',[],'FontSize', 15)
%  
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% 
% txt = ['y = ' num2str(p_Southernmost_Daily_AODc_Wind(1)) 'x + ' num2str(p_Southernmost_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Southernmost_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Southernmost_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Southernmost_Wind_Daily_forAODc))];
% Region = ('AOD_f');
% 
% text(10, 0.45, txt,'FontSize', 15)
% text(10,0.42, R_txt,'FontSize', 15)
% text(10,0.39,N_txt,'FontSize', 15)
% text(10,0.36,p_txt, 'FontSize', 15)
% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% clear txt R_txt N_txt Region
% 
% 
% 
% %%
% 
% set(gcf,'PaperPositionMode','auto')
% print(gcf,'DIURNALWINDS_OPENOCEAN_MAOD_FineAOD_allregions_daily_scatterplot.png','-dpng','-r300');       %  *// 300 dpi
% 
% 
% 
% %% 
% % ANNUAL DATA SCATTERPLOTS
% clear;
% cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_vars
% load('Central_complete_TT_annual.mat')
% load('Southernmost_complete_TT_annual.mat')
% load('Corner_complete_TT_annual.mat')
% 
% %% Scatter plot, annual data:
% % I've kept the variable names the same so that plotting is easy using the
% % same script as above
% 
% clear x y nanVals p yfit R
% x_Corner_Wind_Daily_forAODc = Corner_complete_TT_annual.MASTER_Winds_daytime;
% y_Corner_AODc_Daily = Corner_complete_TT_annual.Corner_AODc_openocean;
% 
% nanVals_Corner_AODc_Wind = ismissing(x_Corner_Wind_Daily_forAODc) | ismissing(y_Corner_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Corner_Wind_Daily_forAODc = x_Corner_Wind_Daily_forAODc(~nanVals_Corner_AODc_Wind);
% y_Corner_AODc_Daily = y_Corner_AODc_Daily(~nanVals_Corner_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Corner_Daily_AODc_Wind = polyfit(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Corner_Daily_AODc_Wind = polyval(p_Corner_Daily_AODc_Wind,x_Corner_Wind_Daily_forAODc);
% 
% [R_Corner_Daily_AODc_Wind, pval_Corner_Daily_AODc_Wind] = corrcoef(x_Corner_Wind_Daily_forAODc, y_Corner_AODc_Daily );
% 
% 
% % Central daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% x_Central_Wind_Daily_forAODc = Central_complete_TT_annual.MASTER_Winds_daytime;
% y_Central_AODc_Daily = Central_complete_TT_annual.Central_AOD_coarse_daily;
% 
% nanVals_Central_AODc_Wind = ismissing(x_Central_Wind_Daily_forAODc) | ismissing(y_Central_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Central_Wind_Daily_forAODc = x_Central_Wind_Daily_forAODc(~nanVals_Central_AODc_Wind);
% y_Central_AODc_Daily = y_Central_AODc_Daily(~nanVals_Central_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Central_Daily_AODc_Wind = polyfit(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Central_Daily_AODc_Wind = polyval(p_Central_Daily_AODc_Wind,x_Central_Wind_Daily_forAODc);
% 
% [R_Central_Daily_AODc_Wind, pval_Central_Daily_AODc_Wind] = corrcoef(x_Central_Wind_Daily_forAODc, y_Central_AODc_Daily );
% 
% 
% % Southernmost daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% clear x y
% % load count.dat
% x_Southernmost_Wind_Daily_forAODc = Southernmost_complete_TT_annual.MASTER_Winds_daytime;
% y_Southernmost_AODc_Daily = Southernmost_complete_TT_annual.Southernmost_AOD_coarse_daily;
%     
% nanVals_Southernmost_AODc_Wind = ismissing(x_Southernmost_Wind_Daily_forAODc) | ismissing(y_Southernmost_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Southernmost_Wind_Daily_forAODc = x_Southernmost_Wind_Daily_forAODc(~nanVals_Southernmost_AODc_Wind);
% y_Southernmost_AODc_Daily = y_Southernmost_AODc_Daily(~nanVals_Southernmost_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Southernmost_Daily_AODc_Wind = polyfit(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Southernmost_Daily_AODc_Wind = polyval(p_Southernmost_Daily_AODc_Wind,x_Southernmost_Wind_Daily_forAODc);
% 
% [R_Southernmost_Daily_AODc_Wind, pval_Southernmost_Daily_AODc_Wind] = corrcoef(x_Southernmost_Wind_Daily_forAODc, y_Southernmost_AODc_Daily );
% 
% %%
% 
% x_Corner_Wind_Daily_forMAOD = Corner_complete_TT_annual.MASTER_Winds_nighttime;
% y_Corner_MAOD_Daily = Corner_complete_TT_annual.MAOD_nocleanair_openocean;
% 
% nanVals_Corner_MAOD_Wind = ismissing(x_Corner_Wind_Daily_forMAOD) | ismissing(y_Corner_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Corner_Wind_Daily_forMAOD = x_Corner_Wind_Daily_forMAOD(~nanVals_Corner_MAOD_Wind);
% y_Corner_MAOD_Daily = y_Corner_MAOD_Daily(~nanVals_Corner_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Corner_Daily_MAOD_Wind = polyfit(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Corner_Daily_MAOD_Wind = polyval(p_Corner_Daily_MAOD_Wind,x_Corner_Wind_Daily_forMAOD);
% 
% [R_Corner_Daily_MAOD_Wind, pval_Corner_Daily_MAOD_Wind] = corrcoef(x_Corner_Wind_Daily_forMAOD, y_Corner_MAOD_Daily );
% 
% 
% % Central yearly scatterplot
% 
% 
% 
% x_Central_Wind_Daily_forMAOD = Central_complete_TT_annual.MASTER_Winds_nighttime;
% y_Central_MAOD_Daily = Central_complete_TT_annual.MAOD_nocleanair_Central_TT_MAOD_Wind_annual_nighttime;
% 
% nanVals_Central_MAOD_Wind = ismissing(x_Central_Wind_Daily_forMAOD) | ismissing(y_Central_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Central_Wind_Daily_forMAOD = x_Central_Wind_Daily_forMAOD(~nanVals_Central_MAOD_Wind);
% y_Central_MAOD_Daily = y_Central_MAOD_Daily(~nanVals_Central_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Central_Daily_MAOD_Wind = polyfit(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Central_Daily = polyval(p_Central_Daily_MAOD_Wind,x_Central_Wind_Daily_forMAOD);
% 
% [R_Central_Daily_MAOD_Wind, pval_Central_Daily_MAOD_Wind] = corrcoef(x_Central_Wind_Daily_forMAOD, y_Central_MAOD_Daily );
% 
% 
% % Southernmost yearly scatterplot
% 
% clear x y
% % load count.dat
% x_Southernmost_Wind_Daily_forMAOD = Southernmost_complete_TT_annual.MASTER_Winds_nighttime;
% y_Southernmost_MAOD_Daily = Southernmost_complete_TT_annual.MAOD_nocleanair_Southernmost_TT_MAOD_Wind_annual_nighttime;
% 
% nanVals_Southernmost_MAOD_Wind = ismissing(x_Southernmost_Wind_Daily_forMAOD) | ismissing(y_Southernmost_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Southernmost_Wind_Daily_forMAOD = x_Southernmost_Wind_Daily_forMAOD(~nanVals_Southernmost_MAOD_Wind);
% y_Southernmost_MAOD_Daily = y_Southernmost_MAOD_Daily(~nanVals_Southernmost_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Southernmost_Daily_MAOD_Wind = polyfit(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Southernmost_Daily_MAOD_Wind = polyval(p_Southernmost_Daily_MAOD_Wind,x_Southernmost_Wind_Daily_forMAOD);
% 
% [R_Southernmost_Daily_MAOD_Wind, pval_Southernmost_Daily_MAOD_Wind] = corrcoef(x_Southernmost_Wind_Daily_forMAOD, y_Southernmost_MAOD_Daily );
% 
% 
% %%
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
% if ~make_it_tight,  clear subplot;  end
% 
% fig = figure; clf;
% 
% %Corner Daily Scatterplot
% subplot(2,3,1)
% scatter(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Corner_Wind_Daily_forMAOD,yfit_Corner_Daily_MAOD_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)
% 
% title('(a) Region A') 
% 
% txt = ['y = ' num2str(p_Corner_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Corner_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Corner_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Corner_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Corner_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% text(7, 0.1, txt,'FontSize', 15)
% text(7, 0.094, R_txt,'FontSize', 15)
% text(7, 0.088,N_txt,'FontSize', 15)
% text(7, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% % Corner Daily Scatterplot:
% subplot(2,3,4)
% scatter(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Corner_Wind_Daily_forAODc,yfit_Corner_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca,'FontSize', 15)
% 
% 
% txt = ['y = ' num2str(p_Corner_Daily_AODc_Wind(1)) 'x + ' num2str(p_Corner_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Corner_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Corner_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Corner_Wind_Daily_forAODc))];
% Region = ('AOD_C');
% 
% 
% 
% text(7, 0.1, txt,'FontSize', 15)
% text(7, 0.094, R_txt,'FontSize', 15)
% text(7, 0.088,N_txt,'FontSize', 15)
% text(7, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% % Central Daily Scatterplot
% subplot(2,3,2)
% scatter(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Central_Wind_Daily_forMAOD,yfit_Central_Daily,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca,'xtick',[],'xticklabel',[], 'ytick',[], 'FontSize', 15)
% 
% title('(b) Region B')
% 
% txt = ['y = ' num2str(p_Central_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Central_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Central_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Central_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Central_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% 
% 
% text(8.2, 0.1, txt,'FontSize', 15)
% text(8.2, 0.094, R_txt,'FontSize', 15)
% text(8.2, 0.088,N_txt,'FontSize', 15)
% text(8.2, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% 
% clear txt R_txt N_txt Region
% 
% % Central Daily Scatterplot
% subplot(2,3,5)
% scatter(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Central_Wind_Daily_forAODc,yfit_Central_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca, 'ytick',[], 'FontSize', 15)
% 
% 
% 
% txt = ['y = ' num2str(p_Central_Daily_AODc_Wind(1)) 'x + ' num2str(p_Central_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Central_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Central_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Central_Wind_Daily_forAODc))];
% Region = ('AOD_C');
% 
% 
% text(8.2, 0.1, txt,'FontSize', 15)
% text(8.2, 0.094, R_txt,'FontSize', 15)
% text(8.2, 0.088,N_txt,'FontSize', 15)
% text(8.2, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% % Southernmost Daily Scatterplot:
% subplot(2,3,3)
% scatter(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily, 20, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Southernmost_Wind_Daily_forMAOD,yfit_Southernmost_Daily_MAOD_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca,'xtick',[],'xticklabel',[],'ytick',[], 'FontSize', 15)
%  
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% title('(c) Region C')
% txt = ['y = ' num2str(p_Southernmost_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Southernmost_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Southernmost_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Southernmost_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Southernmost_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% 
% text(8, 0.1, txt,'FontSize', 15)
% text(8, 0.094, R_txt,'FontSize', 15)
% text(8, 0.088,N_txt,'FontSize', 15)
% text(8, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% % Southernmost Daily Scatterplot:
% subplot(2,3,6)
% scatter(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily, 20, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% 
% hold on
% plot(x_Southernmost_Wind_Daily_forAODc,yfit_Southernmost_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca,'ytick',[],'FontSize', 15)
%  
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% 
% txt = ['y = ' num2str(p_Southernmost_Daily_AODc_Wind(1)) 'x + ' num2str(p_Southernmost_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Southernmost_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Southernmost_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Southernmost_Wind_Daily_forAODc))];
% Region = ('AOD_C');
% 
% 
% 
% text(8, 0.1, txt,'FontSize', 15)
% text(8, 0.094, R_txt,'FontSize', 15)
% text(8, 0.088,N_txt,'FontSize', 15)
% text(8, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% %%
% cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures/Diurnal_winds
% set(gcf,'PaperPositionMode','auto')
% print(gcf,'DIURNALWINDS_OPENOCEAN_MAOD_Coarse_AOD_allregions_annual_scatterplot.png','-dpng','-r300');       %  *// 300 dpi
% 
% 
% 
% 
% %%
% %% Scatter plot, annual data:
% % I've kept the variable names the same so that plotting is easy using the
% % same script as above
% 
% clear x y nanVals p yfit R
% x_Corner_Wind_Daily_forAODc = Corner_complete_TT_annual.MASTER_Winds_daytime;
% y_Corner_AODc_Daily = Corner_complete_TT_annual.Corner_AOD_total_openocean;
% 
% nanVals_Corner_AODc_Wind = ismissing(x_Corner_Wind_Daily_forAODc) | ismissing(y_Corner_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Corner_Wind_Daily_forAODc = x_Corner_Wind_Daily_forAODc(~nanVals_Corner_AODc_Wind);
% y_Corner_AODc_Daily = y_Corner_AODc_Daily(~nanVals_Corner_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Corner_Daily_AODc_Wind = polyfit(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Corner_Daily_AODc_Wind = polyval(p_Corner_Daily_AODc_Wind,x_Corner_Wind_Daily_forAODc);
% 
% [R_Corner_Daily_AODc_Wind, pval_Corner_Daily_AODc_Wind] = corrcoef(x_Corner_Wind_Daily_forAODc, y_Corner_AODc_Daily );
% 
% 
% % Central daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% x_Central_Wind_Daily_forAODc = Central_complete_TT_annual.MASTER_Winds_daytime;
% y_Central_AODc_Daily = Central_complete_TT_annual.Central_AOD_total_daily;
% 
% nanVals_Central_AODc_Wind = ismissing(x_Central_Wind_Daily_forAODc) | ismissing(y_Central_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Central_Wind_Daily_forAODc = x_Central_Wind_Daily_forAODc(~nanVals_Central_AODc_Wind);
% y_Central_AODc_Daily = y_Central_AODc_Daily(~nanVals_Central_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Central_Daily_AODc_Wind = polyfit(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Central_Daily_AODc_Wind = polyval(p_Central_Daily_AODc_Wind,x_Central_Wind_Daily_forAODc);
% 
% [R_Central_Daily_AODc_Wind, pval_Central_Daily_AODc_Wind] = corrcoef(x_Central_Wind_Daily_forAODc, y_Central_AODc_Daily );
% 
% 
% % Southernmost daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% clear x y
% % load count.dat
% x_Southernmost_Wind_Daily_forAODc = Southernmost_complete_TT_annual.MASTER_Winds_daytime;
% y_Southernmost_AODc_Daily = Southernmost_complete_TT_annual.Southernmost_AOD_total_daily;
%     
% nanVals_Southernmost_AODc_Wind = ismissing(x_Southernmost_Wind_Daily_forAODc) | ismissing(y_Southernmost_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Southernmost_Wind_Daily_forAODc = x_Southernmost_Wind_Daily_forAODc(~nanVals_Southernmost_AODc_Wind);
% y_Southernmost_AODc_Daily = y_Southernmost_AODc_Daily(~nanVals_Southernmost_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Southernmost_Daily_AODc_Wind = polyfit(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Southernmost_Daily_AODc_Wind = polyval(p_Southernmost_Daily_AODc_Wind,x_Southernmost_Wind_Daily_forAODc);
% 
% [R_Southernmost_Daily_AODc_Wind, pval_Southernmost_Daily_AODc_Wind] = corrcoef(x_Southernmost_Wind_Daily_forAODc, y_Southernmost_AODc_Daily );
% 
% %%
% 
% x_Corner_Wind_Daily_forMAOD = Corner_complete_TT_annual.MASTER_Winds_nighttime;
% y_Corner_MAOD_Daily = Corner_complete_TT_annual.MAOD_nocleanair_openocean;
% 
% nanVals_Corner_MAOD_Wind = ismissing(x_Corner_Wind_Daily_forMAOD) | ismissing(y_Corner_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Corner_Wind_Daily_forMAOD = x_Corner_Wind_Daily_forMAOD(~nanVals_Corner_MAOD_Wind);
% y_Corner_MAOD_Daily = y_Corner_MAOD_Daily(~nanVals_Corner_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Corner_Daily_MAOD_Wind = polyfit(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Corner_Daily_MAOD_Wind = polyval(p_Corner_Daily_MAOD_Wind,x_Corner_Wind_Daily_forMAOD);
% 
% [R_Corner_Daily_MAOD_Wind, pval_Corner_Daily_MAOD_Wind] = corrcoef(x_Corner_Wind_Daily_forMAOD, y_Corner_MAOD_Daily );
% 
% 
% % Central yearly scatterplot
% 
% 
% 
% x_Central_Wind_Daily_forMAOD = Central_complete_TT_annual.MASTER_Winds_nighttime;
% y_Central_MAOD_Daily = Central_complete_TT_annual.MAOD_nocleanair_Central_TT_MAOD_Wind_annual_nighttime;
% 
% nanVals_Central_MAOD_Wind = ismissing(x_Central_Wind_Daily_forMAOD) | ismissing(y_Central_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Central_Wind_Daily_forMAOD = x_Central_Wind_Daily_forMAOD(~nanVals_Central_MAOD_Wind);
% y_Central_MAOD_Daily = y_Central_MAOD_Daily(~nanVals_Central_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Central_Daily_MAOD_Wind = polyfit(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Central_Daily = polyval(p_Central_Daily_MAOD_Wind,x_Central_Wind_Daily_forMAOD);
% 
% [R_Central_Daily_MAOD_Wind, pval_Central_Daily_MAOD_Wind] = corrcoef(x_Central_Wind_Daily_forMAOD, y_Central_MAOD_Daily );
% 
% 
% % Southernmost yearly scatterplot
% 
% clear x y
% % load count.dat
% x_Southernmost_Wind_Daily_forMAOD = Southernmost_complete_TT_annual.MASTER_Winds_nighttime;
% y_Southernmost_MAOD_Daily = Southernmost_complete_TT_annual.MAOD_nocleanair_Southernmost_TT_MAOD_Wind_annual_nighttime;
% 
% nanVals_Southernmost_MAOD_Wind = ismissing(x_Southernmost_Wind_Daily_forMAOD) | ismissing(y_Southernmost_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Southernmost_Wind_Daily_forMAOD = x_Southernmost_Wind_Daily_forMAOD(~nanVals_Southernmost_MAOD_Wind);
% y_Southernmost_MAOD_Daily = y_Southernmost_MAOD_Daily(~nanVals_Southernmost_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Southernmost_Daily_MAOD_Wind = polyfit(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Southernmost_Daily_MAOD_Wind = polyval(p_Southernmost_Daily_MAOD_Wind,x_Southernmost_Wind_Daily_forMAOD);
% 
% [R_Southernmost_Daily_MAOD_Wind, pval_Southernmost_Daily_MAOD_Wind] = corrcoef(x_Southernmost_Wind_Daily_forMAOD, y_Southernmost_MAOD_Daily );
% 
% 
% %%
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
% if ~make_it_tight,  clear subplot;  end
% 
% fig = figure; clf;
% 
% %Corner Daily Scatterplot
% subplot(2,3,1)
% scatter(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Corner_Wind_Daily_forMAOD,yfit_Corner_Daily_MAOD_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)
% 
% title('(a) Region A') 
% 
% txt = ['y = ' num2str(p_Corner_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Corner_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Corner_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Corner_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Corner_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% text(7, 0.1, txt,'FontSize', 15)
% text(7, 0.094, R_txt,'FontSize', 15)
% text(7, 0.088,N_txt,'FontSize', 15)
% text(7, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% % Corner Daily Scatterplot:
% subplot(2,3,4)
% scatter(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Corner_Wind_Daily_forAODc,yfit_Corner_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.06 0.14])
% xlim([6.5 11])
% set(gca,'FontSize', 15)
% 
% 
% txt = ['y = ' num2str(p_Corner_Daily_AODc_Wind(1)) 'x + ' num2str(p_Corner_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Corner_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Corner_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Corner_Wind_Daily_forAODc))];
% Region = ('AOD_T');
% 
% 
% 
% text(7, 0.1, txt,'FontSize', 15)
% text(7, 0.094, R_txt,'FontSize', 15)
% text(7, 0.088,N_txt,'FontSize', 15)
% text(7, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.135, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% % Central Daily Scatterplot
% subplot(2,3,2)
% scatter(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Central_Wind_Daily_forMAOD,yfit_Central_Daily,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca,'xtick',[],'xticklabel',[], 'ytick',[], 'FontSize', 15)
% 
% title('(b) Region B')
% 
% txt = ['y = ' num2str(p_Central_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Central_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Central_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Central_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Central_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% 
% 
% text(8.2, 0.1, txt,'FontSize', 15)
% text(8.2, 0.094, R_txt,'FontSize', 15)
% text(8.2, 0.088,N_txt,'FontSize', 15)
% text(8.2, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% 
% clear txt R_txt N_txt Region
% 
% % Central Daily Scatterplot
% subplot(2,3,5)
% scatter(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Central_Wind_Daily_forAODc,yfit_Central_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.06 0.14])
% xlim([6.5 11])
% set(gca, 'ytick',[], 'FontSize', 15)
% 
% 
% 
% txt = ['y = ' num2str(p_Central_Daily_AODc_Wind(1)) 'x + ' num2str(p_Central_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Central_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Central_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Central_Wind_Daily_forAODc))];
% Region = ('AOD_T');
% 
% 
% text(8.2, 0.1, txt,'FontSize', 15)
% text(8.2, 0.094, R_txt,'FontSize', 15)
% text(8.2, 0.088,N_txt,'FontSize', 15)
% text(8.2, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.135, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% % Southernmost Daily Scatterplot:
% subplot(2,3,3)
% scatter(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily, 20, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Southernmost_Wind_Daily_forMAOD,yfit_Southernmost_Daily_MAOD_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca,'xtick',[],'xticklabel',[],'ytick',[], 'FontSize', 15)
%  
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% title('(c) Region C')
% txt = ['y = ' num2str(p_Southernmost_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Southernmost_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Southernmost_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Southernmost_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Southernmost_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% 
% text(8, 0.1, txt,'FontSize', 15)
% text(8, 0.094, R_txt,'FontSize', 15)
% text(8, 0.088,N_txt,'FontSize', 15)
% text(8, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% % Southernmost Daily Scatterplot:
% subplot(2,3,6)
% scatter(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily, 20, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% 
% hold on
% plot(x_Southernmost_Wind_Daily_forAODc,yfit_Southernmost_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.06 0.14])
% xlim([6.5 11])
% set(gca,'ytick',[],'FontSize', 15)
%  
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% 
% txt = ['y = ' num2str(p_Southernmost_Daily_AODc_Wind(1)) 'x + ' num2str(p_Southernmost_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Southernmost_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Southernmost_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Southernmost_Wind_Daily_forAODc))];
% Region = ('AOD_T');
% 
% 
% 
% text(8, 0.1, txt,'FontSize', 15)
% text(8, 0.094, R_txt,'FontSize', 15)
% text(8, 0.088,N_txt,'FontSize', 15)
% text(8, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.135, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% %%
% cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures/Diurnal_winds
% set(gcf,'PaperPositionMode','auto')
% print(gcf,'DIURNALWINDS_OPENOCEAN_MAOD_Total_AOD_allregions_annual_scatterplot.png','-dpng','-r300');       %  *// 300 dpi
% 
% 
% %% Scatter plot, annual data:
% % I've kept the variable names the same so that plotting is easy using the
% % same script as above
% 
% % Fine AOD:
% clear x y nanVals p yfit R
% x_Corner_Wind_Daily_forAODc = Corner_complete_TT_annual.MASTER_Winds_daytime;
% y_Corner_AODc_Daily = Corner_complete_TT_annual.Corner_AOD_fine_openocean;
% 
% nanVals_Corner_AODc_Wind = ismissing(x_Corner_Wind_Daily_forAODc) | ismissing(y_Corner_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Corner_Wind_Daily_forAODc = x_Corner_Wind_Daily_forAODc(~nanVals_Corner_AODc_Wind);
% y_Corner_AODc_Daily = y_Corner_AODc_Daily(~nanVals_Corner_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Corner_Daily_AODc_Wind = polyfit(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Corner_Daily_AODc_Wind = polyval(p_Corner_Daily_AODc_Wind,x_Corner_Wind_Daily_forAODc);
% 
% [R_Corner_Daily_AODc_Wind, pval_Corner_Daily_AODc_Wind] = corrcoef(x_Corner_Wind_Daily_forAODc, y_Corner_AODc_Daily );
% 
% 
% % Central daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% x_Central_Wind_Daily_forAODc = Central_complete_TT_annual.MASTER_Winds_daytime;
% y_Central_AODc_Daily = Central_complete_TT_annual.Central_AOD_fine_daily;
% 
% nanVals_Central_AODc_Wind = ismissing(x_Central_Wind_Daily_forAODc) | ismissing(y_Central_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Central_Wind_Daily_forAODc = x_Central_Wind_Daily_forAODc(~nanVals_Central_AODc_Wind);
% y_Central_AODc_Daily = y_Central_AODc_Daily(~nanVals_Central_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Central_Daily_AODc_Wind = polyfit(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Central_Daily_AODc_Wind = polyval(p_Central_Daily_AODc_Wind,x_Central_Wind_Daily_forAODc);
% 
% [R_Central_Daily_AODc_Wind, pval_Central_Daily_AODc_Wind] = corrcoef(x_Central_Wind_Daily_forAODc, y_Central_AODc_Daily );
% 
% 
% % Southernmost daily scatterplot
% 
% 
% % Trying to plot MAOD against windspeed, let's start with corner region.
% % I've already plotted out the histogram. :
% 
% clear x y
% % load count.dat
% x_Southernmost_Wind_Daily_forAODc = Southernmost_complete_TT_annual.MASTER_Winds_daytime;
% y_Southernmost_AODc_Daily = Southernmost_complete_TT_annual.Southernmost_AOD_fine_daily;
%     
% nanVals_Southernmost_AODc_Wind = ismissing(x_Southernmost_Wind_Daily_forAODc) | ismissing(y_Southernmost_AODc_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Southernmost_Wind_Daily_forAODc = x_Southernmost_Wind_Daily_forAODc(~nanVals_Southernmost_AODc_Wind);
% y_Southernmost_AODc_Daily = y_Southernmost_AODc_Daily(~nanVals_Southernmost_AODc_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Southernmost_Daily_AODc_Wind = polyfit(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Southernmost_Daily_AODc_Wind = polyval(p_Southernmost_Daily_AODc_Wind,x_Southernmost_Wind_Daily_forAODc);
% 
% [R_Southernmost_Daily_AODc_Wind, pval_Southernmost_Daily_AODc_Wind] = corrcoef(x_Southernmost_Wind_Daily_forAODc, y_Southernmost_AODc_Daily );
% 
% %%
% 
% x_Corner_Wind_Daily_forMAOD = Corner_complete_TT_annual.MASTER_Winds_nighttime;
% y_Corner_MAOD_Daily = Corner_complete_TT_annual.MAOD_nocleanair_openocean;
% 
% nanVals_Corner_MAOD_Wind = ismissing(x_Corner_Wind_Daily_forMAOD) | ismissing(y_Corner_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Corner_Wind_Daily_forMAOD = x_Corner_Wind_Daily_forMAOD(~nanVals_Corner_MAOD_Wind);
% y_Corner_MAOD_Daily = y_Corner_MAOD_Daily(~nanVals_Corner_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Corner_Daily_MAOD_Wind = polyfit(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Corner_Daily_MAOD_Wind = polyval(p_Corner_Daily_MAOD_Wind,x_Corner_Wind_Daily_forMAOD);
% 
% [R_Corner_Daily_MAOD_Wind, pval_Corner_Daily_MAOD_Wind] = corrcoef(x_Corner_Wind_Daily_forMAOD, y_Corner_MAOD_Daily );
% 
% 
% % Central yearly scatterplot
% 
% 
% 
% x_Central_Wind_Daily_forMAOD = Central_complete_TT_annual.MASTER_Winds_nighttime;
% y_Central_MAOD_Daily = Central_complete_TT_annual.MAOD_nocleanair_Central_TT_MAOD_Wind_annual_nighttime;
% 
% nanVals_Central_MAOD_Wind = ismissing(x_Central_Wind_Daily_forMAOD) | ismissing(y_Central_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Central_Wind_Daily_forMAOD = x_Central_Wind_Daily_forMAOD(~nanVals_Central_MAOD_Wind);
% y_Central_MAOD_Daily = y_Central_MAOD_Daily(~nanVals_Central_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Central_Daily_MAOD_Wind = polyfit(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Central_Daily = polyval(p_Central_Daily_MAOD_Wind,x_Central_Wind_Daily_forMAOD);
% 
% [R_Central_Daily_MAOD_Wind, pval_Central_Daily_MAOD_Wind] = corrcoef(x_Central_Wind_Daily_forMAOD, y_Central_MAOD_Daily );
% 
% 
% % Southernmost yearly scatterplot
% 
% clear x y
% % load count.dat
% x_Southernmost_Wind_Daily_forMAOD = Southernmost_complete_TT_annual.MASTER_Winds_nighttime;
% y_Southernmost_MAOD_Daily = Southernmost_complete_TT_annual.MAOD_nocleanair_Southernmost_TT_MAOD_Wind_annual_nighttime;
% 
% nanVals_Southernmost_MAOD_Wind = ismissing(x_Southernmost_Wind_Daily_forMAOD) | ismissing(y_Southernmost_MAOD_Daily); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Southernmost_Wind_Daily_forMAOD = x_Southernmost_Wind_Daily_forMAOD(~nanVals_Southernmost_MAOD_Wind);
% y_Southernmost_MAOD_Daily = y_Southernmost_MAOD_Daily(~nanVals_Southernmost_MAOD_Wind);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Southernmost_Daily_MAOD_Wind = polyfit(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Southernmost_Daily_MAOD_Wind = polyval(p_Southernmost_Daily_MAOD_Wind,x_Southernmost_Wind_Daily_forMAOD);
% 
% [R_Southernmost_Daily_MAOD_Wind, pval_Southernmost_Daily_MAOD_Wind] = corrcoef(x_Southernmost_Wind_Daily_forMAOD, y_Southernmost_MAOD_Daily );
% 
% 
% %%
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
% if ~make_it_tight,  clear subplot;  end
% 
% fig = figure; clf;
% 
% %Corner Daily Scatterplot
% subplot(2,3,1)
% scatter(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Corner_Wind_Daily_forMAOD,yfit_Corner_Daily_MAOD_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)
% 
% title('(a) Region A') 
% 
% txt = ['y = ' num2str(p_Corner_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Corner_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Corner_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Corner_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Corner_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% text(7, 0.1, txt,'FontSize', 15)
% text(7, 0.094, R_txt,'FontSize', 15)
% text(7, 0.088,N_txt,'FontSize', 15)
% text(7, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% % Corner Daily Scatterplot:
% subplot(2,3,4)
% scatter(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Corner_Wind_Daily_forAODc,yfit_Corner_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.12])
% xlim([6.5 11])
% set(gca,'FontSize', 15)
% 
% 
% txt = ['y = ' num2str(p_Corner_Daily_AODc_Wind(1)) 'x + ' num2str(p_Corner_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Corner_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Corner_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Corner_Wind_Daily_forAODc))];
% Region = ('AOD_f');
% 
% 
% 
% text(7, 0.1, txt,'FontSize', 15)
% text(7, 0.094, R_txt,'FontSize', 15)
% text(7, 0.088,N_txt,'FontSize', 15)
% text(7, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% % Central Daily Scatterplot
% subplot(2,3,2)
% scatter(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Central_Wind_Daily_forMAOD,yfit_Central_Daily,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([6.5 11])
% set(gca,'xtick',[],'xticklabel',[], 'ytick',[], 'FontSize', 15)
% 
% title('(b) Region B')
% 
% txt = ['y = ' num2str(p_Central_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Central_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Central_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Central_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Central_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% 
% 
% text(8.2, 0.1, txt,'FontSize', 15)
% text(8.2, 0.094, R_txt,'FontSize', 15)
% text(8.2, 0.088,N_txt,'FontSize', 15)
% text(8.2, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% 
% clear txt R_txt N_txt Region
% 
% % Central Daily Scatterplot
% subplot(2,3,5)
% scatter(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Central_Wind_Daily_forAODc,yfit_Central_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.12])
% xlim([6.5 11])
% set(gca, 'ytick',[], 'FontSize', 15)
% 
% 
% 
% txt = ['y = ' num2str(p_Central_Daily_AODc_Wind(1)) 'x + ' num2str(p_Central_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Central_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Central_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Central_Wind_Daily_forAODc))];
% Region = ('AOD_f');
% 
% 
% text(8.2, 0.1, txt,'FontSize', 15)
% text(8.2, 0.094, R_txt,'FontSize', 15)
% text(8.2, 0.088,N_txt,'FontSize', 15)
% text(8.2, 0.082, p_txt,'FontSize', 15)
% text(6.6, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% % Southernmost Daily Scatterplot:
% subplot(2,3,3)
% scatter(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily, 20, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% hold on
% plot(x_Southernmost_Wind_Daily_forMAOD,yfit_Southernmost_Daily_MAOD_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0.025 0.125])
% xlim([5 10])
% set(gca,'xtick',[],'xticklabel',[],'ytick',[], 'FontSize', 15)
%  
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% title('(c) Region C')
% txt = ['y = ' num2str(p_Southernmost_Daily_MAOD_Wind(1)) 'x + ' num2str(p_Southernmost_Daily_MAOD_Wind(2))];
% R_txt = ['R = ' num2str(R_Southernmost_Daily_MAOD_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Southernmost_Daily_MAOD_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Southernmost_Wind_Daily_forMAOD))];
% Region = ('MAOD');
% 
% 
% text(7, 0.1, txt,'FontSize', 15)
% text(7, 0.094, R_txt,'FontSize', 15)
% text(7, 0.088,N_txt,'FontSize', 15)
% text(7, 0.082, p_txt,'FontSize', 15)
% text(5.1, 0.12, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% % Southernmost Daily Scatterplot:
% subplot(2,3,6)
% scatter(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily, 20, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
% 
% hold on
% plot(x_Southernmost_Wind_Daily_forAODc,yfit_Southernmost_Daily_AODc_Wind,'-',...
%     'LineWidth', 2)
% 
% ylim([0 0.12])
% xlim([5 10])
% set(gca,'ytick',[],'FontSize', 15)
%  
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% 
% txt = ['y = ' num2str(p_Southernmost_Daily_AODc_Wind(1)) 'x + ' num2str(p_Southernmost_Daily_AODc_Wind(2))];
% R_txt = ['R = ' num2str(R_Southernmost_Daily_AODc_Wind(1,2))];
% p_txt = ['p = ' num2str(pval_Southernmost_Daily_AODc_Wind(1,2))];
% N_txt = ['n = ' num2str(length(x_Southernmost_Wind_Daily_forAODc))];
% Region = ('AOD_f');
% 
% 
% 
% text(7, 0.1, txt,'FontSize', 15)
% text(7, 0.094, R_txt,'FontSize', 15)
% text(7, 0.088,N_txt,'FontSize', 15)
% text(7, 0.082, p_txt,'FontSize', 15)
% text(5.1, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')
% 
% clear txt R_txt N_txt Region
% 
% 
% 
% %%
% cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures/Diurnal_winds
% set(gcf,'PaperPositionMode','auto')
% print(gcf,'DIURNALWINDS_OPENOCEAN_MAOD_Fine_AOD_allregions_annual_scatterplot.png','-dpng','-r300');       %  *// 300 dpi
% 
% 

%%


% OPEN OCEAN PIXELS ONLY 

%% SUBPLOTTED WINTER VS SUMMER CORNER AND CENTRAL REGIONS:

% I toggled through these, changing the variables from
% Corner_AODc_openocean to Corner_AOD_total_openocean, and so on.

% Corner Region
clear x y nanVals x_winter y_winter
x_winter_Corner = Spring_TT_Corner_BSea_daily.MASTER_Winds_daytime;
y_winter_Corner = Spring_TT_Corner_BSea_daily.Corner_AODc_openocean;
x_summer_Corner = Summer_TT_Corner_BSea_daily.MASTER_Winds_daytime;
y_summer_Corner = Summer_TT_Corner_BSea_daily.Corner_AODc_openocean;

% resampling x and y to exclude any entries with NaN Values
nanVals_winter_Corner = ismissing(x_winter_Corner) | ismissing(y_winter_Corner); 
nanVals_summer_Corner = ismissing(x_summer_Corner) | ismissing(y_summer_Corner);

x_winter_Corner = x_winter_Corner(~nanVals_winter_Corner);
y_winter_Corner = y_winter_Corner(~nanVals_winter_Corner);

x_summer_Corner = x_summer_Corner(~nanVals_summer_Corner);
y_summer_Corner = y_summer_Corner(~nanVals_summer_Corner);

% Use polyfit to compute a linear regression that predicts y from x:
p_winter_Corner = polyfit(x_winter_Corner,y_winter_Corner,1);
p_summer_Corner = polyfit(x_summer_Corner, y_summer_Corner,1);

% p(1) is the slope and p(2) is the intercept of the linear predictor. You can also obtain regression coefficients using the Basic Fitting UI.
% Call polyval to use p to predict y, calling the result yfit:
yfit_winter_Corner = polyval(p_winter_Corner,x_winter_Corner);
yfit_summer_Corner = polyval(p_summer_Corner, x_summer_Corner);

[R_winter_Corner, pval_winter_Corner] = corrcoef(x_winter_Corner, y_winter_Corner );
[R_summer_Corner, pval_summer_Corner] = corrcoef(x_summer_Corner, y_summer_Corner);

% 
% % Central Region SPRING vs SUMMER
x_spring_Central = Spring_TT_Central_BSea_daily.MASTER_Winds_daytime;
y_spring_Central = Spring_TT_Central_BSea_daily.Central_AODc_openocean;
x_summer_Central = Summer_TT_Central_BSea_daily.MASTER_Winds_daytime;
y_summer_Central = Summer_TT_Central_BSea_daily.Central_AODc_openocean;

% resampling x and y to exclude any entries with NaN Values
nanVals_winter_Central = ismissing(x_spring_Central) | ismissing(y_spring_Central); 
nanVals_summer_Central = ismissing(x_summer_Central) | ismissing(y_summer_Central);

x_spring_Central = x_spring_Central(~nanVals_winter_Central);
y_spring_Central = y_spring_Central(~nanVals_winter_Central);

x_summer_Central = x_summer_Central(~nanVals_summer_Central);
y_summer_Central = y_summer_Central(~nanVals_summer_Central);

% Use polyfit to compute a linear regression that predicts y from x:
p_spring_Central = polyfit(x_spring_Central,y_spring_Central,1);
p_summer_Central = polyfit(x_summer_Central, y_summer_Central,1);

% p(1) is the slope and p(2) is the intercept of the linear predictor. You can also obtain regression coefficients using the Basic Fitting UI.
% Call polyval to use p to predict y, calling the result yfit:
yfit_spring_Central = polyval(p_spring_Central,x_spring_Central);
yfit_summer_Central = polyval(p_summer_Central, x_summer_Central);

[R_spring_Central, pval_spring_Central] = corrcoef(x_spring_Central, y_spring_Central );
[R_summer_Central, pval_summer_Central] = corrcoef(x_summer_Central, y_summer_Central);


% % Central Region FALL vs SUMMER
x_fall_Central = Fall_TT_Central_BSea_daily.MASTER_Winds_daytime;
y_fall_Central = Fall_TT_Central_BSea_daily.Central_AODc_openocean;

% renaming this to '*_fs' so as not to override above spring vs summer
x_summer_Central_fs = Summer_TT_Central_BSea_daily.MASTER_Winds_daytime;
y_summer_Central_fs = Summer_TT_Central_BSea_daily.Central_AODc_openocean;

% resampling x and y to exclude any entries with NaN Values
nanVals_fall_Central = ismissing(x_fall_Central) | ismissing(y_fall_Central); 
nanVals_summer_Central_fs = ismissing(x_summer_Central_fs) | ismissing(y_summer_Central_fs);

x_fall_Central = x_fall_Central(~nanVals_fall_Central);
y_fall_Central = y_fall_Central(~nanVals_fall_Central);

x_summer_Central_fs = x_summer_Central_fs(~nanVals_summer_Central_fs);
y_summer_Central_fs = y_summer_Central_fs(~nanVals_summer_Central_fs);

% Use polyfit to compute a linear regression that predicts y from x:
p_fall_Central = polyfit(x_fall_Central,y_fall_Central,1);
p_summer_Central_fs = polyfit(x_summer_Central_fs, y_summer_Central_fs,1);

% p(1) is the slope and p(2) is the intercept of the linear predictor. You can also obtain regression coefficients using the Basic Fitting UI.
% Call polyval to use p to predict y, calling the result yfit:
yfit_fall_Central = polyval(p_fall_Central,x_fall_Central);
yfit_summer_Central_fs = polyval(p_summer_Central_fs, x_summer_Central_fs);

[R_fall_Central, pval_fall_Central] = corrcoef(x_fall_Central, y_fall_Central );
[R_summer_Central_fs, pval_summer_Central_fs] = corrcoef(x_summer_Central_fs, y_summer_Central_fs);

%%
olive = rgb('olive green');
green = rgb('forest green');
orange = rgb('orange');
dark_orange = rgb('dark orange');
ice_blue = rgb('periwinkle blue');
royal_blue = rgb('royal blue');
pink = rgb('pink');

%%

make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.11 0.04], [0.04 0.05], [0.1 0.03]);
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);

if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

% Corner Winter VS Summer Scatterplot:
subplot(2,2,3);
scatter(x_winter_Corner,y_winter_Corner,20, ice_blue);
hold on
plot(x_winter_Corner,yfit_winter_Corner,'-',...
    'LineWidth', 2, 'Color',royal_blue)

hold on
scatter(x_summer_Corner, y_summer_Corner, 20, pink,'filled','MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha',.2);
hold on
plot(x_summer_Corner, yfit_summer_Corner, '-',...
    'LineWidth',2, 'Color','r')


ylim([0 0.4])
xlim ([0 25])

legend('Winter data','Winter linear fit', 'Summer data', 'Summer linear fit') 
set(gca,'FontSize', 14) 
% set(gca,'XTick',[])
ylabel('Daily AOD_C average')
% title('Region A') 
xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));

grid off

txt_winter = ['y = ' num2str(round(p_winter_Corner(1),3,'significant')) 'x + ' num2str(round(p_winter_Corner(2),3,'significant'))];
R_txt_winter = ['R = ' num2str(round(R_winter_Corner(1,2),3,'significant'))];
N_txt_winter = ['n = ' num2str(round(length(x_winter_Corner),3,'significant'))];
p_txt_winter = ['p = ' num2str(round(pval_winter_Corner(1,2),3,'significant'))];

txt_summer = ['y = ' num2str(round(p_summer_Corner(1),3,'significant')) 'x + ' num2str(round(p_summer_Corner(2),3,'significant'))];
R_txt_summer = ['R = ' num2str(round(R_summer_Corner(1,2),3,'significant'))];
N_txt_summer = ['n = ' num2str(round(length(x_summer_Corner),3,'significant'))];
p_txt_summer = ['p = ' num2str(round(pval_summer_Corner(1,2),3,'significant'))];
Region = ('(c)');

text(2, 0.38, txt_winter,'FontSize', 14,'Color', royal_blue)
text(2,0.36, R_txt_winter,'FontSize', 14, 'Color', royal_blue)
text(2,0.34,N_txt_winter,'FontSize', 14, 'Color', royal_blue)
text(2,0.32,p_txt_winter, 'FontSize', 14, 'Color', royal_blue)
text(0.25, 0.38, Region, 'FontSize', 14, 'FontWeight', 'bold')

text(2, 0.28, txt_summer,'FontSize', 14, 'Color', 'r')
text(2,0.26, R_txt_summer,'FontSize', 14, 'Color', 'r')
text(2,0.24,N_txt_summer,'FontSize', 14,'Color', 'r')
text(2,0.22, p_txt_summer, 'FontSize',14,'Color','r')

% Central SPRING VS Summer Scatterplot:
subplot(2,2,4);
% 
scatter(x_spring_Central,y_spring_Central,20, olive,'s');%'filled','MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha',.2); 
hold on
plot(x_spring_Central,yfit_spring_Central,'-',...
    'LineWidth', 2, 'Color',green)

hold on
scatter(x_summer_Central, y_summer_Central, 20, pink,'filled','MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha',.2);

hold on
plot(x_summer_Central, yfit_summer_Central, '-',...
    'LineWidth',2, 'Color','r')

hold on
scatter(x_fall_Central,y_fall_Central,20, orange,'d');%'filled','MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha',.2); 
hold on
plot(x_fall_Central,yfit_fall_Central,'-',...
    'LineWidth', 2, 'Color',dark_orange)

legend('Spring data','Spring linear fit', 'Summer data', 'Summer linear fit', 'Fall data', 'Fall linear fit','color','none') 
set(gca,'FontSize', 14)
xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
%  title('Region B') 

grid off
ylim([0 0.4])
xlim ([0 25])

set(gca,'YTick',[])
% set(gca,'XTick',[])

txt_spring = ['y = ' num2str(round(p_spring_Central(1),3,'significant')) 'x + ' num2str(round(p_spring_Central(2),3,'significant'))];
R_txt_spring = ['R = ' num2str(round(R_spring_Central(1,2),3,'significant'))];
N_txt_spring = ['n = ' num2str(round(length(x_spring_Central),3,'significant'))];
p_txt_spring = ['p = ' num2str(round(pval_spring_Central(1,2),3,'significant'))];

txt_fall = ['y = ' num2str(round(p_fall_Central(1),3,'significant')) 'x + ' num2str(round(p_fall_Central(2),3,'significant'))];
R_txt_fall = ['R = ' num2str(round(R_fall_Central(1,2),3,'significant'))];
N_txt_fall = ['n = ' num2str(round(length(x_fall_Central),3,'significant'))];
p_txt_fall = ['p = ' num2str(round(pval_fall_Central(1,2),3,'significant'))];

txt_summer = ['y = ' num2str(round(p_summer_Central(1),3,'significant')) 'x + ' num2str(round(p_summer_Central(2),3,'significant'))];
R_txt_summer = ['R = ' num2str(round(R_summer_Central(1,2),3,'significant'))];
N_txt_summer = ['n = ' num2str(round(length(x_summer_Central),3,'significant'))];
p_txt_summer = ['p = ' num2str(round(pval_summer_Central(1,2),3,'significant'))];
Region = '(d)';

text(2, 0.38, txt_spring,'FontSize', 14, 'Color', green)
text(2,0.36, R_txt_spring,'FontSize', 14, 'Color', green)
text(2,0.34,N_txt_spring,'FontSize', 14, 'Color', green)
text(2, 0.32, p_txt_spring,'FontSize', 14' ,'Color', green)
text(0.25, 0.38, Region, 'FontSize', 14, 'FontWeight', 'bold')

text(2, 0.28, txt_fall,'FontSize', 14, 'Color', dark_orange)
text(2,0.26, R_txt_fall,'FontSize', 14, 'Color', dark_orange)
text(2,0.24,N_txt_fall,'FontSize', 14, 'Color', dark_orange)
text(2,0.22,p_txt_fall,'FontSize', 14,'Color' , dark_orange)

text(16, 0.22, txt_summer,'FontSize', 14, 'Color', 'r')
text(16,0.2, R_txt_summer,'FontSize', 14, 'Color', 'r')
text(16,0.18,N_txt_summer,'FontSize', 14, 'Color', 'r')
text(16,0.16,p_txt_summer,'FontSize', 14,'Color' , 'r')


% Central FALL VS Summer Scatterplot:
% h(3) = subplot(3,2,3);
% 
% % 
% scatter(x_fall_Central,y_fall_Central,15, ice_blue,'filled','MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha',.2); 
% hold on
% plot(x_fall_Central,yfit_fall_Central,'-',...
%     'LineWidth', 2, 'Color',royal_blue)
% hold on
% scatter(x_summer_Central, y_summer_Central, 15, pink,'filled','MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha',.2);
% 
% hold on
% plot(x_summer_Central, yfit_summer_Central, '-',...
%     'LineWidth',2, 'Color','r')
% 
% legend('Fall data','Fall linear fit', 'Summer data', 'Summer linear fit') 
% set(gca,'FontSize', 14)
% % xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
% % ylabel('Daily AOD_C average')
% title('Region B: AOD_C vs. Wind Speed') 
% ylabel('Daily OD average')
% 
% ylim([0 0.4])
% xlim ([0 20])
% grid off 
% 
% txt_fall = ['y = ' num2str(round(p_fall_Central(1)) 'x + ' num2str(round(p_fall_Central(2))];
% R_txt_fall = ['R = ' num2str(round(R_fall_Central(1,2))];
% N_txt_fall = ['n = ' num2str(round(length(x_fall_Central))];
% p_txt_fall = ['p = ' num2str(round(pval_fall_Central(1,2))];
% 
% txt_summer = ['y = ' num2str(round(p_summer_Central(1)) 'x + ' num2str(round(p_summer_Central(2))];
% R_txt_summer = ['R = ' num2str(round(R_summer_Central(1,2))];
% N_txt_summer = ['n = ' num2str(round(length(x_summer_Central))];
% p_txt_summer = ['p = ' num2str(round(pval_summer_Central(1,2))];
% Region = '(c)';
% 
% 
% text(2, 0.35, txt_fall,'FontSize', 14, 'Color', royal_blue)
% text(2,0.32, R_txt_fall,'FontSize', 14, 'Color', royal_blue)
% text(2,0.29,N_txt_fall,'FontSize', 14, 'Color', royal_blue)
% text(2,0.26, p_txt_fall, 'FontSize', 14, 'Color' , royal_blue)
% text(0.25, 0.36, Region, 'FontSize', 14, 'FontWeight', 'bold')
% 
% 
% text(8, 0.35, txt_summer,'FontSize', 14, 'Color', 'r')
% text(8,0.32, R_txt_summer,'FontSize', 14, 'Color', 'r')
% text(8,0.29,N_txt_summer,'FontSize', 14, 'Color', 'r')
% text(8,0.26, p_txt_summer,'FontSize',14,'Color', 'r')


%%
% set(gcf,'PaperPositionMode','auto')
% print(gcf,'NEWWINDS_OPENOCEAN_Coarse_AOD_Summer_vs_seasons_scatterplot.png','-dpng','-r300');       %  *// 300 dpi


%% SUBPLOTTED SUMMER VS WINTER CORNER AND CENTRAL REGIONS:

% Corner Region
clear x y nanVals x_winter y_winter
x_winter_Corner = Winter_TT_Corner_BSea_daily.MASTER_Winds_nighttime;
y_winter_Corner = Winter_TT_Corner_BSea_daily.MAOD_nocleanair_openocean;
x_summer_Corner = Summer_TT_Corner_BSea_daily.MASTER_Winds_nighttime;
y_summer_Corner = Summer_TT_Corner_BSea_daily.MAOD_nocleanair_openocean;

% resampling x and y to exclude any entries with NaN Values
nanVals_winter_Corner = ismissing(x_winter_Corner) | ismissing(y_winter_Corner); 
nanVals_summer_Corner = ismissing(x_summer_Corner) | ismissing(y_summer_Corner);

x_winter_Corner = x_winter_Corner(~nanVals_winter_Corner);
y_winter_Corner = y_winter_Corner(~nanVals_winter_Corner);

x_summer_Corner = x_summer_Corner(~nanVals_summer_Corner);
y_summer_Corner = y_summer_Corner(~nanVals_summer_Corner);

% Use polyfit to compute a linear regression that predicts y from x:
p_winter_Corner = polyfit(x_winter_Corner,y_winter_Corner,1);
p_summer_Corner = polyfit(x_summer_Corner, y_summer_Corner,1);

% p(1) is the slope and p(2) is the intercept of the linear predictor. You can also obtain regression coefficients using the Basic Fitting UI.
% Call polyval to use p to predict y, calling the result yfit:
yfit_winter_Corner = polyval(p_winter_Corner,x_winter_Corner);
yfit_summer_Corner = polyval(p_summer_Corner, x_summer_Corner);

[R_winter_Corner, pval_winter_Corner] = corrcoef(x_winter_Corner, y_winter_Corner );
[R_summer_Corner, pval_summer_Corner] = corrcoef(x_summer_Corner, y_summer_Corner);


% Central Region
clear x y nanVals x_winter y_winter
x_winter_Central = Winter_TT_Central_BSea_daily.MASTER_Winds_nighttime;
y_winter_Central = Winter_TT_Central_BSea_daily.MAOD_nocleanair_openocean;
x_summer_Central = Summer_TT_Central_BSea_daily.MASTER_Winds_nighttime;
y_summer_Central = Summer_TT_Central_BSea_daily.MAOD_nocleanair_openocean;

% resampling x and y to exclude any entries with NaN Values
nanVals_winter_Central = ismissing(x_winter_Central) | ismissing(y_winter_Central); 
nanVals_summer_Central = ismissing(x_summer_Central) | ismissing(y_summer_Central);

x_winter_Central = x_winter_Central(~nanVals_winter_Central);
y_winter_Central = y_winter_Central(~nanVals_winter_Central);

x_summer_Central = x_summer_Central(~nanVals_summer_Central);
y_summer_Central = y_summer_Central(~nanVals_summer_Central);

% Use polyfit to compute a linear regression that predicts y from x:
p_winter_Central = polyfit(x_winter_Central,y_winter_Central,1);
p_summer_Central = polyfit(x_summer_Central, y_summer_Central,1);

% p(1) is the slope and p(2) is the intercept of the linear predictor. You can also obtain regression coefficients using the Basic Fitting UI.
% Call polyval to use p to predict y, calling the result yfit:
yfit_winter_Central = polyval(p_winter_Central,x_winter_Central);
yfit_summer_Central = polyval(p_summer_Central, x_summer_Central);

[R_winter_Central, pval_winter_Central] = corrcoef(x_winter_Central, y_winter_Central );
[R_summer_Central, pval_summer_Central] = corrcoef(x_summer_Central, y_summer_Central);

%%

% ice_blue = rgb('periwinkle blue');
% royal_blue = rgb('royal blue');
% pink = rgb('pink');

%%
% 
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.075 0.05], [0.08 0.05], [0.1 0.03]);
% if ~make_it_tight,  clear subplot;  end

% fig = figure; clf;

% Corner Winter VS Summer Scatterplot:
subplot(2,2,1);


scatter(x_winter_Corner,y_winter_Corner,20, ice_blue);%,'filled') ;
hold on
plot(x_winter_Corner,yfit_winter_Corner,'-',...
    'LineWidth', 2, 'Color',royal_blue)


hold on
scatter(x_summer_Corner, y_summer_Corner, 20, pink,'filled');
hold on
plot(x_summer_Corner, yfit_summer_Corner, '-',...
    'LineWidth',2, 'Color','r')

ylim([0 0.4])
xlim ([0 25])
legend('Winter data','Winter linear fit', 'Summer data', 'Summer linear fit') 
set(gca,'FontSize', 14)
set(gca, 'XTick', [])

ylabel('Daily MAOD average')
title('Region A') 


txt_winter = ['y = ' num2str(round(p_winter_Corner(1),3,'significant')) 'x + ' num2str(round(p_winter_Corner(2),3,'significant'))];
R_txt_winter = ['R = ' num2str(round(R_winter_Corner(1,2),3,'significant'))];
N_txt_winter = ['n = ' num2str(round(length(x_winter_Corner),3,'significant'))];
p_txt_winter = ['p = ' num2str(round(pval_winter_Corner(1,2),3,'significant'))];

txt_summer = ['y = ' num2str(round(p_summer_Corner(1),3,'significant')) 'x + ' num2str(round(p_summer_Corner(2),3,'significant'))];
R_txt_summer = ['R = ' num2str(round(R_summer_Corner(1,2),3,'significant'))];
N_txt_summer = ['n = ' num2str(round(length(x_summer_Corner),3,'significant'))];
p_txt_summer = ['p = ' num2str(round(pval_summer_Corner(1,2),3,'significant'))];

Region = ('(a)');

text(2, 0.38, txt_winter,'FontSize', 14, 'Color', royal_blue)
text(2, 0.36, R_txt_winter,'FontSize', 14, 'Color', royal_blue)
text(2, 0.34, N_txt_winter,'FontSize', 14, 'Color', royal_blue)
text(2, 0.32, p_txt_winter,'FontSize',14,'Color',royal_blue)

text(0.25, 0.38, Region, 'FontSize', 14, 'FontWeight', 'bold')

text(2, 0.28, txt_summer,'FontSize', 14, 'Color', 'r')
text(2,0.26, R_txt_summer,'FontSize', 14, 'Color', 'r')
text(2,0.24,N_txt_summer,'FontSize', 14, 'Color', 'r')
text(2,0.22,p_txt_summer,'FontSize',14,'Color','r')

% Central Winter VS Summer Scatterplot:
subplot(2,2,2);
% pos = get(h,'Position');
% new = mean(cellfun(@(v)v(1),pos(1:2)));
% set(h(4),'Position',[new,pos{end}(2:end)])

scatter(x_winter_Central,y_winter_Central,20, ice_blue)%,'filled') ;
hold on
plot(x_winter_Central,yfit_winter_Central,'-',...
    'LineWidth', 2, 'Color',royal_blue)
hold on
scatter(x_summer_Central, y_summer_Central, 20, pink,'filled');

hold on
plot(x_summer_Central, yfit_summer_Central, '-',...
    'LineWidth',2, 'Color','r')

ylim([0 0.4])
xlim ([0 25])

legend('Winter data','Winter linear fit', 'Summer data', 'Summer linear fit') 
set(gca,'FontSize', 14)
set(gca,'YTick',[])
set(gca, 'XTick', [])
title('Region B')
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));

txt_winter = ['y = ' num2str(round(p_winter_Central(1),3,'significant')) 'x + ' num2str(round(p_winter_Central(2),3,'significant'))];
R_txt_winter = ['R = ' num2str(round(R_winter_Central(1,2),3,'significant'))];
N_txt_winter = ['n = ' num2str(round(length(x_winter_Central),3,'significant'))];
p_txt_winter = ['p = ' num2str(round(pval_winter_Central(1,2),3,'significant'))];

txt_summer = ['y = ' num2str(round(p_summer_Central(1),3,'significant')) 'x + ' num2str(round(p_summer_Central(2),3,'significant'))];
R_txt_summer = ['R = ' num2str(round(R_summer_Central(1,2),3,'significant'))];
N_txt_summer = ['n = ' num2str(round(length(x_summer_Central),3,'significant'))];
p_txt_summer = ['p = ' num2str(round(pval_summer_Central(1,2),3,'significant'))];

Region = '(b)';

text(2, 0.38, txt_winter,'FontSize', 14, 'Color', royal_blue)
text(2, 0.36, R_txt_winter,'FontSize', 14, 'Color', royal_blue)
text(2, 0.34, N_txt_winter,'FontSize', 14, 'Color', royal_blue)
text(2, 0.32, p_txt_winter,'FontSize',14,'Color',royal_blue)

text(0.25, 0.38, Region, 'FontSize', 14, 'FontWeight', 'bold')

text(2, 0.28, txt_summer,'FontSize', 14, 'Color', 'r')
text(2,0.26, R_txt_summer,'FontSize', 14, 'Color', 'r')
text(2,0.24,N_txt_summer,'FontSize', 14, 'Color', 'r')
text(2,0.22,p_txt_summer,'FontSize',14,'Color','r')





%%
set(gcf,'PaperPositionMode','auto')
print(gcf,'DIURNALWINDS_OPENOCEAN_MAOD_AOD_fine_allregions_Summer_vs_Winter_scatterplot.png_V2','-dpng','-r300');      


