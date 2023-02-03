
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

%% Scatter plot, All Regions, daily data: AODt & AODf only
% SUBPLOTTED

%%

x_Corner_Wind_Daily_forAODt = Corner_complete_TT_daily.MASTER_Winds_daytime;
y_Corner_AODt_Daily = Corner_complete_TT_daily.Corner_AOD_total_openocean;

nanVals_Corner_AODt_Wind = ismissing(x_Corner_Wind_Daily_forAODt) | ismissing(y_Corner_AODt_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_Daily_forAODt = x_Corner_Wind_Daily_forAODt(~nanVals_Corner_AODt_Wind);
y_Corner_AODt_Daily = y_Corner_AODt_Daily(~nanVals_Corner_AODt_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_Daily_AODt_Wind = polyfit(x_Corner_Wind_Daily_forAODt,y_Corner_AODt_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_Daily_AODt_Wind = polyval(p_Corner_Daily_AODt_Wind,x_Corner_Wind_Daily_forAODt);

[R_Corner_Daily_AODt_Wind, pval_Corner_Daily_AODt_Wind] = corrcoef(x_Corner_Wind_Daily_forAODt, y_Corner_AODt_Daily );


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_Wind_Daily_forAODt = Central_complete_TT_daily.MASTER_Winds_daytime;
y_Central_AODt_Daily = Central_complete_TT_daily.Central_AOD_total_openocean;

nanVals_Central_AODt_Wind = ismissing(x_Central_Wind_Daily_forAODt) | ismissing(y_Central_AODt_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_Daily_forAODt = x_Central_Wind_Daily_forAODt(~nanVals_Central_AODt_Wind);
y_Central_AODt_Daily = y_Central_AODt_Daily(~nanVals_Central_AODt_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_Daily_AODt_Wind = polyfit(x_Central_Wind_Daily_forAODt,y_Central_AODt_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_Daily_AODt_Wind = polyval(p_Central_Daily_AODt_Wind,x_Central_Wind_Daily_forAODt);

[R_Central_Daily_AODt_Wind, pval_Central_Daily_AODt_Wind] = corrcoef(x_Central_Wind_Daily_forAODt, y_Central_AODt_Daily );


% Southernmost daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

% load count.dat
x_Southernmost_Wind_Daily_forAODt = Southernmost_complete_TT_daily.MASTER_Winds_daytime;
y_Southernmost_AODt_Daily = Southernmost_complete_TT_daily.Southernmost_AOD_total_openocean;

nanVals_Southernmost_AODt_Wind = ismissing(x_Southernmost_Wind_Daily_forAODt) | ismissing(y_Southernmost_AODt_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_Daily_forAODt = x_Southernmost_Wind_Daily_forAODt(~nanVals_Southernmost_AODt_Wind);
y_Southernmost_AODt_Daily = y_Southernmost_AODt_Daily(~nanVals_Southernmost_AODt_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_Daily_AODt_Wind = polyfit(x_Southernmost_Wind_Daily_forAODt,y_Southernmost_AODt_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_Daily_AODt_Wind = polyval(p_Southernmost_Daily_AODt_Wind,x_Southernmost_Wind_Daily_forAODt);

[R_Southernmost_Daily_AODt_Wind, pval_Southernmost_Daily_AODt_Wind] = corrcoef(x_Southernmost_Wind_Daily_forAODt, y_Southernmost_AODt_Daily );

%%

x_Corner_Wind_Daily_forAODf = Corner_complete_TT_daily.MASTER_Winds_daytime;
y_Corner_AODf_Daily = Corner_complete_TT_daily.Corner_AOD_fine_openocean;

nanVals_Corner_AODf_Wind = ismissing(x_Corner_Wind_Daily_forAODf) | ismissing(y_Corner_AODf_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_Daily_forAODf = x_Corner_Wind_Daily_forAODf(~nanVals_Corner_AODf_Wind);
y_Corner_AODf_Daily = y_Corner_AODf_Daily(~nanVals_Corner_AODf_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_Daily_AODf_Wind = polyfit(x_Corner_Wind_Daily_forAODf,y_Corner_AODf_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_Daily_AODf_Wind = polyval(p_Corner_Daily_AODf_Wind,x_Corner_Wind_Daily_forAODf);

[R_Corner_Daily_AODf_Wind, pval_Corner_Daily_AODf_Wind] = corrcoef(x_Corner_Wind_Daily_forAODf, y_Corner_AODf_Daily );


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_Wind_Daily_forAODf = Central_complete_TT_daily.MASTER_Winds_daytime;
y_Central_AODf_Daily = Central_complete_TT_daily.Central_AOD_fine_openocean;

nanVals_Central_AODf_Wind = ismissing(x_Central_Wind_Daily_forAODf) | ismissing(y_Central_AODf_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_Daily_forAODf = x_Central_Wind_Daily_forAODf(~nanVals_Central_AODf_Wind);
y_Central_AODf_Daily = y_Central_AODf_Daily(~nanVals_Central_AODf_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_Daily_AODf_Wind = polyfit(x_Central_Wind_Daily_forAODf,y_Central_AODf_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_Daily_AODf_Wind = polyval(p_Central_Daily_AODf_Wind,x_Central_Wind_Daily_forAODf);

[R_Central_Daily_AODf_Wind, pval_Central_Daily_AODf_Wind] = corrcoef(x_Central_Wind_Daily_forAODf, y_Central_AODf_Daily );


% Southernmost daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_Wind_Daily_forAODf = Southernmost_complete_TT_daily.MASTER_Winds_daytime;
y_Southernmost_AODf_Daily = Southernmost_complete_TT_daily.Southernmost_AOD_fine_openocean;

nanVals_Southernmost_AODf_Wind = ismissing(x_Southernmost_Wind_Daily_forAODf) | ismissing(y_Southernmost_AODf_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_Daily_forAODf = x_Southernmost_Wind_Daily_forAODf(~nanVals_Southernmost_AODf_Wind);
y_Southernmost_AODf_Daily = y_Southernmost_AODf_Daily(~nanVals_Southernmost_AODf_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_Daily_AODf_Wind = polyfit(x_Southernmost_Wind_Daily_forAODf,y_Southernmost_AODf_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_Daily_AODf_Wind = polyval(p_Southernmost_Daily_AODf_Wind,x_Southernmost_Wind_Daily_forAODf);

[R_Southernmost_Daily_AODf_Wind, pval_Southernmost_Daily_AODf_Wind] = corrcoef(x_Southernmost_Wind_Daily_forAODf, y_Southernmost_AODf_Daily );


%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

%Corner Daily Scatterplot
subplot(2,3,1)
scatter(x_Corner_Wind_Daily_forAODt,y_Corner_AODt_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_Daily_forAODt,yfit_Corner_Daily_AODt_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 30])
set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)

title('Region A') 

txt = ['y = ' num2str(round(p_Corner_Daily_AODt_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_AODt_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_AODt_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_AODt_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_Daily_forAODt), 3, 'significant'))];
Region = ('AOD_T');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10, 0.36, p_txt,'FontSize', 15)
text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region



% Corner Daily Scatterplot:
subplot(2,3,4)
scatter(x_Corner_Wind_Daily_forAODf,y_Corner_AODf_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_Daily_forAODf,yfit_Corner_Daily_AODf_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 30])
set(gca,'FontSize', 15)


txt = ['y = ' num2str(round(p_Corner_Daily_AODf_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_AODf_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_AODf_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_AODf_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_Daily_forAODf),3,'significant'))];
Region = ('AOD_f');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10, 0.36, p_txt,'FontSize', 15)
text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% Central Daily Scatterplot
subplot(2,3,2)
scatter(x_Central_Wind_Daily_forAODt,y_Central_AODt_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_Daily_forAODt,yfit_Central_Daily_AODt_Wind,'-',...
    'LineWidth', 2)


ylim([0 0.5])
xlim([0 30])
set(gca,'xtick',[],'xticklabel',[], 'ytick',[], 'FontSize', 15)

title('Region B')

txt = ['y = ' num2str(round(p_Central_Daily_AODt_Wind(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily_AODt_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_AODt_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_AODt_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_Daily_forAODt),3,'significant'))];
Region = ('AOD_T');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10, 0.36, p_txt,'FontSize', 15)

text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region

% Central Daily Scatterplot
subplot(2,3,5)
scatter(x_Central_Wind_Daily_forAODf,y_Central_AODf_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_Daily_forAODf,yfit_Central_Daily_AODf_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 30])
set(gca, 'ytick',[], 'FontSize', 15)



txt = ['y = ' num2str(round(p_Central_Daily_AODf_Wind(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily_AODf_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_AODf_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_AODf_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_Daily_forAODf),3,'significant'))];
Region = ('AOD_f');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10, 0.36, p_txt,'FontSize', 15)

text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% Southernmost Daily Scatterplot:
subplot(2,3,3)
scatter(x_Southernmost_Wind_Daily_forAODt,y_Southernmost_AODt_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Southernmost_Wind_Daily_forAODt,yfit_Southernmost_Daily_AODt_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 30])
set(gca,'xtick',[],'xticklabel',[],'ytick',[], 'FontSize', 15)
 
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
title('Region C')
txt = ['y = ' num2str(round(p_Southernmost_Daily_AODt_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_AODt_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_AODt_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_AODt_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_Daily_forAODt),3,'significant'))];
Region = ('AOD_T');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10,0.36,p_txt, 'FontSize', 15)
text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
clear txt R_txt N_txt Region



% Southernmost Daily Scatterplot:
subplot(2,3,6)
scatter(x_Southernmost_Wind_Daily_forAODf,y_Southernmost_AODf_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;

hold on
plot(x_Southernmost_Wind_Daily_forAODf,yfit_Southernmost_Daily_AODf_Wind,'-',...
    'LineWidth', 2)
ylim([0 0.5])
xlim([0 30])
set(gca,'ytick',[],'FontSize', 15)
 
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));

txt = ['y = ' num2str(round(p_Southernmost_Daily_AODf_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_AODf_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_AODf_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_AODf_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_Daily_forAODf),3,'significant'))];
Region = ('AOD_f');

text(10, 0.45, txt,'FontSize', 15)
text(10,0.42, R_txt,'FontSize', 15)
text(10,0.39,N_txt,'FontSize', 15)
text(10,0.36,p_txt, 'FontSize', 15)
text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
clear txt R_txt N_txt Region



%%

set(gcf,'PaperPositionMode','auto')
print(gcf,'DIURNALWINDS_OPENOCEAN_AODt_AODf_allregions_daily_scatterplot.png','-dpng','-r300');       %  *// 300 dpi


%% 
% ANNUAL DATA SCATTERPLOTS
clear;
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_vars
load('Central_complete_TT_annual.mat')
load('Southernmost_complete_TT_annual.mat')
load('Corner_complete_TT_annual.mat')

%% Scatter plot, annual data:
% I've kept the variable names the same so that plotting is easy using the
% same script as above

clear x y nanVals p yfit R
x_Corner_Wind_Daily_forAODc = Corner_complete_TT_annual.MASTER_Winds_daytime;
y_Corner_AODc_Daily = Corner_complete_TT_annual.Corner_AODc_openocean;

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

x_Central_Wind_Daily_forAODc = Central_complete_TT_annual.MASTER_Winds_daytime;
y_Central_AODc_Daily = Central_complete_TT_annual.Central_AOD_coarse_daily;

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
x_Southernmost_Wind_Daily_forAODc = Southernmost_complete_TT_annual.MASTER_Winds_daytime;
y_Southernmost_AODc_Daily = Southernmost_complete_TT_annual.Southernmost_AOD_coarse_daily;
    
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

x_Corner_Wind_Daily_forMAOD = Corner_complete_TT_annual.MASTER_Winds_nighttime;
y_Corner_MAOD_Daily = Corner_complete_TT_annual.MAOD_nocleanair_openocean;

nanVals_Corner_MAOD_Wind = ismissing(x_Corner_Wind_Daily_forMAOD) | ismissing(y_Corner_MAOD_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_Daily_forMAOD = x_Corner_Wind_Daily_forMAOD(~nanVals_Corner_MAOD_Wind);
y_Corner_MAOD_Daily = y_Corner_MAOD_Daily(~nanVals_Corner_MAOD_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_Daily_MAOD_Wind = polyfit(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_Daily_MAOD_Wind = polyval(p_Corner_Daily_MAOD_Wind,x_Corner_Wind_Daily_forMAOD);

[R_Corner_Daily_MAOD_Wind, pval_Corner_Daily_MAOD_Wind] = corrcoef(x_Corner_Wind_Daily_forMAOD, y_Corner_MAOD_Daily );


% Central yearly scatterplot



x_Central_Wind_Daily_forMAOD = Central_complete_TT_annual.MASTER_Winds_nighttime;
y_Central_MAOD_Daily = Central_complete_TT_annual.MAOD_nocleanair_Central_TT_MAOD_Wind_annual_nighttime;

nanVals_Central_MAOD_Wind = ismissing(x_Central_Wind_Daily_forMAOD) | ismissing(y_Central_MAOD_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_Daily_forMAOD = x_Central_Wind_Daily_forMAOD(~nanVals_Central_MAOD_Wind);
y_Central_MAOD_Daily = y_Central_MAOD_Daily(~nanVals_Central_MAOD_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_Daily_MAOD_Wind = polyfit(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_Daily = polyval(p_Central_Daily_MAOD_Wind,x_Central_Wind_Daily_forMAOD);

[R_Central_Daily_MAOD_Wind, pval_Central_Daily_MAOD_Wind] = corrcoef(x_Central_Wind_Daily_forMAOD, y_Central_MAOD_Daily );


% Southernmost yearly scatterplot

clear x y
% load count.dat
x_Southernmost_Wind_Daily_forMAOD = Southernmost_complete_TT_annual.MASTER_Winds_nighttime;
y_Southernmost_MAOD_Daily = Southernmost_complete_TT_annual.MAOD_nocleanair_Southernmost_TT_MAOD_Wind_annual_nighttime;

nanVals_Southernmost_MAOD_Wind = ismissing(x_Southernmost_Wind_Daily_forMAOD) | ismissing(y_Southernmost_MAOD_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_Daily_forMAOD = x_Southernmost_Wind_Daily_forMAOD(~nanVals_Southernmost_MAOD_Wind);
y_Southernmost_MAOD_Daily = y_Southernmost_MAOD_Daily(~nanVals_Southernmost_MAOD_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_Daily_MAOD_Wind = polyfit(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_Daily_MAOD_Wind = polyval(p_Southernmost_Daily_MAOD_Wind,x_Southernmost_Wind_Daily_forMAOD);

[R_Southernmost_Daily_MAOD_Wind, pval_Southernmost_Daily_MAOD_Wind] = corrcoef(x_Southernmost_Wind_Daily_forMAOD, y_Southernmost_MAOD_Daily );

%% Scatter plot, annual data: AOD_t

clear x y nanVals p yfit R
x_Corner_Wind_Daily_forAODt = Corner_complete_TT_annual.MASTER_Winds_daytime;
y_Corner_AODt_Daily = Corner_complete_TT_annual.Corner_AOD_total_openocean;

nanVals_Corner_AODt_Wind = ismissing(x_Corner_Wind_Daily_forAODt) | ismissing(y_Corner_AODt_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_Daily_forAODt = x_Corner_Wind_Daily_forAODt(~nanVals_Corner_AODt_Wind);
y_Corner_AODt_Daily = y_Corner_AODt_Daily(~nanVals_Corner_AODt_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_Daily_AODt_Wind = polyfit(x_Corner_Wind_Daily_forAODt,y_Corner_AODt_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_Daily_AODt_Wind = polyval(p_Corner_Daily_AODt_Wind,x_Corner_Wind_Daily_forAODt);

[R_Corner_Daily_AODt_Wind, pval_Corner_Daily_AODt_Wind] = corrcoef(x_Corner_Wind_Daily_forAODt, y_Corner_AODt_Daily );


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_Wind_Daily_forAODt = Central_complete_TT_annual.MASTER_Winds_daytime;
y_Central_AODt_Daily = Central_complete_TT_annual.Central_AOD_total_daily;

nanVals_Central_AODt_Wind = ismissing(x_Central_Wind_Daily_forAODt) | ismissing(y_Central_AODt_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_Daily_forAODt = x_Central_Wind_Daily_forAODt(~nanVals_Central_AODt_Wind);
y_Central_AODt_Daily = y_Central_AODt_Daily(~nanVals_Central_AODt_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_Daily_AODt_Wind = polyfit(x_Central_Wind_Daily_forAODt,y_Central_AODt_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_Daily_AODt_Wind = polyval(p_Central_Daily_AODt_Wind,x_Central_Wind_Daily_forAODt);

[R_Central_Daily_AODt_Wind, pval_Central_Daily_AODt_Wind] = corrcoef(x_Central_Wind_Daily_forAODt, y_Central_AODt_Daily );


% Southernmost daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_Wind_Daily_forAODt = Southernmost_complete_TT_annual.MASTER_Winds_daytime;
y_Southernmost_AODt_Daily = Southernmost_complete_TT_annual.Southernmost_AOD_total_daily;
    
nanVals_Southernmost_AODt_Wind = ismissing(x_Southernmost_Wind_Daily_forAODt) | ismissing(y_Southernmost_AODt_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_Daily_forAODt = x_Southernmost_Wind_Daily_forAODt(~nanVals_Southernmost_AODt_Wind);
y_Southernmost_AODt_Daily = y_Southernmost_AODt_Daily(~nanVals_Southernmost_AODt_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_Daily_AODt_Wind = polyfit(x_Southernmost_Wind_Daily_forAODt,y_Southernmost_AODt_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_Daily_AODt_Wind = polyval(p_Southernmost_Daily_AODt_Wind,x_Southernmost_Wind_Daily_forAODt);

[R_Southernmost_Daily_AODt_Wind, pval_Southernmost_Daily_AODt_Wind] = corrcoef(x_Southernmost_Wind_Daily_forAODt, y_Southernmost_AODt_Daily );



%% Scatter plot, annual data: AOD_f

clear x y nanVals p yfit R
x_Corner_Wind_Daily_forAODf = Corner_complete_TT_annual.MASTER_Winds_daytime;
y_Corner_AODf_Daily = Corner_complete_TT_annual.Corner_AOD_fine_openocean;

nanVals_Corner_AODf_Wind = ismissing(x_Corner_Wind_Daily_forAODf) | ismissing(y_Corner_AODf_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_Daily_forAODf = x_Corner_Wind_Daily_forAODf(~nanVals_Corner_AODf_Wind);
y_Corner_AODf_Daily = y_Corner_AODf_Daily(~nanVals_Corner_AODf_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_Daily_AODf_Wind = polyfit(x_Corner_Wind_Daily_forAODf,y_Corner_AODf_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_Daily_AODf_Wind = polyval(p_Corner_Daily_AODf_Wind,x_Corner_Wind_Daily_forAODf);

[R_Corner_Daily_AODf_Wind, pval_Corner_Daily_AODf_Wind] = corrcoef(x_Corner_Wind_Daily_forAODf, y_Corner_AODf_Daily );


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_Wind_Daily_forAODf = Central_complete_TT_annual.MASTER_Winds_daytime;
y_Central_AODf_Daily = Central_complete_TT_annual.Central_AOD_fine_daily;

nanVals_Central_AODf_Wind = ismissing(x_Central_Wind_Daily_forAODf) | ismissing(y_Central_AODf_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_Daily_forAODf = x_Central_Wind_Daily_forAODf(~nanVals_Central_AODf_Wind);
y_Central_AODf_Daily = y_Central_AODf_Daily(~nanVals_Central_AODf_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_Daily_AODf_Wind = polyfit(x_Central_Wind_Daily_forAODf,y_Central_AODf_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_Daily_AODf_Wind = polyval(p_Central_Daily_AODf_Wind,x_Central_Wind_Daily_forAODf);

[R_Central_Daily_AODf_Wind, pval_Central_Daily_AODf_Wind] = corrcoef(x_Central_Wind_Daily_forAODf, y_Central_AODf_Daily );


% Southernmost daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_Wind_Daily_forAODf = Southernmost_complete_TT_annual.MASTER_Winds_daytime;
y_Southernmost_AODf_Daily = Southernmost_complete_TT_annual.Southernmost_AOD_fine_daily;
    
nanVals_Southernmost_AODf_Wind = ismissing(x_Southernmost_Wind_Daily_forAODf) | ismissing(y_Southernmost_AODf_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_Daily_forAODf = x_Southernmost_Wind_Daily_forAODf(~nanVals_Southernmost_AODf_Wind);
y_Southernmost_AODf_Daily = y_Southernmost_AODf_Daily(~nanVals_Southernmost_AODf_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_Daily_AODf_Wind = polyfit(x_Southernmost_Wind_Daily_forAODf,y_Southernmost_AODf_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_Daily_AODf_Wind = polyval(p_Southernmost_Daily_AODf_Wind,x_Southernmost_Wind_Daily_forAODf);

[R_Southernmost_Daily_AODf_Wind, pval_Southernmost_Daily_AODf_Wind] = corrcoef(x_Southernmost_Wind_Daily_forAODf, y_Southernmost_AODf_Daily );

%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

%Corner Daily Scatterplot
subplot(4,3,1)
scatter(x_Corner_Wind_Daily_forMAOD,y_Corner_MAOD_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_Daily_forMAOD,yfit_Corner_Daily_MAOD_Wind,'-',...
    'LineWidth', 2)


ylim([0 0.12])
xlim([5.5 12])  
% set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)
set(gca, 'FontSize', 15)

title('(a) Region A') 

txt = ['y = ' num2str(round(p_Corner_Daily_MAOD_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_MAOD_Wind(2), 3, 'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_MAOD_Wind(1,2), 3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_MAOD_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_Daily_forMAOD),3,'significant'))];
Region = ('MAOD');

text(7, 0.1, txt,'FontSize', 15)
text(7, 0.09, R_txt,'FontSize', 15)
text(7, 0.08,N_txt,'FontSize', 15)
text(7, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold') 

clear txt R_txt N_txt Region



% Corner Daily Scatterplot:
subplot(4,3,4)
scatter(x_Corner_Wind_Daily_forAODc,y_Corner_AODc_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_Daily_forAODc,yfit_Corner_Daily_AODc_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
set(gca,'FontSize', 15)


txt = ['y = ' num2str(round(p_Corner_Daily_AODc_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_AODc_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_AODc_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_AODc_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_Daily_forAODc),3,'significant'))];
Region = ('AOD_C');



text(7, 0.1, txt,'FontSize', 15)
text(7, 0.09, R_txt,'FontSize', 15)
text(7, 0.08,N_txt,'FontSize', 15)
text(7, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold') 

clear txt R_txt N_txt Region



% Corner Daily Scatterplot:
subplot(4,3,7)
scatter(x_Corner_Wind_Daily_forAODt,y_Corner_AODt_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_Daily_forAODt,yfit_Corner_Daily_AODt_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
set(gca,'FontSize', 15)


txt = ['y = ' num2str(round(p_Corner_Daily_AODt_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_AODt_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_AODt_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_AODt_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_Daily_forAODt),3,'significant'))];
Region = ('AOD_T');



text(7, 0.1, txt,'FontSize', 15)
text(7, 0.09, R_txt,'FontSize', 15)
text(7, 0.08,N_txt,'FontSize', 15)
text(7, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')  

clear txt R_txt N_txt Region



% Corner Daily Scatterplot:
subplot(4,3,10)
scatter(x_Corner_Wind_Daily_forAODf,y_Corner_AODf_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_Daily_forAODf,yfit_Corner_Daily_AODf_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
set(gca,'FontSize', 15)


txt = ['y = ' num2str(round(p_Corner_Daily_AODf_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_AODf_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_AODf_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_AODf_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_Daily_forAODf),3,'significant'))];
Region = ('AOD_f');



text(7, 0.1, txt,'FontSize', 15)
text(7, 0.09, R_txt,'FontSize', 15)
text(7, 0.08,N_txt,'FontSize', 15)
text(7, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')  

clear txt R_txt N_txt Region



% Central Daily Scatterplot
subplot(4,3,2)
scatter(x_Central_Wind_Daily_forMAOD,y_Central_MAOD_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_Daily_forMAOD,yfit_Central_Daily,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
% set(gca,'xtick',[],'xticklabel',[], 'ytick',[], 'FontSize', 15)
set(gca, 'ytick',[], 'FontSize', 15)

title('(b) Region B')

txt = ['y = ' num2str(round(p_Central_Daily_MAOD_Wind(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily_MAOD_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_MAOD_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_MAOD_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_Daily_forMAOD),3,'significant'))];
Region = ('MAOD');



text(8.2, 0.1, txt,'FontSize', 15)
text(8.2, 0.09, R_txt,'FontSize', 15)
text(8.2, 0.08,N_txt,'FontSize', 15)
text(8.2, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')  


clear txt R_txt N_txt Region

% Central Daily Scatterplot
subplot(4,3,5)
scatter(x_Central_Wind_Daily_forAODc,y_Central_AODc_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_Daily_forAODc,yfit_Central_Daily_AODc_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
set(gca, 'ytick',[], 'FontSize', 15)



txt = ['y = ' num2str(round(p_Central_Daily_AODc_Wind(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily_AODc_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_AODc_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_AODc_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_Daily_forAODc),3,'significant'))];
Region = ('AOD_C');


text(8.2, 0.1, txt,'FontSize', 15)
text(8.2, 0.09, R_txt,'FontSize', 15)
text(8.2, 0.08,N_txt,'FontSize', 15)
text(8.2, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')  

clear txt R_txt N_txt Region


% Central Daily Scatterplot
subplot(4,3,8)
scatter(x_Central_Wind_Daily_forAODt,y_Central_AODt_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_Daily_forAODt,yfit_Central_Daily_AODt_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
set(gca, 'ytick',[], 'FontSize', 15)



txt = ['y = ' num2str(round(p_Central_Daily_AODt_Wind(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily_AODt_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_AODt_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_AODt_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_Daily_forAODt),3,'significant'))];
Region = ('AOD_T');


text(8.2, 0.1, txt,'FontSize', 15)
text(8.2, 0.09, R_txt,'FontSize', 15)
text(8.2, 0.08,N_txt,'FontSize', 15)
text(8.2, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')  

clear txt R_txt N_txt Region



% Central Daily Scatterplot
subplot(4,3,11)
scatter(x_Central_Wind_Daily_forAODf,y_Central_AODf_Daily,20,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_Daily_forAODf,yfit_Central_Daily_AODf_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
set(gca, 'ytick',[], 'FontSize', 15)



txt = ['y = ' num2str(round(p_Central_Daily_AODf_Wind(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily_AODf_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_AODf_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_AODf_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_Daily_forAODf),3,'significant'))];
Region = ('AOD_f');


text(8.2, 0.1, txt,'FontSize', 15)
text(8.2, 0.09, R_txt,'FontSize', 15)
text(8.2, 0.08,N_txt,'FontSize', 15)
text(8.2, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')  

clear txt R_txt N_txt Region


% Southernmost Daily Scatterplot:
subplot(4,3,3)
scatter(x_Southernmost_Wind_Daily_forMAOD,y_Southernmost_MAOD_Daily, 20, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Southernmost_Wind_Daily_forMAOD,yfit_Southernmost_Daily_MAOD_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
% set(gca,'xtick',[],'xticklabel',[],'ytick',[], 'FontSize', 15)
set(gca, 'ytick',[], 'FontSize', 15)

% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
title('(c) Region C')
txt = ['y = ' num2str(round(p_Southernmost_Daily_MAOD_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_MAOD_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_MAOD_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_MAOD_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_Daily_forMAOD),3,'significant'))];
Region = ('MAOD');


text(8, 0.1, txt,'FontSize', 15)
text(8, 0.09, R_txt,'FontSize', 15)
text(8, 0.08,N_txt,'FontSize', 15)
text(8, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')  

clear txt R_txt N_txt Region



% Southernmost Daily Scatterplot:
subplot(4,3,6)
scatter(x_Southernmost_Wind_Daily_forAODc,y_Southernmost_AODc_Daily, 20, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;

hold on
plot(x_Southernmost_Wind_Daily_forAODc,yfit_Southernmost_Daily_AODc_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
set(gca,'ytick',[],'FontSize', 15)
 
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));

txt = ['y = ' num2str(round(p_Southernmost_Daily_AODc_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_AODc_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_AODc_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_AODc_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_Daily_forAODc),3,'significant'))];
Region = ('AOD_C');



text(8, 0.1, txt,'FontSize', 15)
text(8, 0.09, R_txt,'FontSize', 15)
text(8, 0.08,N_txt,'FontSize', 15)
text(8, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')  

clear txt R_txt N_txt Region



% Southernmost Daily Scatterplot:
subplot(4,3,9)
scatter(x_Southernmost_Wind_Daily_forAODt,y_Southernmost_AODt_Daily, 20, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;

hold on
plot(x_Southernmost_Wind_Daily_forAODt,yfit_Southernmost_Daily_AODt_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
set(gca,'ytick',[],'FontSize', 15)
 
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));

txt = ['y = ' num2str(round(p_Southernmost_Daily_AODt_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_AODt_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_AODt_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_AODt_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_Daily_forAODt),3,'significant'))];
Region = ('AOD_T');



text(8, 0.1, txt,'FontSize', 15)
text(8, 0.09, R_txt,'FontSize', 15)
text(8, 0.08,N_txt,'FontSize', 15)
text(8, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')  

clear txt R_txt N_txt Region


% Southernmost Daily Scatterplot:
subplot(4,3,12)
scatter(x_Southernmost_Wind_Daily_forAODf,y_Southernmost_AODf_Daily, 20, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;

hold on
plot(x_Southernmost_Wind_Daily_forAODf,yfit_Southernmost_Daily_AODf_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.12])
xlim([5.5 12])  
set(gca,'ytick',[],'FontSize', 15)
 
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));

txt = ['y = ' num2str(round(p_Southernmost_Daily_AODf_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_AODf_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_AODf_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_AODf_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_Daily_forAODf),3,'significant'))];
Region = ('AOD_f');



text(8, 0.1, txt,'FontSize', 15)
text(8, 0.09, R_txt,'FontSize', 15)
text(8, 0.08,N_txt,'FontSize', 15)
text(8, 0.07, p_txt,'FontSize', 15)
text(5.8, 0.11, Region, 'FontSize', 15, 'FontWeight', 'bold')  

clear txt R_txt N_txt Region


%%
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures/Diurnal_winds
set(gcf,'PaperPositionMode','auto')
print(gcf,'DIURNALWINDS_OPENOCEAN_Allvars_allregions_annual_scatterplot.png','-dpng','-r300');       %  *// 300 dpi


%%


% OPEN OCEAN PIXELS ONLY 

%% SUBPLOTTED WINTER VS SUMMER CORNER AND CENTRAL REGIONS:

% I toggled through these, changing the variables from
% Corner_AODc_openocean to Corner_AOD_total_openocean, and so on.

% Corner Region
clear x y nanVals x_winter y_winter
x_winter_Corner = Spring_TT_Corner_BSea_daily.MASTER_Winds_daytime;
y_winter_Corner = Spring_TT_Corner_BSea_daily.Corner_AOD_fine_openocean;
x_summer_Corner = Summer_TT_Corner_BSea_daily.MASTER_Winds_daytime;
y_summer_Corner = Summer_TT_Corner_BSea_daily.Corner_AOD_fine_openocean;

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
y_spring_Central = Spring_TT_Central_BSea_daily.Central_AOD_fine_openocean;
x_summer_Central = Summer_TT_Central_BSea_daily.MASTER_Winds_daytime;
y_summer_Central = Summer_TT_Central_BSea_daily.Central_AOD_fine_openocean;

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
y_fall_Central = Fall_TT_Central_BSea_daily.Central_AOD_fine_openocean;

% renaming this to '*_fs' so as not to override above spring vs summer
x_summer_Central_fs = Summer_TT_Central_BSea_daily.MASTER_Winds_daytime;
y_summer_Central_fs = Summer_TT_Central_BSea_daily.Central_AOD_fine_openocean;

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
xlim ([0 30])

legend('Winter data','Winter linear fit', 'Summer data', 'Summer linear fit') 
set(gca,'FontSize', 14) 
% set(gca,'XTick',[])
ylabel('Daily AOD_f average')
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
xlim ([0 30])

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

text(18, 0.2, txt_summer,'FontSize', 14, 'Color', 'r')
text(18,0.18, R_txt_summer,'FontSize', 14, 'Color', 'r')
text(18,0.16,N_txt_summer,'FontSize', 14, 'Color', 'r')
text(18,0.14,p_txt_summer,'FontSize', 14,'Color' , 'r')



%% SUBPLOTTED SUMMER VS WINTER CORNER AND CENTRAL REGIONS:

% Corner Region
clear x y nanVals x_winter y_winter
x_winter_Corner = Spring_TT_Corner_BSea_daily.MASTER_Winds_daytime;
y_winter_Corner = Spring_TT_Corner_BSea_daily.Corner_AOD_total_openocean;
x_summer_Corner = Summer_TT_Corner_BSea_daily.MASTER_Winds_daytime;
y_summer_Corner = Summer_TT_Corner_BSea_daily.Corner_AOD_total_openocean;

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
y_spring_Central = Spring_TT_Central_BSea_daily.Central_AOD_total_openocean;
x_summer_Central = Summer_TT_Central_BSea_daily.MASTER_Winds_daytime;
y_summer_Central = Summer_TT_Central_BSea_daily.Central_AOD_total_openocean;

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
y_fall_Central = Fall_TT_Central_BSea_daily.Central_AOD_total_openocean;

% renaming this to '*_fs' so as not to override above spring vs summer
x_summer_Central_fs = Summer_TT_Central_BSea_daily.MASTER_Winds_daytime;
y_summer_Central_fs = Summer_TT_Central_BSea_daily.Central_AOD_total_openocean;

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

hold on;
% Corner Winter VS Summer Scatterplot:
subplot(2,2,1);
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
xlim ([0 30])

legend('Winter data','Winter linear fit', 'Summer data', 'Summer linear fit') 
set(gca,'FontSize', 14) 
set(gca,'XTick',[])
ylabel('Daily AOD_T average')
% title('Region A') 
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));

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
subplot(2,2,2);
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
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
%  title('Region B') 

grid off
ylim([0 0.4])
xlim ([0 30])

set(gca,'YTick',[])
set(gca,'XTick',[])

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

text(16, 0.11, txt_summer,'FontSize', 14, 'Color', 'r')
text(16,0.09, R_txt_summer,'FontSize', 14, 'Color', 'r')
text(16,0.07,N_txt_summer,'FontSize', 14, 'Color', 'r')
text(16,0.05,p_txt_summer,'FontSize', 14,'Color' , 'r')



%%
set(gcf,'PaperPositionMode','auto')
print(gcf,'DIURNALWINDS_OPENOCEAN_AODT_AODf_allregions_Summer_vs_otherseasons_scatterplot.png','-dpng','-r300');      


