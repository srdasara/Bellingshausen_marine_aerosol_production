%% My attempt to make a seasonal scatterplot (with SST overlayed)! 


clear;
cd /Volumes/'G-DRIVE mobile USB'/SST_matfiles/
load('SST_Seasonal_TT.mat')

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_vars
load('Central_complete_TT_daily.mat')
load('Corner_complete_TT_daily.mat')
load('Southernmost_complete_TT_daily.mat')


%%
% Trying out seasonal time range so that I can then step across 'retime'
% and specify 'quarterly' as the timestep in that function. See below.

SS_season = timerange('2007-03-01', '2018-11-30', 'closed'); %closed includes startime and endtime values
dt = calmonths(3);

Southernmost_TT_SEASONAL = Southernmost_complete_TT_daily(SS_season,:);
Southernmost_TT_SEASONAL = retime(Southernmost_TT_SEASONAL,'regular',@nanmean, 'TimeStep',dt);

Central_TT_SEASONAL = Central_complete_TT_daily(SS_season,:);
Central_TT_SEASONAL = retime(Central_TT_SEASONAL, 'regular', @nanmean, 'TimeStep', dt);

Corner_TT_SEASONAL = Corner_complete_TT_daily(SS_season, :);
Corner_TT_SEASONAL = retime(Corner_TT_SEASONAL,'regular', @nanmean, 'TimeStep',dt);


 %% 
 
 % I'm now going to try the scatterplot with the ice-masked AODc & MAOD variables instead
 
 % Scatter plot, All Regions, daily data:
  % SUBPLOTTED
  
  %%

Corner_SST_SEASONALvals = [Corner_SST_SEASONAL.Corner_SST_daily; NaN; NaN];
Central_SST_SEASONALvals = [Central_SST_SEASONAL.Central_SST_daily;NaN;NaN];
Southernmost_SST_SEASONALvals = [Southernmost_SST_SEASONAL.Southernmost_SST_daily;NaN;NaN];
 

%%

clear Corner_TT_OD_Wind_daily Central_TT_OD_Wind_daily Southernmost_TT_OD_Wind_daily
clear x y nanVals p yfit R

x_Corner_Wind_SEASONAL_forAODc = Corner_TT_SEASONAL.MASTER_Winds_daytime;
y_Corner_AODc_SEASONAL = Corner_TT_SEASONAL.Corner_AODc_openocean;

nanVals_Corner_AODc_Wind = ismissing(x_Corner_Wind_SEASONAL_forAODc) | ismissing(y_Corner_AODc_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_SEASONAL_forAODc = x_Corner_Wind_SEASONAL_forAODc(~nanVals_Corner_AODc_Wind);
y_Corner_AODc_SEASONAL = y_Corner_AODc_SEASONAL(~nanVals_Corner_AODc_Wind);
Corner_SST_AODc = Corner_SST_SEASONALvals(~nanVals_Corner_AODc_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_SEASONAL_AODc_Wind = polyfit(x_Corner_Wind_SEASONAL_forAODc,y_Corner_AODc_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_SEASONAL_AODc_Wind = polyval(p_Corner_SEASONAL_AODc_Wind,x_Corner_Wind_SEASONAL_forAODc);

[R_Corner_SEASONAL_AODc_Wind, pval_Corner_SEASONAL_AODc_Wind] = corrcoef(x_Corner_Wind_SEASONAL_forAODc, y_Corner_AODc_SEASONAL );


% Central Seasonal scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_Wind_SEASONAL_forAODc = Central_TT_SEASONAL.MASTER_Winds_daytime;
y_Central_AODc_SEASONAL = Central_TT_SEASONAL.Central_AODc_openocean;

nanVals_Central_AODc_Wind = ismissing(x_Central_Wind_SEASONAL_forAODc) | ismissing(y_Central_AODc_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_SEASONAL_forAODc = x_Central_Wind_SEASONAL_forAODc(~nanVals_Central_AODc_Wind);
y_Central_AODc_SEASONAL = y_Central_AODc_SEASONAL(~nanVals_Central_AODc_Wind);
Central_SST_AODc = Central_SST_SEASONALvals(~nanVals_Central_AODc_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_SEASONAL_AODc_Wind = polyfit(x_Central_Wind_SEASONAL_forAODc,y_Central_AODc_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_SEASONAL_AODc_Wind = polyval(p_Central_SEASONAL_AODc_Wind,x_Central_Wind_SEASONAL_forAODc);

[R_Central_SEASONAL_AODc_Wind, pval_Central_SEASONAL_AODc_Wind] = corrcoef(x_Central_Wind_SEASONAL_forAODc, y_Central_AODc_SEASONAL );


% Southernmost Seasonal scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_Wind_SEASONAL_forAODc = Southernmost_TT_SEASONAL.MASTER_Winds_daytime;
y_Southernmost_AODc_SEASONAL = Southernmost_TT_SEASONAL.Southernmost_AODc_openocean;

nanVals_Southernmost_AODc_Wind = ismissing(x_Southernmost_Wind_SEASONAL_forAODc) | ismissing(y_Southernmost_AODc_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_SEASONAL_forAODc = x_Southernmost_Wind_SEASONAL_forAODc(~nanVals_Southernmost_AODc_Wind);
y_Southernmost_AODc_SEASONAL = y_Southernmost_AODc_SEASONAL(~nanVals_Southernmost_AODc_Wind);
Southernmost_SST_AODc = Southernmost_SST_SEASONALvals(~nanVals_Southernmost_AODc_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_SEASONAL_AODc_Wind = polyfit(x_Southernmost_Wind_SEASONAL_forAODc,y_Southernmost_AODc_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_SEASONAL_AODc_Wind = polyval(p_Southernmost_SEASONAL_AODc_Wind,x_Southernmost_Wind_SEASONAL_forAODc);

[R_Southernmost_SEASONAL_AODc_Wind, pval_Southernmost_SEASONAL_AODc_Wind] = corrcoef(x_Southernmost_Wind_SEASONAL_forAODc, y_Southernmost_AODc_SEASONAL );

%%
% AODt



clear Corner_TT_OD_Wind_daily Central_TT_OD_Wind_daily Southernmost_TT_OD_Wind_daily
clear x y nanVals p yfit R

x_Corner_Wind_SEASONAL_forAODt = Corner_TT_SEASONAL.MASTER_Winds_daytime;
y_Corner_AODt_SEASONAL = Corner_TT_SEASONAL.Corner_AOD_total_openocean;

nanVals_Corner_AODt_Wind = ismissing(x_Corner_Wind_SEASONAL_forAODt) | ismissing(y_Corner_AODt_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_SEASONAL_forAODt = x_Corner_Wind_SEASONAL_forAODt(~nanVals_Corner_AODt_Wind);
y_Corner_AODt_SEASONAL = y_Corner_AODt_SEASONAL(~nanVals_Corner_AODt_Wind);
Corner_SST_AODt = Corner_SST_SEASONALvals(~nanVals_Corner_AODt_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_SEASONAL_AODt_Wind = polyfit(x_Corner_Wind_SEASONAL_forAODt,y_Corner_AODt_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_SEASONAL_AODt_Wind = polyval(p_Corner_SEASONAL_AODt_Wind,x_Corner_Wind_SEASONAL_forAODt);

[R_Corner_SEASONAL_AODt_Wind, pval_Corner_SEASONAL_AODt_Wind] = corrcoef(x_Corner_Wind_SEASONAL_forAODt, y_Corner_AODt_SEASONAL );


% Central Seasonal scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_Wind_SEASONAL_forAODt = Central_TT_SEASONAL.MASTER_Winds_daytime;
y_Central_AODt_SEASONAL = Central_TT_SEASONAL.Central_AOD_total_openocean;

nanVals_Central_AODt_Wind = ismissing(x_Central_Wind_SEASONAL_forAODt) | ismissing(y_Central_AODt_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_SEASONAL_forAODt = x_Central_Wind_SEASONAL_forAODt(~nanVals_Central_AODt_Wind);
y_Central_AODt_SEASONAL = y_Central_AODt_SEASONAL(~nanVals_Central_AODt_Wind);
Central_SST_AODt = Central_SST_SEASONALvals(~nanVals_Central_AODt_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_SEASONAL_AODt_Wind = polyfit(x_Central_Wind_SEASONAL_forAODt,y_Central_AODt_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_SEASONAL_AODt_Wind = polyval(p_Central_SEASONAL_AODt_Wind,x_Central_Wind_SEASONAL_forAODt);

[R_Central_SEASONAL_AODt_Wind, pval_Central_SEASONAL_AODt_Wind] = corrcoef(x_Central_Wind_SEASONAL_forAODt, y_Central_AODt_SEASONAL );


% Southernmost Seasonal scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_Wind_SEASONAL_forAODt = Southernmost_TT_SEASONAL.MASTER_Winds_daytime;
y_Southernmost_AODt_SEASONAL = Southernmost_TT_SEASONAL.Southernmost_AOD_total_openocean;

nanVals_Southernmost_AODt_Wind = ismissing(x_Southernmost_Wind_SEASONAL_forAODt) | ismissing(y_Southernmost_AODt_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_SEASONAL_forAODt = x_Southernmost_Wind_SEASONAL_forAODt(~nanVals_Southernmost_AODt_Wind);
y_Southernmost_AODt_SEASONAL = y_Southernmost_AODt_SEASONAL(~nanVals_Southernmost_AODt_Wind);
Southernmost_SST_AODt = Southernmost_SST_SEASONALvals(~nanVals_Southernmost_AODt_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_SEASONAL_AODt_Wind = polyfit(x_Southernmost_Wind_SEASONAL_forAODt,y_Southernmost_AODt_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_SEASONAL_AODt_Wind = polyval(p_Southernmost_SEASONAL_AODt_Wind,x_Southernmost_Wind_SEASONAL_forAODt);

[R_Southernmost_SEASONAL_AODt_Wind, pval_Southernmost_SEASONAL_AODt_Wind] = corrcoef(x_Southernmost_Wind_SEASONAL_forAODt, y_Southernmost_AODt_SEASONAL );


%%
% AODf

clear Corner_TT_OD_Wind_daily Central_TT_OD_Wind_daily Southernmost_TT_OD_Wind_daily
clear x y nanVals p yfit R

x_Corner_Wind_SEASONAL_forAODf = Corner_TT_SEASONAL.MASTER_Winds_daytime;
y_Corner_AODf_SEASONAL = Corner_TT_SEASONAL.Corner_AOD_fine_openocean;

nanVals_Corner_AODf_Wind = ismissing(x_Corner_Wind_SEASONAL_forAODf) | ismissing(y_Corner_AODf_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_SEASONAL_forAODf = x_Corner_Wind_SEASONAL_forAODf(~nanVals_Corner_AODf_Wind);
y_Corner_AODf_SEASONAL = y_Corner_AODf_SEASONAL(~nanVals_Corner_AODf_Wind);
Corner_SST_AODf = Corner_SST_SEASONALvals(~nanVals_Corner_AODf_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_SEASONAL_AODf_Wind = polyfit(x_Corner_Wind_SEASONAL_forAODf,y_Corner_AODf_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_SEASONAL_AODf_Wind = polyval(p_Corner_SEASONAL_AODf_Wind,x_Corner_Wind_SEASONAL_forAODf);

[R_Corner_SEASONAL_AODf_Wind, pval_Corner_SEASONAL_AODf_Wind] = corrcoef(x_Corner_Wind_SEASONAL_forAODf, y_Corner_AODf_SEASONAL );


% Central Seasonal scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_Wind_SEASONAL_forAODf = Central_TT_SEASONAL.MASTER_Winds_daytime;
y_Central_AODf_SEASONAL = Central_TT_SEASONAL.Central_AOD_fine_openocean;

nanVals_Central_AODf_Wind = ismissing(x_Central_Wind_SEASONAL_forAODf) | ismissing(y_Central_AODf_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_SEASONAL_forAODf = x_Central_Wind_SEASONAL_forAODf(~nanVals_Central_AODf_Wind);
y_Central_AODf_SEASONAL = y_Central_AODf_SEASONAL(~nanVals_Central_AODf_Wind);
Central_SST_AODf = Central_SST_SEASONALvals(~nanVals_Central_AODf_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_SEASONAL_AODf_Wind = polyfit(x_Central_Wind_SEASONAL_forAODf,y_Central_AODf_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_SEASONAL_AODf_Wind = polyval(p_Central_SEASONAL_AODf_Wind,x_Central_Wind_SEASONAL_forAODf);

[R_Central_SEASONAL_AODf_Wind, pval_Central_SEASONAL_AODf_Wind] = corrcoef(x_Central_Wind_SEASONAL_forAODf, y_Central_AODf_SEASONAL );


% Southernmost Seasonal scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_Wind_SEASONAL_forAODf = Southernmost_TT_SEASONAL.MASTER_Winds_daytime;
y_Southernmost_AODf_SEASONAL = Southernmost_TT_SEASONAL.Southernmost_AOD_fine_openocean;

nanVals_Southernmost_AODf_Wind = ismissing(x_Southernmost_Wind_SEASONAL_forAODf) | ismissing(y_Southernmost_AODf_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_SEASONAL_forAODf = x_Southernmost_Wind_SEASONAL_forAODf(~nanVals_Southernmost_AODf_Wind);
y_Southernmost_AODf_SEASONAL = y_Southernmost_AODf_SEASONAL(~nanVals_Southernmost_AODf_Wind);
Southernmost_SST_AODf = Southernmost_SST_SEASONALvals(~nanVals_Southernmost_AODf_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_SEASONAL_AODf_Wind = polyfit(x_Southernmost_Wind_SEASONAL_forAODf,y_Southernmost_AODf_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_SEASONAL_AODf_Wind = polyval(p_Southernmost_SEASONAL_AODf_Wind,x_Southernmost_Wind_SEASONAL_forAODf);

[R_Southernmost_SEASONAL_AODf_Wind, pval_Southernmost_SEASONAL_AODf_Wind] = corrcoef(x_Southernmost_Wind_SEASONAL_forAODf, y_Southernmost_AODf_SEASONAL );

%%

clear x y nanVals p yfit R
x_Corner_Wind_SEASONAL_forMAOD = Corner_TT_SEASONAL.MASTER_Winds_nighttime;
y_Corner_MAOD_SEASONAL = Corner_TT_SEASONAL.MAOD_nocleanair_openocean;

nanVals_Corner_MAOD_Wind = ismissing(x_Corner_Wind_SEASONAL_forMAOD) | ismissing(y_Corner_MAOD_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_Wind_SEASONAL_forMAOD = x_Corner_Wind_SEASONAL_forMAOD(~nanVals_Corner_MAOD_Wind);
y_Corner_MAOD_SEASONAL = y_Corner_MAOD_SEASONAL(~nanVals_Corner_MAOD_Wind);
Corner_SST_MAOD = Corner_SST_SEASONALvals(~nanVals_Corner_MAOD_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_SEASONAL_MAOD_Wind = polyfit(x_Corner_Wind_SEASONAL_forMAOD,y_Corner_MAOD_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_SEASONAL_MAOD_Wind = polyval(p_Corner_SEASONAL_MAOD_Wind,x_Corner_Wind_SEASONAL_forMAOD);

[R_Corner_SEASONAL_MAOD_Wind, pval_Corner_SEASONAL_MAOD_Wind] = corrcoef(x_Corner_Wind_SEASONAL_forMAOD, y_Corner_MAOD_SEASONAL );


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_Wind_SEASONAL_forMAOD = Central_TT_SEASONAL.MASTER_Winds_nighttime;
y_Central_MAOD_SEASONAL = Central_TT_SEASONAL.MAOD_nocleanair_openocean;

nanVals_Central_MAOD_Wind = ismissing(x_Central_Wind_SEASONAL_forMAOD) | ismissing(y_Central_MAOD_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_Wind_SEASONAL_forMAOD = x_Central_Wind_SEASONAL_forMAOD(~nanVals_Central_MAOD_Wind);
y_Central_MAOD_SEASONAL = y_Central_MAOD_SEASONAL(~nanVals_Central_MAOD_Wind);
Central_SST_MAOD = Central_SST_SEASONALvals(~nanVals_Central_MAOD_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_SEASONAL_MAOD_Wind = polyfit(x_Central_Wind_SEASONAL_forMAOD,y_Central_MAOD_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_SEASONAL = polyval(p_Central_SEASONAL_MAOD_Wind,x_Central_Wind_SEASONAL_forMAOD);

[R_Central_SEASONAL_MAOD_Wind, pval_Central_SEASONAL_MAOD_Wind] = corrcoef(x_Central_Wind_SEASONAL_forMAOD, y_Central_MAOD_SEASONAL );


% Southernmost daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_Wind_SEASONAL_forMAOD = Southernmost_TT_SEASONAL.MASTER_Winds_nighttime;
y_Southernmost_MAOD_SEASONAL = Southernmost_TT_SEASONAL.MAOD_nocleanair_openocean;

nanVals_Southernmost_MAOD_Wind = ismissing(x_Southernmost_Wind_SEASONAL_forMAOD) | ismissing(y_Southernmost_MAOD_SEASONAL); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_Wind_SEASONAL_forMAOD = x_Southernmost_Wind_SEASONAL_forMAOD(~nanVals_Southernmost_MAOD_Wind);
y_Southernmost_MAOD_SEASONAL = y_Southernmost_MAOD_SEASONAL(~nanVals_Southernmost_MAOD_Wind);
Southernmost_SST_MAOD = Southernmost_SST_SEASONALvals(~nanVals_Southernmost_MAOD_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_SEASONAL_MAOD_Wind = polyfit(x_Southernmost_Wind_SEASONAL_forMAOD,y_Southernmost_MAOD_SEASONAL,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_SEASONAL_MAOD_Wind = polyval(p_Southernmost_SEASONAL_MAOD_Wind,x_Southernmost_Wind_SEASONAL_forMAOD);

[R_Southernmost_SEASONAL_MAOD_Wind, pval_Southernmost_SEASONAL_MAOD_Wind] = corrcoef(x_Southernmost_Wind_SEASONAL_forMAOD, y_Southernmost_MAOD_SEASONAL );


%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.08 0.05], [0.1 0.06]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

%Corner SEASONAL Scatterplot
subplot(4,3,1)
scatter(x_Corner_Wind_SEASONAL_forMAOD,y_Corner_MAOD_SEASONAL,40,Corner_SST_MAOD,...
    'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_SEASONAL_forMAOD,yfit_Corner_SEASONAL_MAOD_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca, 'FontSize', 15)

title('Region A') 

txt = ['y = ' num2str(round(p_Corner_SEASONAL_MAOD_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_SEASONAL_MAOD_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_SEASONAL_MAOD_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_SEASONAL_MAOD_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_SEASONAL_forMAOD),3,'significant'))];
Region = ('MAOD');

text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% colorbar
caxis([-1.5, 1.5])
cmocean('balance')

% Corner SEASONAL Scatterplot:
subplot(4,3,4)
scatter(x_Corner_Wind_SEASONAL_forAODc,y_Corner_AODc_SEASONAL,40,Corner_SST_AODc,...
'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_SEASONAL_forAODc,yfit_Corner_SEASONAL_AODc_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca,'FontSize', 15)


txt = ['y = ' num2str(round(p_Corner_SEASONAL_AODc_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_SEASONAL_AODc_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_SEASONAL_AODc_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_SEASONAL_AODc_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_SEASONAL_forAODc),3,'significant'))];
Region = ('AOD_C');

text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region

% colorbar
caxis([-1.5, 1.5])
cmocean('balance')

% Corner SEASONAL Scatterplot:
subplot(4,3,7)
scatter(x_Corner_Wind_SEASONAL_forAODt,y_Corner_AODt_SEASONAL,40,Corner_SST_AODt,...
'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_SEASONAL_forAODt,yfit_Corner_SEASONAL_AODt_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca,'FontSize', 15)


txt = ['y = ' num2str(round(p_Corner_SEASONAL_AODt_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_SEASONAL_AODt_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_SEASONAL_AODt_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_SEASONAL_AODt_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_SEASONAL_forAODt),3,'significant'))];
Region = ('AOD_T');

text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region

% colorbar
caxis([-1.5, 1.5])
cmocean('balance')

% Corner SEASONAL Scatterplot:
subplot(4,3,10)
scatter(x_Corner_Wind_SEASONAL_forAODf,y_Corner_AODf_SEASONAL,40,Corner_SST_AODf,...
'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_Wind_SEASONAL_forAODf,yfit_Corner_SEASONAL_AODf_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca,'FontSize', 15)


txt = ['y = ' num2str(round(p_Corner_SEASONAL_AODf_Wind(1),3,'significant')) 'x + ' num2str(round(p_Corner_SEASONAL_AODf_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_SEASONAL_AODf_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_SEASONAL_AODf_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_Wind_SEASONAL_forAODf),3,'significant'))];
Region = ('AOD_f');

text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region

% colorbar
caxis([-1.5, 1.5])
cmocean('balance')

% Central SEASONAL Scatterplot
subplot(4,3,2)
scatter(x_Central_Wind_SEASONAL_forMAOD,y_Central_MAOD_SEASONAL,40,Central_SST_MAOD,...
'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_SEASONAL_forMAOD,yfit_Central_SEASONAL,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca,'ytick',[], 'FontSize', 15)

title('Region B')

txt = ['y = ' num2str(round(p_Central_SEASONAL_MAOD_Wind(1),3,'significant')) 'x + ' num2str(round(p_Central_SEASONAL_MAOD_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_SEASONAL_MAOD_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_SEASONAL_MAOD_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_SEASONAL_forMAOD),3,'significant'))];
Region = ('MAOD');


text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region

% colorbar
caxis([-1.5, 1.5])
cmocean('balance')


% Central SEASONAL Scatterplot
subplot(4,3,5)
scatter(x_Central_Wind_SEASONAL_forAODc,y_Central_AODc_SEASONAL,40,Central_SST_AODc,...
    'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_SEASONAL_forAODc,yfit_Central_SEASONAL_AODc_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca, 'ytick',[], 'FontSize', 15)

txt = ['y = ' num2str(round(p_Central_SEASONAL_AODc_Wind(1),3,'significant')) 'x + ' num2str(round(p_Central_SEASONAL_AODf_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_SEASONAL_AODc_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_SEASONAL_AODc_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_SEASONAL_forAODc),3,'significant'))];
Region = ('AOD_C');

text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region
% colorbar
caxis([-1.5, 1.5])
cmocean('balance')



% Central SEASONAL Scatterplot
subplot(4,3,8)
scatter(x_Central_Wind_SEASONAL_forAODt,y_Central_AODt_SEASONAL,40,Central_SST_AODt,...
    'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_SEASONAL_forAODt,yfit_Central_SEASONAL_AODt_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca, 'ytick',[], 'FontSize', 15)



txt = ['y = ' num2str(round(p_Central_SEASONAL_AODt_Wind(1),3,'significant')) 'x + ' num2str(round(p_Central_SEASONAL_AODt_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_SEASONAL_AODt_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_SEASONAL_AODt_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_SEASONAL_forAODt),3,'significant'))];
Region = ('AOD_T');


text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region
% colorbar
caxis([-1.5, 1.5])
cmocean('balance')



% Central SEASONAL Scatterplot
subplot(4,3,11)
scatter(x_Central_Wind_SEASONAL_forAODf,y_Central_AODf_SEASONAL,40,Central_SST_AODf,...
    'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_Wind_SEASONAL_forAODf,yfit_Central_SEASONAL_AODf_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca, 'ytick',[], 'FontSize', 15)



txt = ['y = ' num2str(round(p_Central_SEASONAL_AODf_Wind(1),3,'significant')) 'x + ' num2str(round(p_Central_SEASONAL_AODf_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_SEASONAL_AODf_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_SEASONAL_AODf_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_Wind_SEASONAL_forAODf),3,'significant'))];
Region = ('AOD_f');


text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region
% colorbar
caxis([-1.5, 1.5])
cmocean('balance')

% Southernmost SEASONAL Scatterplot:
subplot(4,3,3)
scatter(x_Southernmost_Wind_SEASONAL_forMAOD,y_Southernmost_MAOD_SEASONAL, 40,Southernmost_SST_MAOD,...
    'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Southernmost_Wind_SEASONAL_forMAOD,yfit_Southernmost_SEASONAL_MAOD_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca,'ytick',[], 'FontSize', 15)
 
% xlabel(sprintf('SEASONAL Wind Speed average (m s^{-1})'));
title('Region C')
txt = ['y = ' num2str(round(p_Southernmost_SEASONAL_MAOD_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_SEASONAL_MAOD_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_SEASONAL_MAOD_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_SEASONAL_MAOD_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_SEASONAL_forMAOD),3,'significant'))];
Region = ('MAOD');

text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')
clear txt R_txt N_txt Region

% colorbar
caxis([-1.5, 1.5])
cmocean('balance')


% Southernmost SEASONAL Scatterplot:
subplot(4,3,6)
scatter(x_Southernmost_Wind_SEASONAL_forAODc,y_Southernmost_AODc_SEASONAL, 40,Southernmost_SST_AODc,...
    'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;

hold on
plot(x_Southernmost_Wind_SEASONAL_forAODc,yfit_Southernmost_SEASONAL_AODc_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca,'ytick',[],'FontSize', 15)
 
% xlabel(sprintf('SEASONAL Wind Speed average (m s^{-1})'));

txt = ['y = ' num2str(round(p_Southernmost_SEASONAL_AODc_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_SEASONAL_AODc_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_SEASONAL_AODc_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_SEASONAL_AODc_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_SEASONAL_forAODc),3,'significant'))];
Region = ('AOD_C');


text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region

% colorbar
caxis([-1.5, 1.5])
cmocean('balance')





% Southernmost SEASONAL Scatterplot:
subplot(4,3,9)
scatter(x_Southernmost_Wind_SEASONAL_forAODt,y_Southernmost_AODt_SEASONAL, 40,Southernmost_SST_AODt,...
    'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;

hold on
plot(x_Southernmost_Wind_SEASONAL_forAODt,yfit_Southernmost_SEASONAL_AODt_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca,'ytick',[],'FontSize', 15)
 
% xlabel(sprintf('SEASONAL Wind Speed average (m s^{-1})'));

txt = ['y = ' num2str(round(p_Southernmost_SEASONAL_AODt_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_SEASONAL_AODt_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_SEASONAL_AODt_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_SEASONAL_AODt_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_SEASONAL_forAODt),3,'significant'))];
Region = ('AOD_T');


text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region
% colorbar
caxis([-1.5,1.5])
cmocean('balance')


% Southernmost SEASONAL Scatterplot:
subplot(4,3,12)
scatter(x_Southernmost_Wind_SEASONAL_forAODf,y_Southernmost_AODf_SEASONAL, 40,Southernmost_SST_AODf,...
    'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;

hold on
plot(x_Southernmost_Wind_SEASONAL_forAODf,yfit_Southernmost_SEASONAL_AODf_Wind,'-',...
    'LineWidth', 2)

ylim([0 0.2])
xlim([3 16])
set(gca,'ytick',[],'FontSize', 15)
 
% xlabel(sprintf('SEASONAL Wind Speed average (m s^{-1})'));

txt = ['y = ' num2str(round(p_Southernmost_SEASONAL_AODf_Wind(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_SEASONAL_AODf_Wind(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_SEASONAL_AODf_Wind(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_SEASONAL_AODf_Wind(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_Wind_SEASONAL_forAODf),3,'significant'))];
Region = ('AOD_f');


text(4, 0.16, txt,'FontSize', 15)
text(4,0.145, R_txt,'FontSize', 15)
text(4,0.13,N_txt,'FontSize', 15)
text(4, 0.115, p_txt,'FontSize', 15)
text(3.4, 0.19, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% colorbar
caxis([-1.5,1.5])
cmocean('balance')


    
hp4 = get(subplot(4,3,9),'Position');
h = colorbar('Position', [0.95  0.25  0.02  0.6]);
h.FontWeight = 'bold';
h.FontSize = 15;


%%

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures/Diurnal_winds
set(gcf,'PaperPositionMode','auto')
print(gcf,'ALLoneFig_SST_MAOD_allAOD_allregions_SEASONAL_scatterplot.png','-dpng','-r400');       %  *// 400 dpi








