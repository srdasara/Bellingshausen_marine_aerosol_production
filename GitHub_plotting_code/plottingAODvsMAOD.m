%% Plots of three different regions in the Bellingshausen Sea

% I'm comparing coarse-mode AOD and MAOD for Lynn
%%


cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_vars

% Scatterplot Code:
clear;
load('Southernmost_complete_TT_daily.mat')
load('Corner_complete_TT_daily.mat')
load('Central_complete_TT_daily.mat')

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures/Diurnal_winds

%%

clear x y nanVals p yfit R
x_Corner_MAOD_Daily = Corner_complete_TT_daily.MAOD_nocleanair_openocean;
y_Corner_AODc_Daily = Corner_complete_TT_daily.Corner_AOD_fine_openocean;

nanVals_Corner_AODc_Wind = ismissing(x_Corner_MAOD_Daily) | ismissing(y_Corner_AODc_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_MAOD_Daily = x_Corner_MAOD_Daily(~nanVals_Corner_AODc_Wind);
y_Corner_AODc_Daily = y_Corner_AODc_Daily(~nanVals_Corner_AODc_Wind);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_Daily_MAOD_AODc = polyfit(x_Corner_MAOD_Daily,y_Corner_AODc_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_Daily_MAOD_AODc = polyval(p_Corner_Daily_MAOD_AODc,x_Corner_MAOD_Daily);

[R_Corner_Daily_MAOD_AODc, pval_Corner_Daily_MAOD_AODc] = corrcoef(x_Corner_MAOD_Daily, y_Corner_AODc_Daily );


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_MAOD_Daily = Central_complete_TT_daily.MAOD_nocleanair_openocean;
y_Central_AODc_Daily = Central_complete_TT_daily.Central_AOD_fine_openocean;

nanVals_Central_MAOD_AODc = ismissing(x_Central_MAOD_Daily) | ismissing(y_Central_AODc_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_MAOD_Daily = x_Central_MAOD_Daily(~nanVals_Central_MAOD_AODc);
y_Central_AODc_Daily = y_Central_AODc_Daily(~nanVals_Central_MAOD_AODc);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_Daily_MAOD_AODc = polyfit(x_Central_MAOD_Daily,y_Central_AODc_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_Daily_MAOD_AODc = polyval(p_Central_Daily_MAOD_AODc,x_Central_MAOD_Daily);

[R_Central_Daily_MAOD_AODc, pval_Central_Daily_MAOD_AODc] = corrcoef(x_Central_MAOD_Daily, y_Central_AODc_Daily );


% Southernmost daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_MAOD_Daily = Southernmost_complete_TT_daily.MAOD_nocleanair_openocean;
y_Southernmost_AODc_Daily = Southernmost_complete_TT_daily.Southernmost_AOD_fine_openocean;

nanVals_Southernmost_MAOD_AODc = ismissing(x_Southernmost_MAOD_Daily) | ismissing(y_Southernmost_AODc_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_MAOD_Daily = x_Southernmost_MAOD_Daily(~nanVals_Southernmost_MAOD_AODc);
y_Southernmost_AODc_Daily = y_Southernmost_AODc_Daily(~nanVals_Southernmost_MAOD_AODc);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_Daily_MAOD_AODc = polyfit(x_Southernmost_MAOD_Daily,y_Southernmost_AODc_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_Daily_MAOD_AODc = polyval(p_Southernmost_Daily_MAOD_AODc,x_Southernmost_MAOD_Daily);

[R_Southernmost_Daily_MAOD_AODc, pval_Southernmost_Daily_MAOD_AODc] = corrcoef(x_Southernmost_MAOD_Daily, y_Southernmost_AODc_Daily );


%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

%Corner Daily Scatterplot
subplot(3,3,1)
scatter(x_Corner_MAOD_Daily,y_Corner_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_MAOD_Daily,yfit_Corner_Daily_MAOD_AODc,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 0.4])
set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)

title('(a) Region A') 

txt = ['y = ' num2str(round(p_Corner_Daily_MAOD_AODc(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_MAOD_AODc(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_MAOD_AODc(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_MAOD_AODc(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_MAOD_Daily),3,'significant'))];
Region = ('AOD_C');

text(0.15, 0.4, txt,'FontSize', 15)
text(0.15,0.35, R_txt,'FontSize', 15)
text(0.15,0.3,N_txt,'FontSize', 15)
text(0.15,0.25,p_txt, 'FontSize', 15)
text(0.01, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region



% Central Daily Scatterplot
subplot(3,3,2)
scatter(x_Central_MAOD_Daily,y_Central_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_MAOD_Daily,yfit_Central_Daily_MAOD_AODc,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 0.4])
% ylabel('Marine Aerosol Optical Depth (MAOD)', 'FontSize', 18);
set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)
set(gca,'ytick',[],'yticklabel',[], 'FontSize', 15)

title('(b) Region B')

txt = ['y = ' num2str(round(p_Central_Daily_MAOD_AODc(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily_MAOD_AODc(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_MAOD_AODc(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_MAOD_AODc(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_MAOD_Daily),3,'significant'))];
Region = ('AOD_C');

text(0.15, 0.4, txt,'FontSize', 15)
text(0.15,0.35, R_txt,'FontSize', 15)
text(0.15,0.3,N_txt,'FontSize', 15)
text(0.15,0.25,p_txt, 'FontSize', 15)
text(0.01, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% Southernmost Daily Scatterplot:
subplot(3,3,3)
scatter(x_Southernmost_MAOD_Daily,y_Southernmost_AODc_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Southernmost_MAOD_Daily,yfit_Southernmost_Daily_MAOD_AODc,'-',...
    'LineWidth', 2)


ylim([0 0.5])
xlim([0 0.4])
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[], 'FontSize', 15)

% xlabel(sprintf('coarse-mode Aerosol Optical Depth (AOD_C)'),'FontSize', 18);
title('(c) Region C')
txt = ['y = ' num2str(round(p_Southernmost_Daily_MAOD_AODc(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_MAOD_AODc(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_MAOD_AODc(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_MAOD_AODc(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_MAOD_Daily),3,'significant'))];
Region = ('AOD_C');

text(0.15, 0.4, txt,'FontSize', 15)
text(0.15,0.35, R_txt,'FontSize', 15)
text(0.15,0.3,N_txt,'FontSize', 15)
text(0.15,0.25,p_txt, 'FontSize', 15)
text(0.01, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
clear txt R_txt N_txt Region

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
if ~make_it_tight,  clear subplot;  end

%%
hold on;

%Corner Daily Scatterplot
subplot(3,3,4)
scatter(x_Corner_MAOD_Daily,y_Corner_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_MAOD_Daily,yfit_Corner_Daily_MAOD_AODc,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 0.4])
set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)

% title('(a) Region A') 

txt = ['y = ' num2str(round(p_Corner_Daily_MAOD_AODc(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_MAOD_AODc(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_MAOD_AODc(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_MAOD_AODc(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_MAOD_Daily),3,'significant'))];
Region = ('AOD_T');
ylabel('Marine Aerosol Optical Depth (MAOD)', 'FontSize', 18);

text(0.15, 0.4, txt,'FontSize', 15)
text(0.15,0.35, R_txt,'FontSize', 15)
text(0.15,0.3,N_txt,'FontSize', 15)
text(0.15,0.25,p_txt, 'FontSize', 15)
text(0.01, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region



% Central Daily Scatterplot
subplot(3,3,5)
scatter(x_Central_MAOD_Daily,y_Central_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_MAOD_Daily,yfit_Central_Daily_MAOD_AODc,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 0.4])
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[], 'FontSize', 15)

% title('(b) Region B')

txt = ['y = ' num2str(round(p_Central_Daily_MAOD_AODc(1),3,'significant')) 'x + ' num2str(p_Central_Daily_MAOD_AODc(2))];
R_txt = ['R = ' num2str(R_Central_Daily_MAOD_AODc(1,2))];
p_txt = ['p = ' num2str(pval_Central_Daily_MAOD_AODc(1,2))];
N_txt = ['n = ' num2str(length(x_Central_MAOD_Daily))];
Region = ('AOD_T');

text(0.15, 0.4, txt,'FontSize', 15)
text(0.15,0.35, R_txt,'FontSize', 15)
text(0.15,0.3,N_txt,'FontSize', 15)
text(0.15,0.25,p_txt, 'FontSize', 15)
text(0.01, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% Southernmost Daily Scatterplot:
subplot(3,3,6)
scatter(x_Southernmost_MAOD_Daily,y_Southernmost_AODc_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Southernmost_MAOD_Daily,yfit_Southernmost_Daily_MAOD_AODc,'-',...
    'LineWidth', 2)

set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[], 'FontSize', 15)
ylim([0 0.5])
xlim([0 0.4])
set(gca,'FontSize', 15)
 
% xlabel(sprintf('coarse-mode Aerosol Optical Depth (AOD_C)'),'FontSize', 18);
% title('(c) Region C')
txt = ['y = ' num2str(p_Southernmost_Daily_MAOD_AODc(1)) 'x + ' num2str(p_Southernmost_Daily_MAOD_AODc(2))];
R_txt = ['R = ' num2str(R_Southernmost_Daily_MAOD_AODc(1,2))];
p_txt = ['p = ' num2str(pval_Southernmost_Daily_MAOD_AODc(1,2))];
N_txt = ['n = ' num2str(length(x_Southernmost_MAOD_Daily))];
Region = ('AOD_T');

text(0.15, 0.4, txt,'FontSize', 15)
text(0.15,0.35, R_txt,'FontSize', 15)
text(0.15,0.3,N_txt,'FontSize', 15)
text(0.15,0.25,p_txt, 'FontSize', 15)
text(0.01, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
clear txt R_txt N_txt Region


%%
hold on;

%Corner Daily Scatterplot
subplot(3,3,7)
scatter(x_Corner_MAOD_Daily,y_Corner_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_MAOD_Daily,yfit_Corner_Daily_MAOD_AODc,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 0.4])
set(gca,'FontSize', 15)

% title('(a) Region A') 

txt = ['y = ' num2str(p_Corner_Daily_MAOD_AODc(1)) 'x + ' num2str(p_Corner_Daily_MAOD_AODc(2))];
R_txt = ['R = ' num2str(R_Corner_Daily_MAOD_AODc(1,2))];
p_txt = ['p = ' num2str(pval_Corner_Daily_MAOD_AODc(1,2))];
N_txt = ['n = ' num2str(length(x_Corner_MAOD_Daily))];
Region = ('AOD_f');

text(0.15, 0.4, txt,'FontSize', 15)
text(0.15,0.35, R_txt,'FontSize', 15)
text(0.15,0.3,N_txt,'FontSize', 15)
text(0.15,0.25,p_txt, 'FontSize', 15)
text(0.01, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
% xlabel(sprintf('AOD_C'), 'FontSize', 18);

clear txt R_txt N_txt Region



% Central Daily Scatterplot
subplot(3,3,8)
scatter(x_Central_MAOD_Daily,y_Central_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_MAOD_Daily,yfit_Central_Daily_MAOD_AODc,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 0.4])
% ylabel('Marine Aerosol Optical Depth (MAOD)', 'FontSize', 18);
set(gca,'ytick',[],'yticklabel',[], 'FontSize', 15)

% title('(b) Region B')
xlabel(sprintf('Aerosol Optical Depth (AOD)'),'FontSize', 18);

txt = ['y = ' num2str(p_Central_Daily_MAOD_AODc(1)) 'x + ' num2str(p_Central_Daily_MAOD_AODc(2))];
R_txt = ['R = ' num2str(R_Central_Daily_MAOD_AODc(1,2))];
p_txt = ['p = ' num2str(pval_Central_Daily_MAOD_AODc(1,2))];
N_txt = ['n = ' num2str(length(x_Central_MAOD_Daily))];
Region = ('AOD_f');

text(0.15, 0.4, txt,'FontSize', 15)
text(0.15,0.35, R_txt,'FontSize', 15)
text(0.15,0.3,N_txt,'FontSize', 15)
text(0.15,0.25,p_txt, 'FontSize', 15)
text(0.01, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

% text(0.7, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% Southernmost Daily Scatterplot:
subplot(3,3,9)
scatter(x_Southernmost_MAOD_Daily,y_Southernmost_AODc_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Southernmost_MAOD_Daily,yfit_Southernmost_Daily_MAOD_AODc,'-',...
    'LineWidth', 2)

set(gca,'ytick',[],'yticklabel',[], 'FontSize', 15)

ylim([0 0.5])
xlim([0 0.4])
set(gca,'FontSize', 15)
 
% title('(c) Region C')
txt = ['y = ' num2str(p_Southernmost_Daily_MAOD_AODc(1)) 'x + ' num2str(p_Southernmost_Daily_MAOD_AODc(2))];
R_txt = ['R = ' num2str(R_Southernmost_Daily_MAOD_AODc(1,2))];
p_txt = ['p = ' num2str(pval_Southernmost_Daily_MAOD_AODc(1,2))];
N_txt = ['n = ' num2str(length(x_Southernmost_MAOD_Daily))];
Region = ('AOD_f');

text(0.15, 0.4, txt,'FontSize', 15)
text(0.15,0.35, R_txt,'FontSize', 15)
text(0.15,0.3,N_txt,'FontSize', 15)
text(0.15,0.25,p_txt, 'FontSize', 15)
text(0.01, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
clear txt R_txt N_txt Region

%%

set(gcf,'PaperPositionMode','auto')
print(gcf,'MAOD_vs_allAOD.png','-dpng','-r300');       %  *// 300 dpi


