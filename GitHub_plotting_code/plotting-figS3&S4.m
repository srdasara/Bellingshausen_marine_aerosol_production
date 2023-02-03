
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_vars

load('Corner_complete_TT_daily.mat')

%%
figure;clf;

subplot(2,1,1)
Wind_allseasons_nighttime = histogram(Corner_complete_TT_daily.MASTER_Winds_nighttime,'EdgeAlpha',0.25);
Wind_allseasons_nighttime.Normalization = 'probability';
m_winds_night = median(Corner_complete_TT_daily.MASTER_Winds_nighttime,'omitnan');
xline([m_winds_night],'-',{'Median = 10.04 m/s'}, 'Color', 'k', 'FontSize', 20,'LineWidth', 1)
title('Nighttime')
set(gca, 'FontSize', 16);
% getting the count
N_nighttime = length(Corner_complete_TT_daily.MASTER_Winds_nighttime(~isnan(Corner_complete_TT_daily.MASTER_Winds_nighttime)));
N_nighttime = ['N for nighttime winds = ' num2str(N_nighttime)];
text(17, 0.09, N_nighttime, 'FontSize', 16)

subplot(2,1,2)
Wind_allseasons_daytime = histogram(Corner_complete_TT_daily.MASTER_Winds_daytime,'EdgeAlpha',0.25);
Wind_allseasons_daytime.Normalization = 'probability';
m_winds_day = median(Corner_complete_TT_daily.MASTER_Winds_daytime,'omitnan');
xline([m_winds_day],'-',{'Median = 9.98 m/s'}, 'Color', 'k', 'FontSize', 20,'LineWidth', 1)
title('Daytime')
set(gca, 'FontSize', 16);
% getting the count
N_daytime = length(Corner_complete_TT_daily.MASTER_Winds_daytime(~isnan(Corner_complete_TT_daily.MASTER_Winds_daytime)));
N_daytime = ['N for daytime winds = ' num2str(N_daytime)];
text(17, 0.09, N_daytime, 'FontSize', 16)

%%
% cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures/Diurnal_winds
% set(gcf,'PaperPositionMode','auto')
% print(gcf,'Diurnal_WindSpeedDistribution_Histogram','-dpng','-r300');       %  *// 300 dpi


%%
% High all seasons winds in corner spatial region

% Daytime is only for AOD, nighttime is only for MAOD
High_winds_allseasons_index_nighttime = Corner_complete_TT_daily.MASTER_Winds_nighttime >= 10.04;
Low_winds_allseasons_index_nighttime = Corner_complete_TT_daily.MASTER_Winds_nighttime < 10.04;

High_winds_allseasons_index_daytime = Corner_complete_TT_daily.MASTER_Winds_daytime >= 9.98;
Low_winds_allseasons_index_daytime = Corner_complete_TT_daily.MASTER_Winds_daytime < 9.98;


High_Winds_allseasons_Corner_TT_nighttime = Corner_complete_TT_daily(High_winds_allseasons_index_nighttime,:);
Low_Winds_allseasons_Corner_TT_nighttime = Corner_complete_TT_daily(Low_winds_allseasons_index_nighttime,:);


High_Winds_allseasons_Corner_TT_daytime = Corner_complete_TT_daily(High_winds_allseasons_index_daytime,:);
Low_Winds_allseasons_Corner_TT_daytime = Corner_complete_TT_daily(Low_winds_allseasons_index_daytime,:);


%%

% Anderson Darling test to examine whether populations are normal: 
% If 0, fails to reject the null hypothesis at the default 5% significance
% level. (i.e., population is normal)
% If 1, population is not normal
% 
% adtest(High_Winds_Winter_Corner_daily_TT.Corner_AOD_coarse_daily) % not normal
% adtest(Low_Winds_Winter_Corner_daily_TT.Corner_AOD_coarse_daily) % yes normal

adtest(High_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean) % not normal
adtest(Low_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean) % not normal
adtest(High_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean)
adtest(Low_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean)

% Wilcoxon Rank Sum Test

% high wind vs low wind
% [p_Ww_winter, h_Ww_winter] = ranksum(High_Winds_Winter_Corner_daily_TT.Corner_AOD_coarse_daily,...
%     Low_Winds_Winter_Corner_daily_TT.Corner_AOD_coarse_daily); % significant

[p_Ww_allseasons_AODc, h_Ww_allseasons_AODc] = ranksum(High_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean,...
    Low_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean);

[p_Ww_allseasons_MAOD, h_Ww_allseasons_MAOD] = ranksum(High_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean,...
    Low_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean);


%% Plotting of histogram:

make_it_tight = false;
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.05], [0.08 0.1], [0.1 0.03]);
if ~make_it_tight,  clear subplot;  end
 
% 

Num_HighWindsCorner_allseasons_MAOD = length(High_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean(~isnan(High_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean)));
Num_LowWindsCorner_allseasons_MAOD = length(Low_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean(~isnan(Low_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean)));

Num_HighWindsCorner_allseasons_AODc = length(High_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean(~isnan(High_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean)));
Num_LowWindsCorner_allseasons_AODc = length(Low_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean(~isnan(Low_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean)));

figure;clf;

subplot(2,1,1)
High_Wind_allseasons = histogram(High_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean,100);
hold on
Low_Wind_allseasons = histogram(Low_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean,50);
High_Wind_allseasons.Normalization = 'probability';
Low_Wind_allseasons.Normalization = 'probability';
% subtitle('MAOD')
legend('Wind Speeds >= 10.04 m s^{-1}', 'Wind Speeds < 10.04 m s^{-1}');

hold on;
m_highwinds = median(High_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean,'omitnan');
m_lowwinds = median(Low_Winds_allseasons_Corner_TT_nighttime.MAOD_nocleanair_openocean,'omitnan');
% xline([m_highwinds m_lowwinds],'-',{'High Winds median','Low Winds median'}, 'Color', 'b', 'FontSize', 12,'LineWidth', 1)
xline([m_highwinds m_lowwinds],'-', 'Color', 'b', 'FontSize', 12,'LineWidth', 1)

legend('Wind Speeds >= 10.04 m s^{-1}', 'Wind Speeds < 10.04 m s^{-1}');
ylim([0 0.15])
xlim([0 0.5])
ax.XTick = [];

median_highwinds = ['median OD for winds >= 10.04 m s^{-1} = ' num2str(m_highwinds)];
median_lowwinds = ['median OD for winds < 10.04 m s^{-1} = ' num2str(m_lowwinds)];


P_txt = ['p = ' num2str(p_Ww_allseasons_MAOD)];
N_txt_highwinds = ['N for winds >= 10.04 m s^{-1} = ' num2str(Num_HighWindsCorner_allseasons_MAOD)];
N_txt_lowwinds = ['N for winds < 10.04 m s^{-1} = ' num2str(Num_LowWindsCorner_allseasons_MAOD)];
text(0.2, 0.09, median_highwinds, 'FontSize',16)
text(0.2, 0.075, median_lowwinds, 'FontSize',16)

text(0.2 ,0.06, P_txt,'FontSize', 16)
text(0.2, 0.045, N_txt_highwinds, 'FontSize', 16)
text(0.2 , 0.03, N_txt_lowwinds, 'FontSize', 16)
set(gca, 'FontSize', 16);

type_txt = 'MAOD (nighttime wind speed)';
text(0.15 , 0.125, type_txt, 'FontSize', 18, 'FontWeight', 'bold')

clear type_txt P_txt N_txt_highwinds N_txt_lowwinds m_highwinds m_lowwinds

subplot(2,1,2)
High_Wind_allseasons = histogram(High_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean,100);
hold on
Low_Wind_allseasons = histogram(Low_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean,70);
High_Wind_allseasons.Normalization = 'probability';
Low_Wind_allseasons.Normalization = 'probability';
% subtitle('AOD_c')

clear m_highwinds m_lowwinds

hold on;
m_highwinds = median(High_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean, 'omitnan');
m_lowwinds = median(Low_Winds_allseasons_Corner_TT_daytime.Corner_AODc_openocean, 'omitnan');
% xline([m_highwinds m_lowwinds],'-',{'High Winds median','Low Winds median'}, 'Color', 'b', 'FontSize', 12,'LineWidth', 1)
xline([m_highwinds m_lowwinds],'-', 'Color', 'b', 'FontSize', 12,'LineWidth', 1)

legend('Wind Speeds >= 9.98 m s^{-1}', 'Wind Speeds < 9.98 m s^{-1}');
ylim([0 0.15])
xlim([0 0.5])
ax.XTick = [];


P_txt = ['p = ' num2str(p_Ww_allseasons_AODc)];
N_txt_highwinds = ['N for winds >= 9.98 m s^{-1} = ' num2str(Num_HighWindsCorner_allseasons_AODc)];
N_txt_lowwinds = ['N for winds < 9.98 m s^{-1} = ' num2str(Num_LowWindsCorner_allseasons_AODc)];
median_highwinds = ['median OD for winds >= 9.98 m s^{-1} = ' num2str(m_highwinds)];
median_lowwinds = ['median OD for winds < 9.98 m s^{-1} = ' num2str(m_lowwinds)];

text(0.2, 0.09, median_highwinds, 'FontSize',16)
text(0.2, 0.075, median_lowwinds, 'FontSize',16)
text(0.2 ,0.06, P_txt,'FontSize', 16)
text(0.2, 0.045, N_txt_highwinds, 'FontSize', 16)
text(0.2 , 0.03, N_txt_lowwinds, 'FontSize', 16)

set(gca, 'FontSize', 16);

type_txt = 'AOD_C (daytime wind speed)';
text(0.15 , 0.125, type_txt, 'FontSize', 18, 'FontWeight', 'bold')


sgt = sgtitle('Examining Optical Depth (OD) in high and low wind speed conditions in Region A');
sgt.FontSize = 19;
sgt.FontWeight = 'bold';


%%
set(gcf,'PaperPositionMode','auto')
print(gcf,'DIURNALWINDS_OPENOCEAN_Histogram_MAOD_CoarseAOD_HighWindvsLowWind_6ms_CornerRegion','-dpng','-r300');       %  *// 300 dpi
