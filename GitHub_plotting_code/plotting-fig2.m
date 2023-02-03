% Plotting all variables in bellingshausen Sea across three different
% regions

%% Plots of three different regions in the Bellingshausen Sea
clearvars -except WINDS_TT_MINUTELY Total_timetable_SO_MOD_NEW Total_timetable_SO_Depol_Ratio_NEW
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Variables/Chlorophyll_Data_matfiles

load('Chlorophyll_regional_matfiles.mat')

% clearvars -except Central_aqua_daily_correct Corner_aqua_daily_correct Southernmost_aqua_daily_correct Total_timetable_SO_MOD_NEW Total_timetable_SO_Depol_Ratio_NEW 


% cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication
%%
% % Load windspeed timetable
% 
% load('Total_timetable_amsrmf.mat')

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Variables/New_WindSpeed_WindDirection_data/Srishti_AOD
load('WINDS_TT_MINUTELY.mat') % not the NoNAN variable


%% Pulling MAOD & Ice timetables.

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PhD_Phase_Two_SouthernOcean/Variables/2017_2020_vars

%
load('Total_timetable_SO_Depol_Ratio_NEW.mat') % this is with all data from 2020 (earlier it was excluding nov and dec months lol)
load('Total_timetable_SO_MOD_NEW.mat')% this is with all data from 2020 (earlier it was excluding nov and dec months lol)

% Redoing timeseries plots with lower pixel count threshold on AOD. This
% will be the most important figure, and i can go back and filter on a
% need-be basis.
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Variables/AOD_matfiles
load('AOD_timetables_2pixelthreshold.mat')

%%
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Figures/NEW_WINDS_Figs  

%%

% clearvars -except Total_timetable_SO_MOD_NEW Total_timetable_SO_Depol_Ratio_NEW Total_timetable_amsrmf
% Central_AOD_TT 
% Ice_TT = retime(Total_timetable_SO_Depol_Ratio_NEW, 'daily', @nanmean);
 Central_TT = addvars(Central_TT_BSea_daily, Central_aqua_daily_correct);
% writetimetable(Central_TT,'Central_TT.csv')

%% Bellingshausen Sea Region:

% This is the corner region without much chlorophyll change as seen in
% seasonal climatologies from chapter one. 

%  lat_BSea = [-60*ones(1,21) -63*ones(1,21) -60 -63]; % Outline of a box
%  lon_BSea = [-100:-80 -80:-1:-100 -100 -100];

% Wind:
Corner_BSea_Wind = WINDS_TT_MINUTELY.MASTER_Latitude >= -63.0...
    & WINDS_TT_MINUTELY.MASTER_Latitude <= -60.0...
    & WINDS_TT_MINUTELY.MASTER_Longitude >= -100.0...
    & WINDS_TT_MINUTELY.MASTER_Longitude <= -80.0; 
    
Corner_BSea_TT_Wind = WINDS_TT_MINUTELY(Corner_BSea_Wind,:);
Corner_BSea_TT_Wind = unique(Corner_BSea_TT_Wind);

% Ice:
Corner_BSea_Ice = Total_timetable_SO_Depol_Ratio_NEW.Total_Latitude_Ice >= -63.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Latitude_Ice <= -60.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Longitude_Ice >= -100.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Longitude_Ice <= -80.0; 
    
Corner_BSea_TT_Ice = Total_timetable_SO_Depol_Ratio_NEW(Corner_BSea_Ice,:);
Corner_BSea_TT_Ice = unique(Corner_BSea_TT_Ice);

% MOD: 
Corner_BSea_MOD = Total_timetable_SO_MOD_NEW.Total_Latitude_Surface >= -63.0...
    & Total_timetable_SO_MOD_NEW.Total_Latitude_Surface <= -60.0...
    & Total_timetable_SO_MOD_NEW.Total_Longitude_Surface >= -100.0...
    & Total_timetable_SO_MOD_NEW.Total_Longitude_Surface <= -80.0; 

Corner_BSea_TT_MOD = Total_timetable_SO_MOD_NEW(Corner_BSea_MOD,:);
Corner_BSea_TT_MOD = unique(Corner_BSea_TT_MOD);


% Pull out windspeed values less than 4 m /s and add it back in so that you are only using
% appropriate ratio values. 
% 
% Wind_speed_correct = Corner_BSea_TT_Wind.Total_windamsrMF;
% 
% Wind_speed_correct(Wind_speed_correct <4 ) = NaN;
% Corner_BSea_TT_Wind = addvars(Corner_BSea_TT_Wind, Wind_speed_correct);


% Pull out color ratio values and add it back in so that you are only using
% appropriate ratio values. 

Color_ratio_correct = Corner_BSea_TT_MOD.Color_Ratio;

Color_ratio_correct(Color_ratio_correct<0 | Color_ratio_correct>2) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
Corner_BSea_TT_MOD = addvars(Corner_BSea_TT_MOD, Color_ratio_correct);
 

% Remove outlier MAOD values
MAOD_filtered = Corner_BSea_TT_MOD.CMOD_Surface;
MAOD_filtered(MAOD_filtered<0 | MAOD_filtered > 0.5 | MAOD_filtered == 0) = NaN; 

MAOD_nocleanair = Corner_BSea_TT_MOD.CMOD_Surface;
MAOD_nocleanair(MAOD_nocleanair<0 |  MAOD_nocleanair == 0) = NaN; 
Corner_BSea_TT_MOD = addvars(Corner_BSea_TT_MOD, MAOD_filtered, MAOD_nocleanair);
 

% Remove outlier ice values
Ice_correct = Corner_BSea_TT_Ice.Total_Surface_532_Integrated_Depolarization_Ratio ;

Ice_correct(Ice_correct<0 | Ice_correct>2) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
Corner_BSea_TT_Ice = addvars(Corner_BSea_TT_Ice, Ice_correct);
 

% Daily
Corner_BSea_TT_Ice_daily = retime(Corner_BSea_TT_Ice, 'daily', @nanmean);
Corner_BSea_TT_MOD_daily = retime(Corner_BSea_TT_MOD, 'daily',@nanmean);
Corner_BSea_TT_amsrmf_daily = retime(Corner_BSea_TT_Wind, 'daily', @nanmean);

% Daily average and standard deviation: 
Corner_DAILY_AVG_MAOD = mean(Corner_BSea_TT_MOD_daily.MAOD_nocleanair, 'omitnan');
Corner_DAILY_AVG_STD_MAOD = std(Corner_BSea_TT_MOD_daily.MAOD_nocleanair, 'omitnan');

% 0.0989 +/- 0.0955

Corner_DAILY_AVG_AODc = mean(Corner_AOD_TT.Corner_AOD_coarse_daily, 'omitnan');
Corner_DAILY_AVG_STD_AODc = std(Corner_AOD_TT.Corner_AOD_coarse_daily, 'omitnan');
%0.0879 +/- 0.0780

% annually:

Corner_BSea_TT_Ice_yearly = retime(Corner_BSea_TT_Ice, 'yearly', @nanmean);
Corner_BSea_TT_MOD_yearly = retime(Corner_BSea_TT_MOD, 'yearly',@nanmean);
Corner_BSea_TT_amsrmf_yearly  = retime(Corner_BSea_TT_Wind, 'yearly', @nanmean);


% Annual average and standard deviation: 
Corner_ANNUAL_AVG_MAOD = mean(Corner_BSea_TT_MOD_yearly.MAOD_nocleanair, 'omitnan');
Corner_ANNUAL_AVG_STD_MAOD = std(Corner_BSea_TT_MOD_yearly.MAOD_nocleanair, 'omitnan');
%0.0864 +/- 0.0059


Corner_AOD_TT_yearly = retime(Corner_AOD_TT,'yearly', @nanmean);
Corner_ANNUAL_AVG_AODc = mean(Corner_AOD_TT_yearly.Corner_AOD_coarse_daily, 'omitnan');
Corner_ANNUAL_AVG_STD_AODc = std(Corner_AOD_TT_yearly.Corner_AOD_coarse_daily, 'omitnan');
%0.0880 +/- 0.0088

% Synchronizing these daily timetables so they span the same days:
SS = timerange('2007-01-01', '2018-12-31', 'closed'); %closed includes startime and endtime values

Corner_timetable_all_daily   = synchronize(Corner_BSea_TT_MOD_daily, Corner_BSea_TT_Ice_daily, Corner_BSea_TT_amsrmf_daily);
Corner_timetable_all_yearly  = synchronize(Corner_BSea_TT_MOD_yearly, Corner_BSea_TT_Ice_yearly, Corner_BSea_TT_amsrmf_yearly);

% And then create new timetables with all of these variables:
Corner_TT_BSea_daily   = Corner_timetable_all_daily(SS,:);
Corner_TT_BSea_yearly  = Corner_timetable_all_yearly(SS,:);


%%

% This is the central box region


%  lat_BSea_chl = [-66*ones(1,21) -69*ones(1,21) -66 -69]; % Outline of a box
%  lon_BSea_chl = [-87:-67 -67:-1:-87 -87 -87];
% 
%  m_line(lon_BSea_chl,lat_BSea_chl,'linewi',3,'color','k');     % Area outline ...
% 

% Wind:
Central_BSea_Wind = WINDS_TT_MINUTELY.MASTER_Latitude >= -69.0...
    & WINDS_TT_MINUTELY.MASTER_Latitude <= -66.0...
    & WINDS_TT_MINUTELY.MASTER_Longitude >= -87.0...
    & WINDS_TT_MINUTELY.MASTER_Longitude <= -67.0; 
    
Central_BSea_TT_Wind = WINDS_TT_MINUTELY(Central_BSea_Wind,:);
Central_BSea_TT_Wind = unique(Central_BSea_TT_Wind);

% Ice:
Central_BSea_Ice = Total_timetable_SO_Depol_Ratio_NEW.Total_Latitude_Ice >= -69.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Latitude_Ice <= -66.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Longitude_Ice >= -87.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Longitude_Ice <= -67.0; 
    
Central_BSea_TT_Ice = Total_timetable_SO_Depol_Ratio_NEW(Central_BSea_Ice,:);
Central_BSea_TT_Ice = unique(Central_BSea_TT_Ice);

% MOD: 
Central_BSea_MOD = Total_timetable_SO_MOD_NEW.Total_Latitude_Surface >= -69.0...
    & Total_timetable_SO_MOD_NEW.Total_Latitude_Surface <= -66.0...
    & Total_timetable_SO_MOD_NEW.Total_Longitude_Surface >= -87.0...
    & Total_timetable_SO_MOD_NEW.Total_Longitude_Surface <= -67.0; 

Central_BSea_TT_MOD = Total_timetable_SO_MOD_NEW(Central_BSea_MOD,:);
Central_BSea_TT_MOD = unique(Central_BSea_TT_MOD);



% Pull out windspeed values less than 4 m /s and add it back in so that you are only using
% appropriate ratio values. 
% 
% Wind_speed_correct = Central_BSea_TT_Wind.Total_windamsrMF;
% 
% Wind_speed_correct(Wind_speed_correct<4) = NaN;
% Central_BSea_TT_Wind = addvars(Central_BSea_TT_Wind, Wind_speed_correct);

% Pull out color ratio values and add it back in so that you are only using
% appropriate ratio values. 

Color_ratio_correct = Central_BSea_TT_MOD.Color_Ratio;

Color_ratio_correct(Color_ratio_correct<0 | Color_ratio_correct>2) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
Central_BSea_TT_MOD = addvars(Central_BSea_TT_MOD, Color_ratio_correct);
 

% Remove outlier MAOD values
MAOD_filtered = Central_BSea_TT_MOD.CMOD_Surface;
MAOD_filtered(MAOD_filtered<0 | MAOD_filtered > 0.5 | MAOD_filtered == 0) = NaN; 

MAOD_nocleanair = Central_BSea_TT_MOD.CMOD_Surface;
MAOD_nocleanair(MAOD_nocleanair<0 |  MAOD_nocleanair == 0) = NaN; 
Central_BSea_TT_MOD = addvars(Central_BSea_TT_MOD, MAOD_filtered, MAOD_nocleanair);
  

% Remove outlier ice values
Ice_correct = Central_BSea_TT_Ice.Total_Surface_532_Integrated_Depolarization_Ratio ;

Ice_correct(Ice_correct<0 | Ice_correct>2) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
Central_BSea_TT_Ice = addvars(Central_BSea_TT_Ice, Ice_correct);
 

% Daily,
Central_BSea_TT_Ice_daily = retime(Central_BSea_TT_Ice, 'daily', @nanmean);
Central_BSea_TT_MOD_daily = retime(Central_BSea_TT_MOD, 'daily',@nanmean);
Central_BSea_TT_amsrmf_daily = retime(Central_BSea_TT_Wind, 'daily', @nanmean);

% Daily average and standard deviation: 
Central_DAILY_AVG_MAOD = mean(Central_BSea_TT_MOD_daily.MAOD_nocleanair, 'omitnan');
Central_DAILY_AVG_STD_MAOD = std(Central_BSea_TT_MOD_daily.MAOD_nocleanair, 'omitnan');

% 0.0638 +/- 0.0721

Central_DAILY_AVG_AODc = mean(Central_AOD_TT.Central_AOD_coarse_daily, 'omitnan');
Central_DAILY_AVG_STD_AODc = std(Central_AOD_TT.Central_AOD_coarse_daily, 'omitnan');
% 0.0848 +/- 0.0744


% annually:
Central_BSea_TT_Ice_yearly = retime(Central_BSea_TT_Ice, 'yearly', @nanmean);
Central_BSea_TT_MOD_yearly = retime(Central_BSea_TT_MOD, 'yearly',@nanmean);
Central_BSea_TT_amsrmf_yearly  = retime(Central_BSea_TT_Wind, 'yearly', @nanmean);

% Annual average and standard deviation: 
Central_ANNUAL_AVG_MAOD = mean(Central_BSea_TT_MOD_yearly.MAOD_nocleanair, 'omitnan');
Central_ANNUAL_AVG_STD_MAOD = std(Central_BSea_TT_MOD_yearly.MAOD_nocleanair, 'omitnan');
%0.0526 +/- 0.0086


Central_AOD_TT_yearly = retime(Central_AOD_TT,'yearly', @nanmean);
Central_ANNUAL_AVG_AODc = mean(Central_AOD_TT_yearly.Central_AOD_coarse_daily, 'omitnan');
Central_ANNUAL_AVG_STD_AODc = std(Central_AOD_TT_yearly.Central_AOD_coarse_daily, 'omitnan');
%0.0850 +/- 0.0107



% Synchronizing these daily timetables so they span the same days:
SS = timerange('2007-01-01', '2018-12-31', 'closed'); %closed includes startime and endtime values

Central_timetable_all_daily   = synchronize(Central_BSea_TT_MOD_daily, Central_BSea_TT_Ice_daily, Central_BSea_TT_amsrmf_daily);
Central_timetable_all_yearly  = synchronize(Central_BSea_TT_MOD_yearly, Central_BSea_TT_Ice_yearly, Central_BSea_TT_amsrmf_yearly);

% And then create new timetables with all of these variables:
Central_TT_BSea_daily   = Central_timetable_all_daily(SS,:);
Central_TT_BSea_yearly  = Central_timetable_all_yearly(SS,:);

%%


% % Open ocean: 
% openocean = Central_TT_BSea_daily.Ice_correct ;
% openocean(openocean > 0.15) = NaN;
% Central_TT_BSea_daily = addvars(Central_TT_BSea_daily, openocean);
% 
% 
% MAOD_nocleanair_openocean = Central_TT_BSea_daily.MAOD_nocleanair;
% MAOD_nocleanair_openocean(isnan(Central_TT_BSea_daily.openocean)) = NaN;
% Central_TT_BSea_daily = addvars(Central_TT_BSea_daily, MAOD_nocleanair_openocean);

%%

% this is the southernmost box

%  lat_BSea_extra = [-73*ones(1,21) -70*ones(1,21) -70 -73]; % Outline of a box
%  lon_BSea_extra = [-95:-75 -75:-1:-95 -95 -95];
%  m_line(lon_BSea_extra,lat_BSea_extra,'linewi',3,'color','magenta');     % Area outline ...


% Wind:
Southernmost_BSea_Wind = WINDS_TT_MINUTELY.MASTER_Latitude >= -73.0...
    & WINDS_TT_MINUTELY.MASTER_Latitude <= -70.0...
    & WINDS_TT_MINUTELY.MASTER_Longitude >= -95.0...
    & WINDS_TT_MINUTELY.MASTER_Longitude <= -75.0; 
    
Southernmost_BSea_TT_Wind = WINDS_TT_MINUTELY(Southernmost_BSea_Wind,:);
Southernmost_BSea_TT_Wind = unique(Southernmost_BSea_TT_Wind);

% Ice:
Southernmost_BSea_Ice = Total_timetable_SO_Depol_Ratio_NEW.Total_Latitude_Ice >= -73.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Latitude_Ice <= -70.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Longitude_Ice >= -95.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Longitude_Ice <= -75.0; 
    
Southernmost_BSea_TT_Ice = Total_timetable_SO_Depol_Ratio_NEW(Southernmost_BSea_Ice,:);
Southernmost_BSea_TT_Ice = unique(Southernmost_BSea_TT_Ice);

% MOD: 
Southernmost_BSea_MOD = Total_timetable_SO_MOD_NEW.Total_Latitude_Surface >= -73.0...
    & Total_timetable_SO_MOD_NEW.Total_Latitude_Surface <= -70.0...
    & Total_timetable_SO_MOD_NEW.Total_Longitude_Surface >= -95.0...
    & Total_timetable_SO_MOD_NEW.Total_Longitude_Surface <= -75.0; 

Southernmost_BSea_TT_MOD = Total_timetable_SO_MOD_NEW(Southernmost_BSea_MOD,:);
Southernmost_BSea_TT_MOD = unique(Southernmost_BSea_TT_MOD);


% Pull out windspeed values less than 4 m /s and add it back in so that you are only using
% appropriate ratio values. 
% 
% Wind_speed_correct = Southernmost_BSea_TT_Wind.Total_windamsrMF;
% 
% Wind_speed_correct(Wind_speed_correct<4) = NaN;
% Southernmost_BSea_TT_Wind = addvars(Southernmost_BSea_TT_Wind, Wind_speed_correct);


% Pull out color ratio values and add it back in so that you are only using
% appropriate ratio values. 

Color_ratio_correct = Southernmost_BSea_TT_MOD.Color_Ratio;

Color_ratio_correct(Color_ratio_correct<0 | Color_ratio_correct>2) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
Southernmost_BSea_TT_MOD = addvars(Southernmost_BSea_TT_MOD, Color_ratio_correct);
 

% Remove outlier MAOD values
MAOD_filtered = Southernmost_BSea_TT_MOD.CMOD_Surface;
MAOD_filtered(MAOD_filtered<0 | MAOD_filtered > 0.5 | MAOD_filtered == 0) = NaN; 

MAOD_nocleanair = Southernmost_BSea_TT_MOD.CMOD_Surface;
MAOD_nocleanair(MAOD_nocleanair<0 |  MAOD_nocleanair == 0) = NaN; 
Southernmost_BSea_TT_MOD = addvars(Southernmost_BSea_TT_MOD, MAOD_filtered, MAOD_nocleanair);
  

% Remove outlier ice values
Ice_correct = Southernmost_BSea_TT_Ice.Total_Surface_532_Integrated_Depolarization_Ratio ;

Ice_correct(Ice_correct<0 | Ice_correct>2) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
Southernmost_BSea_TT_Ice = addvars(Southernmost_BSea_TT_Ice, Ice_correct);
 


% Daily, 
Southernmost_BSea_TT_Ice_daily = retime(Southernmost_BSea_TT_Ice, 'daily', @nanmean);
Southernmost_BSea_TT_MOD_daily = retime(Southernmost_BSea_TT_MOD, 'daily',@nanmean);
Southernmost_BSea_TT_amsrmf_daily = retime(Southernmost_BSea_TT_Wind, 'daily', @nanmean);

% Daily average and standard deviation: 
Southernmost_DAILY_AVG_MAOD = mean(Southernmost_BSea_TT_MOD_daily.MAOD_nocleanair, 'omitnan');
Southernmost_DAILY_AVG_STD_MAOD = std(Southernmost_BSea_TT_MOD_daily.MAOD_nocleanair, 'omitnan');

%0.0461 +/- 0.0545

Southernmost_DAILY_AVG_AODc = mean(Southernmost_AOD_TT.Southernmost_AOD_coarse_daily, 'omitnan');
Southernmost_DAILY_AVG_STD_AODc = std(Southernmost_AOD_TT.Southernmost_AOD_coarse_daily, 'omitnan');

%0.0680 +/- 0.0422

% annually:

Southernmost_BSea_TT_Ice_yearly = retime(Southernmost_BSea_TT_Ice, 'yearly', @nanmean);
Southernmost_BSea_TT_MOD_yearly = retime(Southernmost_BSea_TT_MOD, 'yearly',@nanmean);
Southernmost_BSea_TT_amsrmf_yearly  = retime(Southernmost_BSea_TT_Wind, 'yearly', @nanmean);


% Annual average and standard deviation: 
Southernmost_ANNUAL_AVG_MAOD = mean(Southernmost_BSea_TT_MOD_yearly.MAOD_nocleanair, 'omitnan');
Southernmost_ANNUAL_AVG_STD_MAOD = std(Southernmost_BSea_TT_MOD_yearly.MAOD_nocleanair, 'omitnan');

%0.0342 +/- 0.0048

Southernmost_AOD_TT_yearly = retime(Southernmost_AOD_TT,'yearly', @nanmean);
Southernmost_ANNUAL_AVG_AODc = mean(Southernmost_AOD_TT_yearly.Southernmost_AOD_coarse_daily, 'omitnan');
Southernmost_ANNUAL_AVG_STD_AODc = std(Southernmost_AOD_TT_yearly.Southernmost_AOD_coarse_daily, 'omitnan');

%0.0682 +/- 0.0074

% Synchronizing these daily timetables so they span the same days:
SS = timerange('2007-01-01', '2018-12-31', 'closed'); %closed includes startime and endtime values

Southernmost_timetable_all_daily   = synchronize(Southernmost_BSea_TT_MOD_daily, Southernmost_BSea_TT_Ice_daily, Southernmost_BSea_TT_amsrmf_daily);
Southernmost_timetable_all_yearly  = synchronize(Southernmost_BSea_TT_MOD_yearly, Southernmost_BSea_TT_Ice_yearly, Southernmost_BSea_TT_amsrmf_yearly);

% And then create new timetables with all of these variables:
Southernmost_TT_BSea_daily  = Southernmost_timetable_all_daily(SS,:);
Southernmost_TT_BSea_yearly = Southernmost_timetable_all_yearly(SS,:);

%%
times_daily =  Corner_TT_BSea_daily.Total_Profile_Time_New_Surface;
% times_8day = Southernmost_TT_BSea_8day.Total_Profile_Time_New_Surface;
% times_monthly = Corner_TT_BSea_monthly.Total_Profile_Time_New_Surface;
%%

% Average annual mean:



%%
 %%%%%% ----- SUBPLOTTED? ------ %%%%%%%%%
 % doc addaxis If needed
 
 apricot = rgb('apricot');
 orange = rgb('orange');
 red = rgb('scarlet');
 gray = rgb('grey');
 black = rgb('black');
 grass_green = rgb('true green');
 ice_blue = rgb('periwinkle blue');
 royal_blue = rgb('royal blue');
 pink = rgb('light pink');
 olive = rgb('deep brown');
 
 %%
 
                         %%%%% DAILY TIMESERIES, CORNER SPATIAL REGIONS %%%%%%

%                          
%         MAOD_plus_std = plus(Corner_TT_BSea_daily.MAOD_correct, Corner_TT_BSea_daily.MAOD_std);  
%         MAOD_minus_std = minus(Corner_TT_BSea_daily.MAOD_correct,Corner_TT_BSea_daily.MAOD_std); 

        %%
        
        % moving mean effort:
        
        MAOD_Corner_movmean = movmean(Corner_TT_BSea_daily.MAOD_nocleanair, 30,'omitnan');
        AODc_Corner_movmean = movmean(Corner_AOD_TT.Corner_AOD_coarse_daily, 30, 'omitnan');
        CHL_Corner_movmean = movmean(Corner_aqua_daily_correct, 30, 'omitnan');
%         Color_Corner_movmean = movmean(Corner_TT_BSea_daily.Color_ratio_correct,30, 'omitnan');
        Ice_Corner_movmean = movmean(Corner_TT_BSea_daily.Ice_correct, 30,'omitnan');        
        Wind_Corner_movmean = movmean(Corner_TT_BSea_daily.MASTER_Winds, 30,'omitnan');
        
        
        MAOD_Central_movmean = movmean(Central_TT_BSea_daily.MAOD_nocleanair, 30,'omitnan');
        AODc_Central_movmean = movmean(Central_AOD_TT.Central_AOD_coarse_daily, 30, 'omitnan');
        CHL_Central_movmean = movmean(Central_aqua_daily_correct, 30, 'omitnan');
%         Color_Central_movmean = movmean(Central_TT_BSea_daily.Color_ratio_correct,30, 'omitnan');
        Ice_Central_movmean = movmean(Central_TT_BSea_daily.Ice_correct, 30,'omitnan');
        Wind_Central_movmean = movmean(Central_TT_BSea_daily.MASTER_Winds, 30,'omitnan');
        
        
        
        MAOD_Southernmost_movmean = movmean(Southernmost_TT_BSea_daily.MAOD_nocleanair, 30,'omitnan');
        AODc_Southernmost_movmean = movmean(Southernmost_AOD_TT.Southernmost_AOD_coarse_daily, 30,'omitnan');
        CHL_Southernmost_movmean = movmean(Southernmost_aqua_daily_correct, 30, 'omitnan');
%         Color_Southernmost_movmean = movmean(Southernmost_TT_BSea_daily.Color_ratio_correct,30, 'omitnan');
        Ice_Southernmost_movmean = movmean(Southernmost_TT_BSea_daily.Ice_correct, 30,'omitnan');
        Wind_Southernmost_movmean = movmean(Southernmost_TT_BSea_daily.MASTER_Winds, 30,'omitnan');
        
        
 %%
 
 make_it_tight = true;
 subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.02], [0.1 0.025], [0.1 0.05]);
 if ~make_it_tight,  clear subplot;  end
 x = times_daily;
 
 figure; clf;
 
 % MAOD
 ax(1) = subplot(5,3,1);
 aa_splot(x, Corner_TT_BSea_daily.MAOD_nocleanair, '-',...
     'linewidth', 0.5, ...
     'Color', gray);
 
 hold on 

 aa_splot(x, MAOD_Corner_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', black);
 
 
 xlim([x(1) x(end)]) 
 ylim([0 0.3])
 
 ylabel('MAOD')
 title('Region A')
% l(1)= legend('CMOD', 'Position',[0.76 0.92 0.1 0.0349]);

legend('MAOD')
%  legend('Position', [0.105 0.94 0.16 0.0349])
 
  grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
%  xtickformat('MM-yyyy')
 
 
 % MAOD
 ax(2) = subplot(5,3,2);
 aa_splot(x, Central_TT_BSea_daily.MAOD_nocleanair, '-',...
     'linewidth', 0.5, ...
     'Color', gray);
 
 hold on 

 aa_splot(x, MAOD_Central_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', black);
 
 
 xlim([x(1) x(end)]) 
 ylim([0 0.3])
 
 ylabel('MAOD')

legend('MAOD')
title('Region B')
%  legend('Position', [0.105 0.94 0.16 0.0349])
 
  grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
 
 
 
 % MAOD
 ax(3) = subplot(5,3,3);
 aa_splot(x, Southernmost_TT_BSea_daily.MAOD_nocleanair, '-',...
     'linewidth', 0.5, ...
     'Color', gray);
 
 hold on 

 aa_splot(x, MAOD_Southernmost_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', black);
 
 
 xlim([x(1) x(end)]) 
 ylim([0 0.3])
 
 ylabel('MAOD')
title('Region C')
legend('MAOD')
%  legend('Position', [0.105 0.94 0.16 0.0349])
 
  grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
 
 
 
 
 
 % AODc
 ax(4) = subplot(5,3,4);
 aa_splot(x, Corner_AOD_TT.Corner_AOD_coarse_daily, '-',...
     'linewidth', 0.5, ...
     'Color', gray);
 hold on 

 aa_splot(x, AODc_Corner_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', black);
 
 xlim([x(1) x(end)]) 
 ylim([0 0.3])
 
 ylabel('AOD_c')

legend('AOD_c')
%  legend('Position', [0.105 0.94 0.16 0.0349])
 
  grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
%  xtickformat('MM-yyyy')
 
 
 
 % AODc
 ax(5) = subplot(5,3,5);
 aa_splot(x, Central_AOD_TT.Central_AOD_coarse_daily, '-',...
     'linewidth', 0.5, ...
     'Color', gray);
 hold on 

 aa_splot(x, AODc_Central_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', black);
 
 xlim([x(1) x(end)]) 
 ylim([0 0.3])
 
 ylabel('AOD_c')
% l(1)= legend('CMOD', 'Position',[0.76 0.92 0.1 0.0349]);

legend('AOD_c')
%  legend('Position', [0.105 0.94 0.16 0.0349])
 
  grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
%  xtickformat('MM-yyyy')
 
 
 
 
 % AODc
 ax(6) = subplot(5,3,6);
 aa_splot(x, Southernmost_AOD_TT.Southernmost_AOD_coarse_daily, '-',...
     'linewidth', 0.5, ...
     'Color', gray);
 hold on 

 aa_splot(x, AODc_Southernmost_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', black);
 
 xlim([x(1) x(end)]) 
 ylim([0 0.3])
 
 ylabel('AOD_c')
% l(1)= legend('CMOD', 'Position',[0.76 0.92 0.1 0.0349]);

legend('AOD_c')
%  legend('Position', [0.105 0.94 0.16 0.0349])
 
  grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
%  xtickformat('MM-yyyy')
 
 
 
 % CHL
 ax(7) = subplot(5,3,7);
 aa_splot(x, Corner_aqua_daily_correct, '-',...
     'linewidth', 0.5, ...
     'Color', grass_green);
 hold on 

 aa_splot(x, CHL_Corner_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', olive);
 
 xlim([x(1) x(end)]) 
 ylim([0 5])
 
 ylabel('mg m^{-3}')
% l(1)= legend('CMOD', 'Position',[0.76 0.92 0.1 0.0349]);

legend('chl-\ita')
%  legend('Position', [0.105 0.94 0.16 0.0349])
 
  grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
%  xtickformat('MM-yyyy')
 
 
 AX=findall(0,'type','axes');
 set(AX, 'FontSize', 10)

 
 
 % CHL
 ax(8) = subplot(5,3,8);
 aa_splot(x, Central_aqua_daily_correct, '-',...
     'linewidth', 0.5, ...
     'Color', grass_green);
 hold on 

 aa_splot(x, CHL_Central_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', olive);
 
 xlim([x(1) x(end)]) 
 ylim([0 5])
 
 ylabel('mg m^{-3}')
% l(1)= legend('CMOD', 'Position',[0.76 0.92 0.1 0.0349]);

legend('chl-\ita')
%  legend('Position', [0.105 0.94 0.16 0.0349])
 
  grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
%  xtickformat('MM-yyyy')
 
 
 % CHL
 ax(9) = subplot(5,3,9);
 aa_splot(x, Southernmost_aqua_daily_correct, '-',...
     'linewidth', 0.5, ...
     'Color', grass_green);
 hold on 

 aa_splot(x, CHL_Southernmost_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', olive);
 
 xlim([x(1) x(end)]) 
 ylim([0 5])
 
 ylabel('mg m^{-3}')
% l(1)= legend('CMOD', 'Position',[0.76 0.92 0.1 0.0349]);

legend('chl-\ita')
%  legend('Position', [0.105 0.94 0.16 0.0349])
 
  grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
%  xtickformat('MM-yyyy')
 
 
 

 % ICE
 ax(10) = subplot(5,3,10);
 aa_splot(x,  Corner_TT_BSea_daily.Ice_correct,'-',...
     'linewidth', 0.5, ...
     'MarkerSize', 4,...
     'MarkerFaceColor', ice_blue,...
     'MarkerEdgeColor', ice_blue,...
     'Color', ice_blue);
  xlim([x(1)  x(end)]) 
  
hold on 

 aa_splot(x, Ice_Corner_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', royal_blue);
 
   ylim([0 1]) 
 
 ylabel('sea ice',...
     'FontName','Helvetica Neue');%         xlabel('Months')
 legend( 'ice')
 
 
 
 grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
%  xtickformat('MM-yyyy')
 
 
 
 % ICE
 ax(11) = subplot(5,3,11);
 aa_splot(x,  Central_TT_BSea_daily.Ice_correct,'-',...
     'linewidth', 0.5, ...
     'MarkerSize', 4,...
     'MarkerFaceColor', ice_blue,...
     'MarkerEdgeColor', ice_blue,...
     'Color', ice_blue);
  xlim([x(1)  x(end)]) 
  
hold on 

 aa_splot(x, Ice_Central_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', royal_blue);
 
   ylim([0 1]) 
 
 ylabel('\delta',...
     'FontName','Helvetica Neue');%         xlabel('Months')
 legend( 'ice')
 
 
 
 grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
%  xtickformat('MM-yyyy')
 
 
 % ICE
 ax(12) = subplot(5,3,12);
 aa_splot(x,  Southernmost_TT_BSea_daily.Ice_correct,'-',...
     'linewidth', 0.5, ...
     'MarkerSize', 4,...
     'MarkerFaceColor', ice_blue,...
     'MarkerEdgeColor', ice_blue,...
     'Color', ice_blue);
  xlim([x(1)  x(end)]) 
  
hold on 

 aa_splot(x, Ice_Southernmost_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', royal_blue);
 
   ylim([0 1]) 
 
 ylabel('\delta',...
     'FontName','Helvetica Neue');%         xlabel('Months')
 legend( 'ice')
 
 
 
 grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
%  xtickformat('MM-yyyy')
 
 
 
 
 %
 % WIND
 %
 ax(13) = subplot(5,3,13);
 aa_splot(x, Corner_TT_BSea_daily.MASTER_Winds, '-',...
     'linewidth', 0.5, ...
     'MarkerSize', 4,...
     'MarkerFaceColor', pink,...
     'MarkerEdgeColor', pink,...
     'Color', pink);
 hold on 

 aa_splot(x, Wind_Corner_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', red);
 
 ylabel('m s^{-1}');
 legend('wind')
 
 xlim([x(1)  x(end)]) 
   ylim([0 20])
 
 grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
 xtickformat('MM-yyyy')
 
 
 
 xtickangle(90)


 
 %
 % WIND
 %
 ax(14) = subplot(5,3,14);
 aa_splot(x, Central_TT_BSea_daily.MASTER_Winds, '-',...
     'linewidth', 0.5, ...
     'MarkerSize', 4,...
     'MarkerFaceColor', pink,...
     'MarkerEdgeColor', pink,...
     'Color', pink);
 hold on 

 aa_splot(x, Wind_Central_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', red);
 
 ylabel('m s^{-1}');
 legend('wind')
 
  xtickangle(90)

 xlim([x(1)  x(end)]) 
   ylim([0 20])
 
 grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
 xtickformat('MM-yyyy')
 

%
 % WIND
 %
 ax(15) = subplot(5,3,15);
 aa_splot(x, Southernmost_TT_BSea_daily.MASTER_Winds, '-',...
     'linewidth', 0.5, ...
     'MarkerSize', 4,...
     'MarkerFaceColor', pink,...
     'MarkerEdgeColor', pink,...
     'Color', pink);
 hold on 

 aa_splot(x, Wind_Southernmost_movmean, '-',...
     'linewidth', 1.5, ...
     'Color', red);
 
 ylabel('m s^{-1}');
 legend('wind')
 
 xlim([x(1)  x(end)]) 
   ylim([0 20])
 
 grid on
 set(gca, 'XTick', (times_daily(1) : calmonths(6) : times_daily(end)) );
 xtickformat('MM-yyyy')
 
 
 AX=findall(0,'type','axes');
 set(AX, 'FontSize', 14)
 
 xtickangle(90)
[ax(1:12).XTickLabel] = deal([]);
a = [2,3,5,6,8,9,11,12,14,15];
[ax(a).YTickLabel] = deal([]);
[ax(a).YLabel] = deal([]);



%%
set(gcf,'PaperPositionMode','auto')
print(gcf,'V2_NEWWINDS_ALLVARS_ALLREGIONS_TimeSeries_daily','-dpng','-r300');       %  *// 300 dpi

