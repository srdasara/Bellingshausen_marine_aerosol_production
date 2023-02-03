%%


cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Variables/Chlorophyll_Data_matfiles

load('Chlorophyll_regional_matfiles.mat')



%%


% cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication
% 
% % Load windspeed timetable
% 
% load('Total_timetable_amsrmf.mat')

% cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Variables/New_WindSpeed_WindDirection_data/Srishti_AOD
% load('WINDS_TT_MINUTELY.mat') % not the NoNAN variable

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Variables/AOD_matfiles
load('AOD_timetables_2pixelthreshold.mat')

% Pulling MAOD timetables.

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PhD_Phase_Two_SouthernOcean/Variables/2017_2020_vars

load('Total_timetable_SO_MOD_NEW.mat')% this is with all data from 2020 (earlier it was excluding nov and dec months lol)
load('Total_timetable_SO_Depol_Ratio_NEW.mat') % this is with all data from 2020 (earlier it was excluding nov and dec months lol)

%% Folder that figures are saved in: 
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures




%% Bellingshausen Sea Region:

% This is the corner region without much chlorophyll change as seen in
% seasonal climatologies from chapter one. 

%  lat_BSea = [-60*ones(1,21) -63*ones(1,21) -60 -63]; % Outline of a box
%  lon_BSea = [-100:-80 -80:-1:-100 -100 -100];

% Wind:
% Corner_BSea_Wind = WINDS_TT_MINUTELY.MASTER_Latitude >= -63.0...
%     & WINDS_TT_MINUTELY.MASTER_Latitude <= -60.0...
%     & WINDS_TT_MINUTELY.MASTER_Longitude >= -100.0...
%     & WINDS_TT_MINUTELY.MASTER_Longitude <= -80.0; 
%     
% Corner_BSea_TT_Wind = WINDS_TT_MINUTELY(Corner_BSea_Wind,:);
% Corner_BSea_TT_Wind = unique(Corner_BSea_TT_Wind);

% % Ice:
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

% Wind_speed_correct = Corner_BSea_TT_Wind.Total_windamsrMF;
% 
% Wind_speed_correct(Wind_speed_correct <4 ) = NaN;
% Corner_BSea_TT_Wind = addvars(Corner_BSea_TT_Wind, Wind_speed_correct);


% Pull out color ratio values and add it back in so that you are only using
% appropriate ratio values. 

% Color_ratio_correct = Corner_BSea_TT_MOD.Color_Ratio;
% 
% Color_ratio_correct(Color_ratio_correct<0 | Color_ratio_correct>2) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
% Corner_BSea_TT_MOD = addvars(Corner_BSea_TT_MOD, Color_ratio_correct);
%  

% Remove outlier MAOD values
MAOD_filtered = Corner_BSea_TT_MOD.CMOD_Surface;
MAOD_filtered(MAOD_filtered<0 | MAOD_filtered > 0.5 | MAOD_filtered == 0) = NaN; 

MAOD_nocleanair = Corner_BSea_TT_MOD.CMOD_Surface;
MAOD_nocleanair(MAOD_nocleanair<0 |  MAOD_nocleanair == 0) = NaN; 
Corner_BSea_TT_MOD = addvars(Corner_BSea_TT_MOD, MAOD_filtered, MAOD_nocleanair);
 

% Remove outlier ice values
Ice_correct = Corner_BSea_TT_Ice.Total_Surface_532_Integrated_Depolarization_Ratio ;

Ice_correct(Ice_correct<0 | Ice_correct>5) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
Corner_BSea_TT_Ice = addvars(Corner_BSea_TT_Ice, Ice_correct);



% Daily
Corner_BSea_TT_Ice_daily = retime(Corner_BSea_TT_Ice, 'daily', @nanmean);
Corner_BSea_TT_MOD_daily = retime(Corner_BSea_TT_MOD, 'daily',@nanmean);
% Corner_BSea_TT_amsrmf_daily = retime(Corner_BSea_TT_Wind, 'daily', @nanmean);

% annually:

% Corner_BSea_TT_Ice_yearly = retime(Corner_BSea_TT_Ice, 'yearly', @nanmean);
% Corner_BSea_TT_MOD_yearly = retime(Corner_BSea_TT_MOD, 'yearly',@nanmean);
% Corner_BSea_TT_amsrmf_yearly  = retime(Corner_BSea_TT_Wind, 'yearly', @nanmean);



% Synchronizing these daily timetables so they span the same days:
SS = timerange('2007-01-01', '2018-12-31', 'closed'); %closed includes startime and endtime values

Corner_timetable_all_daily   = synchronize(Corner_BSea_TT_MOD_daily, Corner_BSea_TT_Ice_daily);
% 
% And then create new timetables with all of these variables:
Corner_TT_BSea_daily  = Corner_timetable_all_daily(SS,:);


% Open ocean: 
openocean = Corner_TT_BSea_daily.Ice_correct ;
openocean(openocean > 0.15) = NaN;
Corner_TT_BSea_daily = addvars(Corner_TT_BSea_daily, openocean);


MAOD_nocleanair_openocean = Corner_TT_BSea_daily.MAOD_nocleanair;
MAOD_nocleanair_openocean(isnan(Corner_TT_BSea_daily.openocean)) = NaN;
Corner_TT_BSea_daily = addvars(Corner_TT_BSea_daily, MAOD_nocleanair_openocean);


%%

% This is the central box region


%  lat_BSea_chl = [-66*ones(1,21) -69*ones(1,21) -66 -69]; % Outline of a box
%  lon_BSea_chl = [-87:-67 -67:-1:-87 -87 -87];
% 
%  m_line(lon_BSea_chl,lat_BSea_chl,'linewi',3,'color','k');     % Area outline ...
% 

% Wind:
% Central_BSea_Wind = WINDS_TT_MINUTELY.MASTER_Latitude >= -69.0...
%     & WINDS_TT_MINUTELY.MASTER_Latitude <= -66.0...
%     & WINDS_TT_MINUTELY.MASTER_Longitude >= -87.0...
%     & WINDS_TT_MINUTELY.MASTER_Longitude <= -67.0; 
%     
% Central_BSea_TT_Wind = WINDS_TT_MINUTELY(Central_BSea_Wind,:);
% Central_BSea_TT_Wind = unique(Central_BSea_TT_Wind);

% Ice:
Central_BSea_Ice = Total_timetable_SO_Depol_Ratio_NEW.Total_Latitude_Ice >= -69.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Latitude_Ice <= -66.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Longitude_Ice >= -87.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Longitude_Ice <= -67.0; 
    
Central_BSea_TT_Ice = Total_timetable_SO_Depol_Ratio_NEW(Central_BSea_Ice,:);
Central_BSea_TT_Ice = unique(Central_BSea_TT_Ice);
% 
% % MOD: 
Central_BSea_MOD = Total_timetable_SO_MOD_NEW.Total_Latitude_Surface >= -69.0...
    & Total_timetable_SO_MOD_NEW.Total_Latitude_Surface <= -66.0...
    & Total_timetable_SO_MOD_NEW.Total_Longitude_Surface >= -87.0...
    & Total_timetable_SO_MOD_NEW.Total_Longitude_Surface <= -67.0; 

Central_BSea_TT_MOD = Total_timetable_SO_MOD_NEW(Central_BSea_MOD,:);
Central_BSea_TT_MOD = unique(Central_BSea_TT_MOD);



% Pull out windspeed values less than 4 m /s and add it back in so that you are only using
% % appropriate ratio values. 
% 
% Wind_speed_correct = Central_BSea_TT_Wind.Total_windamsrMF;
% 
% Wind_speed_correct(Wind_speed_correct<4) = NaN;
% Central_BSea_TT_Wind = addvars(Central_BSea_TT_Wind, Wind_speed_correct);
% 
% % Pull out color ratio values and add it back in so that you are only using
% appropriate ratio values. 

% Color_ratio_correct = Central_BSea_TT_MOD.Color_Ratio;
% 
% Color_ratio_correct(Color_ratio_correct<0 | Color_ratio_correct>2) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
% Central_BSea_TT_MOD = addvars(Central_BSea_TT_MOD, Color_ratio_correct);
%  
% 
% % Remove outlier MAOD values
MAOD_filtered = Central_BSea_TT_MOD.CMOD_Surface;
MAOD_filtered(MAOD_filtered<0 | MAOD_filtered > 0.5 | MAOD_filtered == 0) = NaN; 

MAOD_nocleanair = Central_BSea_TT_MOD.CMOD_Surface;
MAOD_nocleanair(MAOD_nocleanair<0 |  MAOD_nocleanair == 0) = NaN; 
Central_BSea_TT_MOD = addvars(Central_BSea_TT_MOD, MAOD_filtered, MAOD_nocleanair);
%   
% 
% % Remove outlier ice values
Ice_correct = Central_BSea_TT_Ice.Total_Surface_532_Integrated_Depolarization_Ratio ;

Ice_correct(Ice_correct<0 | Ice_correct>5) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
Central_BSea_TT_Ice = addvars(Central_BSea_TT_Ice, Ice_correct);


% % Daily,
Central_BSea_TT_Ice_daily = retime(Central_BSea_TT_Ice, 'daily', @nanmean);
Central_BSea_TT_MOD_daily = retime(Central_BSea_TT_MOD, 'daily',@nanmean);
% Central_BSea_TT_amsrmf_daily = retime(Central_BSea_TT_Wind, 'daily', @nanmean);


% annually:
% Central_BSea_TT_Ice_yearly = retime(Central_BSea_TT_Ice, 'yearly', @nanmean);
% Central_BSea_TT_MOD_yearly = retime(Central_BSea_TT_MOD, 'yearly',@nanmean);
% Central_BSea_TT_amsrmf_yearly  = retime(Central_BSea_TT_Wind, 'yearly', @nanmean);


% Synchronizing these daily timetables so they span the same days:
SS = timerange('2007-01-01', '2018-12-31', 'closed'); %closed includes startime and endtime values

Central_timetable_all_daily   = synchronize(Central_BSea_TT_MOD_daily,Central_BSea_TT_Ice_daily);
% Central_timetable_all_yearly  = synchronize(Central_BSea_TT_MOD_yearly, Central_BSea_TT_Ice_yearly, Central_BSea_TT_amsrmf_yearly);

% And then create new timetables with all of these variables:
Central_TT_BSea_daily   = Central_timetable_all_daily(SS,:);


% Open ocean: 
openocean = Central_TT_BSea_daily.Ice_correct ;
openocean(openocean > 0.15) = NaN;
Central_TT_BSea_daily = addvars(Central_TT_BSea_daily, openocean);


MAOD_nocleanair_openocean = Central_TT_BSea_daily.MAOD_nocleanair;
MAOD_nocleanair_openocean(isnan(Central_TT_BSea_daily.openocean)) = NaN;
Central_TT_BSea_daily = addvars(Central_TT_BSea_daily, MAOD_nocleanair_openocean);


%%

% this is the southernmost box

%  lat_BSea_extra = [-73*ones(1,21) -70*ones(1,21) -70 -73]; % Outline of a box
%  lon_BSea_extra = [-95:-75 -75:-1:-95 -95 -95];
%  m_line(lon_BSea_extra,lat_BSea_extra,'linewi',3,'color','magenta');     % Area outline ...


% Wind:
% Southernmost_BSea_Wind = WINDS_TT_MINUTELY.MASTER_Latitude >= -73.0...
%     & WINDS_TT_MINUTELY.MASTER_Latitude <= -70.0...
%     & WINDS_TT_MINUTELY.MASTER_Longitude >= -95.0...
%     & WINDS_TT_MINUTELY.MASTER_Longitude <= -75.0; 
    
% Southernmost_BSea_TT_Wind = WINDS_TT_MINUTELY(Southernmost_BSea_Wind,:);
% Southernmost_BSea_TT_Wind = unique(Southernmost_BSea_TT_Wind);

% Ice:
Southernmost_BSea_Ice = Total_timetable_SO_Depol_Ratio_NEW.Total_Latitude_Ice >= -73.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Latitude_Ice <= -70.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Longitude_Ice >= -95.0...
    & Total_timetable_SO_Depol_Ratio_NEW.Total_Longitude_Ice <= -75.0; 
    
Southernmost_BSea_TT_Ice = Total_timetable_SO_Depol_Ratio_NEW(Southernmost_BSea_Ice,:);
Southernmost_BSea_TT_Ice = unique(Southernmost_BSea_TT_Ice);

% % MOD: 
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
% 
% Color_ratio_correct = Southernmost_BSea_TT_MOD.Color_Ratio;
% 
% Color_ratio_correct(Color_ratio_correct<0 | Color_ratio_correct>2) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
% Southernmost_BSea_TT_MOD = addvars(Southernmost_BSea_TT_MOD, Color_ratio_correct);
%  

% Remove outlier MAOD values
MAOD_filtered = Southernmost_BSea_TT_MOD.CMOD_Surface;
MAOD_filtered(MAOD_filtered<0 | MAOD_filtered > 0.5 | MAOD_filtered == 0) = NaN; 

MAOD_nocleanair = Southernmost_BSea_TT_MOD.CMOD_Surface;
MAOD_nocleanair(MAOD_nocleanair<0 |  MAOD_nocleanair == 0) = NaN; 
Southernmost_BSea_TT_MOD = addvars(Southernmost_BSea_TT_MOD, MAOD_filtered, MAOD_nocleanair);
  

% % Remove outlier ice values
Ice_correct = Southernmost_BSea_TT_Ice.Total_Surface_532_Integrated_Depolarization_Ratio ;

Ice_correct(Ice_correct<0 | Ice_correct>5) = NaN; % I might need to verify with Jay later to see if these ratio bounds are appropriate. 
Southernmost_BSea_TT_Ice = addvars(Southernmost_BSea_TT_Ice, Ice_correct);
 


% Daily, 
Southernmost_BSea_TT_Ice_daily = retime(Southernmost_BSea_TT_Ice, 'daily', @nanmean);
Southernmost_BSea_TT_MOD_daily = retime(Southernmost_BSea_TT_MOD, 'daily',@nanmean);
% Southernmost_BSea_TT_amsrmf_daily = retime(Southernmost_BSea_TT_Wind, 'daily', @nanmean);


% annually:

% Southernmost_BSea_TT_Ice_yearly = retime(Southernmost_BSea_TT_Ice, 'yearly', @nanmean);
% Southernmost_BSea_TT_MOD_yearly = retime(Southernmost_BSea_TT_MOD, 'yearly',@nanmean);
% Southernmost_BSea_TT_amsrmf_yearly  = retime(Southernmost_BSea_TT_Wind, 'yearly', @nanmean);


% Synchronizing these daily timetables so they span the same days:
SS = timerange('2007-01-01', '2018-12-31', 'closed'); %closed includes startime and endtime values

Southernmost_timetable_all_daily   = synchronize(Southernmost_BSea_TT_MOD_daily, Southernmost_BSea_TT_Ice_daily);
% Southernmost_timetable_all_yearly  = synchronize(Southernmost_BSea_TT_MOD_yearly, Southernmost_BSea_TT_Ice_yearly, Southernmost_BSea_TT_amsrmf_yearly);

% And then create new timetables with all of these variables:
Southernmost_TT_BSea_daily  = Southernmost_timetable_all_daily(SS,:);


% Open ocean: 
openocean = Southernmost_TT_BSea_daily.Ice_correct ;
openocean(openocean > 0.15) = NaN;
Southernmost_TT_BSea_daily = addvars(Southernmost_TT_BSea_daily, openocean);


MAOD_nocleanair_openocean = Southernmost_TT_BSea_daily.MAOD_nocleanair;
MAOD_nocleanair_openocean(isnan(Southernmost_TT_BSea_daily.openocean)) = NaN;
Southernmost_TT_BSea_daily = addvars(Southernmost_TT_BSea_daily, MAOD_nocleanair_openocean);



%%
times_daily =  Corner_TT_BSea_daily.Total_Profile_Time_New_Surface;
% times_8day = Southernmost_TT_BSea_8day.Total_Profile_Time_New_Surface;
% times_monthly = Corner_TT_BSea_monthly.Total_Profile_Time_New_Surface;


%% 
Corner_TT_OD_Wind_daily = synchronize(Corner_AOD_TT, Corner_TT_BSea_daily);
Central_TT_OD_Wind_daily = synchronize(Central_AOD_TT, Central_TT_BSea_daily); 
Southernmost_TT_OD_Wind_daily = synchronize(Southernmost_AOD_TT,Southernmost_TT_BSea_daily);


Corner_AODc_openocean = Corner_TT_OD_Wind_daily.Corner_AOD_coarse_daily;
Corner_AODc_openocean(isnan(Corner_TT_OD_Wind_daily.openocean)) = NaN; 
Corner_TT_OD_Wind_daily = addvars(Corner_TT_OD_Wind_daily, Corner_AODc_openocean);


Central_AODc_openocean = Central_TT_OD_Wind_daily.Central_AOD_coarse_daily;
Central_AODc_openocean(isnan(Central_TT_OD_Wind_daily.openocean)) = NaN; 
Central_TT_OD_Wind_daily = addvars(Central_TT_OD_Wind_daily, Central_AODc_openocean);


Southernmost_AODc_openocean = Southernmost_TT_OD_Wind_daily.Southernmost_AOD_coarse_daily;
Southernmost_AODc_openocean(isnan(Southernmost_TT_OD_Wind_daily.openocean)) = NaN; 
Southernmost_TT_OD_Wind_daily = addvars(Southernmost_TT_OD_Wind_daily, Southernmost_AODc_openocean);











% Scatterplot Code:

%% Also checking winter vs summer as in Dror et al 2018

% Corner spatial region, daily, subsetting
Winter_TR_Corner_daily_index = month(Corner_TT_OD_Wind_daily.times_2007_2018) >=6 &...
    month(Corner_TT_OD_Wind_daily.times_2007_2018) <=8;

Spring_TR_Corner_daily_index = month(Corner_TT_OD_Wind_daily.times_2007_2018) >=9 &...
    month(Corner_TT_OD_Wind_daily.times_2007_2018) <=11;

Summer_TR_Corner_daily_index = month(Corner_TT_OD_Wind_daily.times_2007_2018) >=1 &...
    month(Corner_TT_OD_Wind_daily.times_2007_2018) <=2 | ...
    month(Corner_TT_OD_Wind_daily.times_2007_2018) == 12;

Fall_TR_Corner_daily_index = month(Corner_TT_OD_Wind_daily.times_2007_2018) >=3 &...
    month(Corner_TT_OD_Wind_daily.times_2007_2018) <=5;

Winter_TT_Corner_BSea_daily = Corner_TT_OD_Wind_daily(Winter_TR_Corner_daily_index,:);
Spring_TT_Corner_BSea_daily = Corner_TT_OD_Wind_daily(Spring_TR_Corner_daily_index,:);
Summer_TT_Corner_BSea_daily = Corner_TT_OD_Wind_daily(Summer_TR_Corner_daily_index,:);
Fall_TT_Corner_BSea_daily  = Corner_TT_OD_Wind_daily(Fall_TR_Corner_daily_index,:);


% Central spatial region, daily, subsetting
Winter_TR_Central_daily_index = month(Central_TT_OD_Wind_daily.times_2007_2018) >=6 &...
    month(Central_TT_OD_Wind_daily.times_2007_2018) <=8;

Spring_TR_Central_daily_index = month(Central_TT_OD_Wind_daily.times_2007_2018) >=9 &...
    month(Central_TT_OD_Wind_daily.times_2007_2018) <=11;

Summer_TR_Central_daily_index = month(Central_TT_OD_Wind_daily.times_2007_2018) >=1 &...
    month(Central_TT_OD_Wind_daily.times_2007_2018) <=2 | ...
    month(Central_TT_OD_Wind_daily.times_2007_2018) == 12;

Fall_TR_Central_daily_index = month(Central_TT_OD_Wind_daily.times_2007_2018) >=3 &...
    month(Central_TT_OD_Wind_daily.times_2007_2018) <=5;

Winter_TT_Central_BSea_daily = Central_TT_OD_Wind_daily(Winter_TR_Central_daily_index,:);
Spring_TT_Central_BSea_daily = Central_TT_OD_Wind_daily(Spring_TR_Central_daily_index,:);
Summer_TT_Central_BSea_daily = Central_TT_OD_Wind_daily(Summer_TR_Central_daily_index,:);
Fall_TT_Central_BSea_daily  = Central_TT_OD_Wind_daily(Fall_TR_Central_daily_index,:);


% Southernmost spatial region, daily, subsetting
Winter_TR_Southernmost_daily_index = month(Southernmost_TT_OD_Wind_daily.times_2007_2018) >=6 &...
    month(Southernmost_TT_OD_Wind_daily.times_2007_2018) <=8;

Spring_TR_Southernmost_daily_index = month(Southernmost_TT_OD_Wind_daily.times_2007_2018) >=9 &...
    month(Southernmost_TT_OD_Wind_daily.times_2007_2018) <=11;

Summer_TR_Southernmost_daily_index = month(Southernmost_TT_OD_Wind_daily.times_2007_2018) >=1 &...
    month(Southernmost_TT_OD_Wind_daily.times_2007_2018) <=2 | ...
    month(Southernmost_TT_OD_Wind_daily.times_2007_2018) == 12;

Fall_TR_Southernmost_daily_index = month(Southernmost_TT_OD_Wind_daily.times_2007_2018) >=3 &...
    month(Southernmost_TT_OD_Wind_daily.times_2007_2018) <=5;

Winter_TT_Southernmost_BSea_daily = Southernmost_TT_OD_Wind_daily(Winter_TR_Southernmost_daily_index,:);
Spring_TT_Southernmost_BSea_daily = Southernmost_TT_OD_Wind_daily(Spring_TR_Southernmost_daily_index,:);
Summer_TT_Southernmost_BSea_daily = Southernmost_TT_OD_Wind_daily(Summer_TR_Southernmost_daily_index,:);
Fall_TT_Southernmost_BSea_daily  = Southernmost_TT_OD_Wind_daily(Fall_TR_Southernmost_daily_index,:);


                                    %% Scatter plot, All Regions, daily data:
                                    % SUBPLOTTED

                                    
                                    %% 

% SCATTERPLOTS BETWEEN ICE & COLOR RATIO AND BETWEEN CHL & COLOR
% RATIO 
%  
% I want the true values in these plots
%                        
% CHL & COLOR RATIO Data: 

% x_Corner_CHL = Corner_aqua_daily_correct;
% % x_Corner_ICE = Corner_TT_BSea_daily.Ice_correct;
% y_Corner_CR = Corner_TT_BSea_daily.Color_ratio_correct;
% 
% nanVals_Corner = ismissing(x_Corner_CHL) | ismissing(y_Corner_CR); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Corner_CHL = x_Corner_CHL(~nanVals_Corner);
% y_Corner_CR = y_Corner_CR(~nanVals_Corner);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Corner = polyfit(x_Corner_CHL,y_Corner_CR,1);
% 
% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Corner = polyval(p_Corner,x_Corner_CHL);
% 
% R_Corner = corrcoef(x_Corner_CHL, y_Corner_CR);


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :
% 
% x_Central_CHL = Central_aqua_daily_correct;
% %x_Central_ICE = Central_TT_BSea_daily.Ice_correct;
% y_Central_CR = Central_TT_BSea_daily.Color_ratio_correct;
% 
% nanVals_Central = ismissing(x_Central_CHL) | ismissing(y_Central_CR); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Central_CHL = x_Central_CHL(~nanVals_Central);
% y_Central_CR = y_Central_CR(~nanVals_Central);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% p_Central = polyfit(x_Central_CHL,y_Central_CR,1);

% % Call polyval to use p to predict y, calling the result yfit:
% yfit_Central = polyval(p_Central,x_Central_CHL);
% 
% R_Central = corrcoef(x_Central_CHL, y_Central_CR );
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
%  x_Southernmost_CHL = Southernmost_aqua_daily_correct;
% %x_Southernmost_ICE = Southernmost_TT_BSea_daily.Ice_correct;
% % y_Southernmost_CR = Southernmost_TT_BSea_daily.Color_ratio_correct;
% 
% nanVals_Southernmost = ismissing(x_Southernmost_CHL) | ismissing(y_Southernmost_CR); % indices of values that are NaN in x or y
% % resampling x and y to exclude any entries with NaN Values
% x_Southernmost_CHL = x_Southernmost_CHL(~nanVals_Southernmost);
% % y_Southernmost_CR = y_Southernmost_CR(~nanVals_Southernmost);
% 
% % Use polyfit to compute a linear regression that predicts y from x:
% % p_Southernmost = polyfit(x_Southernmost_CHL,y_Southernmost_CR,1);
% % 
% % % Call polyval to use p to predict y, calling the result yfit:
% % yfit_Southernmost = polyval(p_Southernmost,x_Southernmost_CHL);
% % 
% % R_Southernmost = corrcoef(x_Southernmost_CHL, y_Southernmost_CR );

                                    %%
clear x y nanVals p yfit R
x_Corner_CHL_Daily_forAODc = Corner_aqua_daily_correct;
y_Corner_AODc_Daily = Corner_TT_OD_Wind_daily.Corner_AODc_openocean;

nanVals_Corner_AODc_CHL = ismissing(x_Corner_CHL_Daily_forAODc) | ismissing(y_Corner_AODc_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_CHL_Daily_forAODc = x_Corner_CHL_Daily_forAODc(~nanVals_Corner_AODc_CHL);
y_Corner_AODc_Daily = y_Corner_AODc_Daily(~nanVals_Corner_AODc_CHL);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_Daily_AODc_CHL = polyfit(x_Corner_CHL_Daily_forAODc,y_Corner_AODc_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_Daily_AODc_CHL = polyval(p_Corner_Daily_AODc_CHL,x_Corner_CHL_Daily_forAODc);

[R_Corner_Daily_AODc_CHL, pval_Corner_Daily_AODc_CHL] = corrcoef(x_Corner_CHL_Daily_forAODc, y_Corner_AODc_Daily );


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_CHL_Daily_forAODc = Central_aqua_daily_correct;
y_Central_AODc_Daily = Central_TT_OD_Wind_daily.Central_AODc_openocean;

nanVals_Central_AODc_CHL = ismissing(x_Central_CHL_Daily_forAODc) | ismissing(y_Central_AODc_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_CHL_Daily_forAODc = x_Central_CHL_Daily_forAODc(~nanVals_Central_AODc_CHL);
y_Central_AODc_Daily = y_Central_AODc_Daily(~nanVals_Central_AODc_CHL);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_Daily_AODc_CHL = polyfit(x_Central_CHL_Daily_forAODc,y_Central_AODc_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_Daily_AODc_CHL = polyval(p_Central_Daily_AODc_CHL,x_Central_CHL_Daily_forAODc);

[R_Central_Daily_AODc_CHL, pval_Central_Daily_AODc_CHL] = corrcoef(x_Central_CHL_Daily_forAODc, y_Central_AODc_Daily );


% Southernmost daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_CHL_Daily_forAODc = Southernmost_aqua_daily_correct;
y_Southernmost_AODc_Daily = Southernmost_TT_OD_Wind_daily.Southernmost_AODc_openocean;

nanVals_Southernmost_AODc_CHL = ismissing(x_Southernmost_CHL_Daily_forAODc) | ismissing(y_Southernmost_AODc_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_CHL_Daily_forAODc = x_Southernmost_CHL_Daily_forAODc(~nanVals_Southernmost_AODc_CHL);
y_Southernmost_AODc_Daily = y_Southernmost_AODc_Daily(~nanVals_Southernmost_AODc_CHL);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_Daily_AODc_CHL = polyfit(x_Southernmost_CHL_Daily_forAODc,y_Southernmost_AODc_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_Daily_AODc_CHL = polyval(p_Southernmost_Daily_AODc_CHL,x_Southernmost_CHL_Daily_forAODc);

[R_Southernmost_Daily_AODc_CHL, pval_Southernmost_Daily_AODc_CHL] = corrcoef(x_Southernmost_CHL_Daily_forAODc, y_Southernmost_AODc_Daily );

%%


clear x y nanVals p yfit R
x_Corner_CHL_Daily_forMAOD = Corner_aqua_daily_correct;
y_Corner_MAOD_Daily = Corner_TT_OD_Wind_daily.MAOD_nocleanair_openocean;

nanVals_Corner_MAOD_CHL = ismissing(x_Corner_CHL_Daily_forMAOD) | ismissing(y_Corner_MAOD_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Corner_CHL_Daily_forMAOD = x_Corner_CHL_Daily_forMAOD(~nanVals_Corner_MAOD_CHL);
y_Corner_MAOD_Daily = y_Corner_MAOD_Daily(~nanVals_Corner_MAOD_CHL);

% Use polyfit to compute a linear regression that predicts y from x:
p_Corner_Daily_MAOD_CHL = polyfit(x_Corner_CHL_Daily_forMAOD,y_Corner_MAOD_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Corner_Daily_MAOD_CHL = polyval(p_Corner_Daily_MAOD_CHL,x_Corner_CHL_Daily_forMAOD);

[R_Corner_Daily_MAOD_CHL, pval_Corner_Daily_MAOD_CHL] = corrcoef(x_Corner_CHL_Daily_forMAOD, y_Corner_MAOD_Daily);


% Central daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

x_Central_CHL_Daily_forMAOD = Central_aqua_daily_correct;
y_Central_MAOD_Daily = Central_TT_OD_Wind_daily.MAOD_nocleanair_openocean;

nanVals_Central_MAOD_CHL = ismissing(x_Central_CHL_Daily_forMAOD) | ismissing(y_Central_MAOD_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Central_CHL_Daily_forMAOD = x_Central_CHL_Daily_forMAOD(~nanVals_Central_MAOD_CHL);
y_Central_MAOD_Daily = y_Central_MAOD_Daily(~nanVals_Central_MAOD_CHL);

% Use polyfit to compute a linear regression that predicts y from x:
p_Central_Daily = polyfit(x_Central_CHL_Daily_forMAOD,y_Central_MAOD_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Central_Daily = polyval(p_Central_Daily,x_Central_CHL_Daily_forMAOD);

[R_Central_Daily_MAOD_CHL, pval_Central_Daily_MAOD_CHL] = corrcoef(x_Central_CHL_Daily_forMAOD, y_Central_MAOD_Daily );


% Southernmost daily scatterplot


% Trying to plot MAOD against windspeed, let's start with corner region.
% I've already plotted out the histogram. :

clear x y
% load count.dat
x_Southernmost_CHL_Daily_forMAOD = Southernmost_aqua_daily_correct;
y_Southernmost_MAOD_Daily = Southernmost_TT_OD_Wind_daily.MAOD_nocleanair_openocean;

nanVals_Southernmost_MAOD_CHL = ismissing(x_Southernmost_CHL_Daily_forMAOD) | ismissing(y_Southernmost_MAOD_Daily); % indices of values that are NaN in x or y
% resampling x and y to exclude any entries with NaN Values
x_Southernmost_CHL_Daily_forMAOD = x_Southernmost_CHL_Daily_forMAOD(~nanVals_Southernmost_MAOD_CHL);
y_Southernmost_MAOD_Daily = y_Southernmost_MAOD_Daily(~nanVals_Southernmost_MAOD_CHL);

% Use polyfit to compute a linear regression that predicts y from x:
p_Southernmost_Daily_MAOD_CHL = polyfit(x_Southernmost_CHL_Daily_forMAOD,y_Southernmost_MAOD_Daily,1);

% Call polyval to use p to predict y, calling the result yfit:
yfit_Southernmost_Daily_MAOD_CHL = polyval(p_Southernmost_Daily_MAOD_CHL,x_Southernmost_CHL_Daily_forMAOD);

[R_Southernmost_Daily_MAOD_CHL, pval_Southernmost_Daily_MAOD_CHL] = corrcoef(x_Southernmost_CHL_Daily_forMAOD, y_Southernmost_MAOD_Daily );


%%
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.015], [0.08 0.05], [0.1 0.03]);
if ~make_it_tight,  clear subplot;  end

fig = figure; clf;

%Corner Daily Scatterplot
subplot(2,3,1)
scatter(x_Corner_CHL_Daily_forMAOD,y_Corner_MAOD_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_CHL_Daily_forMAOD,yfit_Corner_Daily_MAOD_CHL,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 1.5])
set(gca,'xtick',[],'xticklabel',[], 'FontSize', 15)

title('Region A') 

txt = ['y = ' num2str(round(p_Corner_Daily_MAOD_CHL(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_MAOD_CHL(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_MAOD_CHL(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_MAOD_CHL(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_CHL_Daily_forMAOD),3,'significant'))];
Region = ('MAOD');

text(0.7, 0.45, txt,'FontSize', 15)
text(0.7,0.42, R_txt,'FontSize', 15)
text(0.7,0.39,N_txt,'FontSize', 15)
text(0.7, 0.36, p_txt,'FontSize', 15)
text(0.05, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region



% Corner Daily Scatterplot:
subplot(2,3,4)
scatter(x_Corner_CHL_Daily_forAODc,y_Corner_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Corner_CHL_Daily_forAODc,yfit_Corner_Daily_AODc_CHL,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 1.5])
set(gca,'FontSize', 15)


txt = ['y = ' num2str(round(p_Corner_Daily_AODc_CHL(1),3,'significant')) 'x + ' num2str(round(p_Corner_Daily_AODc_CHL(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Corner_Daily_AODc_CHL(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Corner_Daily_AODc_CHL(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Corner_CHL_Daily_forAODc),3,'significant'))];
Region = ('AOD_C');

text(0.7, 0.45, txt,'FontSize', 15)
text(0.7,0.42, R_txt,'FontSize', 15)
text(0.7,0.39,N_txt,'FontSize', 15)
text(0.7, 0.36, p_txt,'FontSize', 15)
text(0.05, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% Central Daily Scatterplot
subplot(2,3,2)
scatter(x_Central_CHL_Daily_forMAOD,y_Central_MAOD_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_CHL_Daily_forMAOD,yfit_Central_Daily,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 1.5])
set(gca,'xtick',[],'xticklabel',[], 'ytick',[], 'FontSize', 15)

title('Region B')

txt = ['y = ' num2str(round(p_Central_Daily(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_MAOD_CHL(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_MAOD_CHL(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_CHL_Daily_forMAOD),3,'significant'))];
Region = ('MAOD');

text(0.7, 0.45, txt,'FontSize', 15)
text(0.7,0.42, R_txt,'FontSize', 15)
text(0.7,0.39,N_txt,'FontSize', 15)
text(0.7, 0.36, p_txt,'FontSize', 15)

text(0.05, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region

% Central Daily Scatterplot
subplot(2,3,5)
scatter(x_Central_CHL_Daily_forAODc,y_Central_AODc_Daily,15,'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Central_CHL_Daily_forAODc,yfit_Central_Daily_AODc_CHL,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 1.5])
set(gca, 'ytick',[], 'FontSize', 15)
xlabel('Daily chlorophyll-{\ita} average (mg m^{-3})', 'FontSize', 20);



txt = ['y = ' num2str(round(p_Central_Daily(1),3,'significant')) 'x + ' num2str(round(p_Central_Daily(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Central_Daily_AODc_CHL(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Central_Daily_AODc_CHL(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Central_CHL_Daily_forAODc),3,'significant'))];
Region = ('AOD_C');

text(0.7, 0.45, txt,'FontSize', 15)
text(0.7,0.42, R_txt,'FontSize', 15)
text(0.7,0.39,N_txt,'FontSize', 15)
text(0.7, 0.36, p_txt,'FontSize', 15)

text(0.05, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')

clear txt R_txt N_txt Region


% Southernmost Daily Scatterplot:
subplot(2,3,3)
scatter(x_Southernmost_CHL_Daily_forMAOD,y_Southernmost_MAOD_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;
hold on
plot(x_Southernmost_CHL_Daily_forMAOD,yfit_Southernmost_Daily_MAOD_CHL,'-',...
    'LineWidth', 2)

ylim([0 0.5])
xlim([0 1.5])
set(gca,'xtick',[],'xticklabel',[],'ytick',[], 'FontSize', 15)
 
% xlabel(sprintf('Daily Wind Speed average (m s^{-1})'));
title('Region C')
txt = ['y = ' num2str(round(p_Southernmost_Daily_MAOD_CHL(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_MAOD_CHL(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_MAOD_CHL(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_MAOD_CHL(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_CHL_Daily_forMAOD),3,'significant'))];
Region = ('MAOD');

text(0.7, 0.45, txt,'FontSize', 15)
text(0.7,0.42, R_txt,'FontSize', 15)
text(0.7,0.39,N_txt,'FontSize', 15)
text(0.7,0.36,p_txt, 'FontSize', 15)
text(0.05, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
clear txt R_txt N_txt Region



% Southernmost Daily Scatterplot:
subplot(2,3,6)
scatter(x_Southernmost_CHL_Daily_forAODc,y_Southernmost_AODc_Daily, 15, 'filled','MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 0.5) ;

hold on
plot(x_Southernmost_CHL_Daily_forAODc,yfit_Southernmost_Daily_AODc_CHL,'-',...
    'LineWidth', 2)
ylim([0 0.5])
xlim([0 1.5])
set(gca,'ytick',[],'FontSize', 15)
 

txt = ['y = ' num2str(round(p_Southernmost_Daily_AODc_CHL(1),3,'significant')) 'x + ' num2str(round(p_Southernmost_Daily_AODc_CHL(2),3,'significant'))];
R_txt = ['R = ' num2str(round(R_Southernmost_Daily_AODc_CHL(1,2),3,'significant'))];
p_txt = ['p = ' num2str(round(pval_Southernmost_Daily_AODc_CHL(1,2),3,'significant'))];
N_txt = ['n = ' num2str(round(length(x_Southernmost_CHL_Daily_forAODc),3,'significant'))];
Region = ('AOD_C');

text(0.7, 0.45, txt,'FontSize', 15)
text(0.7,0.42, R_txt,'FontSize', 15)
text(0.7,0.39,N_txt,'FontSize', 15)
text(0.7,0.36,p_txt, 'FontSize', 15)
text(0.05, 0.46, Region, 'FontSize', 15, 'FontWeight', 'bold')
clear txt R_txt N_txt Region



%%
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PaperTwo_PhD_Causal_Analysis/Revisions/Revision_Figures/Diurnal_winds
set(gcf,'PaperPositionMode','auto')
print(gcf,'NEWWINDS_OPENOCEAN_MAOD_CoarseAOD_CHL_allregions_daily_scatterplot.png','-dpng','-r300');       %  *// 300 dpi


