
%% 

% Here I'm plotting Seasonal Climatologies for all variables, coarse-mode
% AOD, MAOD, Chl-a, sea ice, wind speed, and SST

%%

clear;
cd /Volumes/'G-DRIVE mobile USB'/SST_matfiles/
load('Master_SST_2013_2018.mat');

%% Pulling MAOD & Ice timetables.

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Variables/Chlorophyll_Data_matfiles
load('Master_chl_a_fullres_2007_2018.mat')
load('Longitude_subset.mat')
load('Latitude_subset.mat')
load('times_2007_2018.mat')

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Variables/AOD_matfiles
load('Master_AOD_coarse_2007_2018.mat')
% load('Master_AOD_fine_2007_2018.mat')
% load('Master_AOD_total_2007_2018.mat')
load('Master_AOD_total_pixelcounts_2007_2018.mat')
load('Dimensions_AOD.mat')

%%

t1 = datetime(2013,01,01);
t2 = datetime(2018,12,31);
times_2013_2018 = t1:caldays(1):t2; 
times_2013_2018 = times_2013_2018';
clear t1 t2

[winter_x, ~] = find(times_2013_2018.Month >= 6 & times_2013_2018.Month <= 8);
% Can check if this worked by typing 'times(winter_x)' in command window
times_2013_2018(winter_x)
SST_winter = Master_SST_2013_2018; 
SST_winter = SST_winter(:,:,winter_x);



% For Spring: 
[spring_x, ~] = find(times_2013_2018.Month >= 9 & times_2013_2018.Month <= 11);
times_2013_2018(spring_x)
SST_spring = Master_SST_2013_2018; 
SST_spring = SST_spring(:,:, spring_x); 


% For Summer: 

[summer_x, ~] = find(times_2013_2018.Month >= 1 & times_2013_2018.Month <=2 | times_2013_2018.Month == 12);

times_2013_2018(summer_x)
SST_summer = Master_SST_2013_2018; 
SST_summer = SST_summer(:,:, summer_x); 


% For Fall: 

[fall_x, ~] = find(times_2013_2018.Month >= 3 & times_2013_2018.Month <= 5); 
times_2013_2018(fall_x) 
SST_fall = Master_SST_2013_2018; 
SST_fall = SST_fall(:,:, fall_x); 

SST_winter_mean = mean(SST_winter,3 ,'omitnan');
SST_spring_mean = mean(SST_spring, 3 , 'omitnan'); 
SST_summer_mean = mean(SST_summer, 3, 'omitnan'); 
SST_fall_mean = mean(SST_fall, 3, 'omitnan'); 

% To have more memory otherwise MATLAB crashes
clear Master_SST_2013_2018

%% Plotting of seasonal climatologies: 

Aqua_Lat = (linspace(Latitude_subset(1), Latitude_subset(end), 40)) ; 
Aqua_Lon = (linspace(Longitude_subset(1), Longitude_subset(end), 46)); 


fun = @(block_struct) nanmean(block_struct.data, [ 1 2]); 

SST_winter_mean_one_degree_res = blockproc(SST_winter_mean,...
    [length(SST_winter_mean(:,1)) ./ length(Aqua_Lon) length(SST_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

SST_spring_mean_one_degree_res = blockproc(SST_spring_mean,...
    [length(SST_winter_mean(:,1)) ./ length(Aqua_Lon) length(SST_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

SST_summer_mean_one_degree_res = blockproc(SST_summer_mean,...
    [length(SST_winter_mean(:,1)) ./ length(Aqua_Lon) length(SST_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

SST_fall_mean_one_degree_res = blockproc(SST_fall_mean,...
    [length(SST_winter_mean(:,1)) ./ length(Aqua_Lon) length(SST_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 
% 


SST_summer_mean_rot = flipud(rot90(SST_summer_mean_one_degree_res)); 
SST_spring_mean_rot = flipud(rot90(SST_spring_mean_one_degree_res)); 
SST_fall_mean_rot   = flipud(rot90(SST_fall_mean_one_degree_res)); 
SST_winter_mean_rot = flipud(rot90(SST_winter_mean_one_degree_res)); 


%%

[winter_x, winter_y] = find(times_2007_2018.Month >= 6 & times_2007_2018.Month <= 8);
% Can check if this worked by typing 'times(winter_x)' in command window
times_2007_2018(winter_x)
Chl_a_winter = Master_chl_a_fullres_2007_2018; 
Chl_a_winter = Chl_a_winter(:,:,winter_x);



% For Spring: 
[spring_x, spring_y] = find(times_2007_2018.Month >= 9 & times_2007_2018.Month <= 11);
times_2007_2018(spring_x)
Chl_a_spring = Master_chl_a_fullres_2007_2018; 
Chl_a_spring = Chl_a_spring(:,:, spring_x); 


% For Summer: 

[summer_x, summer_y] = find(times_2007_2018.Month >= 1 & times_2007_2018.Month <=2 | times_2007_2018.Month == 12);

times_2007_2018(summer_x)
Chl_a_summer = Master_chl_a_fullres_2007_2018; 
Chl_a_summer = Chl_a_summer(:,:, summer_x); 


% For Fall: 

[fall_x, fall_y] = find(times_2007_2018.Month >= 3 & times_2007_2018.Month <= 5); 
times_2007_2018(fall_x) 
Chl_a_fall = Master_chl_a_fullres_2007_2018; 
Chl_a_fall = Chl_a_fall(:,:, fall_x); 


Chl_a_winter_mean = mean(Chl_a_winter,3 ,'omitnan');
Chl_a_spring_mean = mean(Chl_a_spring, 3 , 'omitnan'); 
Chl_a_summer_mean = mean(Chl_a_summer, 3, 'omitnan'); 
Chl_a_fall_mean = mean(Chl_a_fall, 3, 'omitnan'); 

%% adjust resolution

fun = @(block_struct) nanmean(block_struct.data, [ 1 2]); 

Chl_a_winter_mean_one_degree_res = blockproc(Chl_a_winter_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

Chl_a_spring_mean_one_degree_res = blockproc(Chl_a_spring_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

Chl_a_summer_mean_one_degree_res = blockproc(Chl_a_summer_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 

Chl_a_fall_mean_one_degree_res = blockproc(Chl_a_fall_mean,...
    [length(Chl_a_winter_mean(:,1)) ./ length(Aqua_Lon) length(Chl_a_winter_mean(1,:)) ./ length(Aqua_Lat)], fun); 
% 


Chl_a_summer_mean_rot = flipud(rot90(Chl_a_summer_mean_one_degree_res)); 
Chl_a_spring_mean_rot = flipud(rot90(Chl_a_spring_mean_one_degree_res)); 
Chl_a_fall_mean_rot   = flipud(rot90(Chl_a_fall_mean_one_degree_res)); 
Chl_a_winter_mean_rot = flipud(rot90(Chl_a_winter_mean_one_degree_res)); 



%%

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Variables/New_WindSpeed_WindDirection_data/Srishti_AOD
load('WINDS_TT_NoNAN_minutely.mat')


cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/PhD_Phase_Two_SouthernOcean/Variables/2017_2020_vars

load('Total_timetable_SO_Depol_Ratio_NEW.mat') % this is with all data from 2020 (earlier it was excluding nov and dec months lol)
load('Total_timetable_SO_MOD_NEW.mat')% this is with all data from 2020 (earlier it was excluding nov and dec months lol)

% cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/New_mat_files_CHL_A_DATA_MONTHLY
% 
% load('Aqua_Longitude_Subset_BellingshausenSea.mat') 
% load('Aqua_Latitude_Subset_BellingshausenSea.mat')

Aqua_Lat = (linspace(Latitude_subset(1), Latitude_subset(end), 40)) ; 
Aqua_Lon = (linspace(Longitude_subset(1), Longitude_subset(end), 46)); 
Aqua_Lat = flip(Aqua_Lat); % this is for plotting purposes because pcolor flips image on y axis. 

%% This is just to check how many boxes to divide Aqua_Lat and Aqua_Lon into
nAqua_Lat = Latitude_subset(1) : -1 : Latitude_subset(end) ; 
nAqua_Lon = Longitude_subset(1) : 1 : Longitude_subset(end); 

nAqua_Lat = flip(nAqua_Lat) ; 
% clear Latitude_subset Longitude_subset
%%
cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis


%%

MAOD_nocleanair = Total_timetable_SO_MOD_NEW.CMOD_Surface;
MAOD_nocleanair(MAOD_nocleanair<0 |  MAOD_nocleanair == 0) = NaN; 
Total_timetable_SO_MOD_NEW = addvars(Total_timetable_SO_MOD_NEW, MAOD_nocleanair);
 
%%
SS = timerange('2007-01-01', '2018-12-31', 'closed'); %closed includes startime and endtime values

% And then create new timetables with all of these variables:
TT_WINDS   = WINDS_TT_NEW_minutely(SS,:);
TT_MAOD_CR = Total_timetable_SO_MOD_NEW(SS,:);
TT_Ice = Total_timetable_SO_Depol_Ratio_NEW(SS,:);

% cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/Updated_Code_Processing_PhdPhaseOne/Srishti/Analysis_and_Vars_For_Publication/Revisions_Figures

% clearvars -except Latitude_subset Longitude_subset Aqua_Lat Aqua_Lon times TT_amsrmf TT_MAOD_CR TT_Ice

%% CMOD, Winds, & Ice Now 

   % CALIPSO didn't start having data until June of 2006. This winter season will be shorter. 
   
%     TR_winter_2006 = timerange('2006-06-01', '2006-09-01'); % Winter: June 1st to September 1st
    TR_winter_2007 = timerange('2007-06-01', '2007-09-01'); 
    TR_winter_2008 = timerange('2008-06-01', '2008-09-01');
    TR_winter_2009 = timerange('2009-06-01', '2009-09-01');
    TR_winter_2010 = timerange('2010-06-01', '2010-09-01'); 
    TR_winter_2011 = timerange('2011-06-01', '2011-09-01');
    TR_winter_2012 = timerange('2012-06-01', '2012-09-01');
    TR_winter_2013 = timerange('2013-06-01', '2013-09-01');
    TR_winter_2014 = timerange('2014-06-01', '2014-09-01');
    TR_winter_2015 = timerange('2015-06-01', '2015-09-01');
    TR_winter_2016 = timerange('2016-06-01', '2016-09-01');
    TR_winter_2017 = timerange('2017-06-01', '2017-09-01');
    TR_winter_2018 = timerange('2018-06-01', '2018-09-01'); 
    
%     TR_spring_2006 = timerange('2006-09-01', '2006-12-01'); % Spring: Sept 1st to Dec 1st
    TR_spring_2007 = timerange('2007-09-01', '2007-12-01');
    TR_spring_2008 = timerange('2008-09-01', '2008-12-01'); 
    TR_spring_2009 = timerange('2009-09-01', '2009-12-01'); 
    TR_spring_2010 = timerange('2010-09-01', '2010-12-01'); 
    TR_spring_2011 = timerange('2011-09-01', '2011-12-01');
    TR_spring_2012 = timerange('2012-09-01', '2012-12-01'); 
    TR_spring_2013 = timerange('2013-09-01', '2013-12-01'); 
    TR_spring_2014 = timerange('2014-09-01', '2014-12-01');
    TR_spring_2015 = timerange('2015-09-01', '2015-12-01'); 
    TR_spring_2016 = timerange('2016-09-01', '2016-12-01'); 
    TR_spring_2017 = timerange('2017-09-01', '2017-12-01');
    TR_spring_2018 = timerange('2018-09-01', '2018-12-01'); 
    
%     TR_summer_2006 = timerange('2006-12-01', '2007-03-01'); % Summer: Dec 1st to March 1st
    TR_summer_2007 = timerange('2007-12-01', '2008-03-01');
    TR_summer_2008 = timerange('2008-12-01', '2009-03-01'); 
    TR_summer_2009 = timerange('2009-12-01', '2010-03-01'); 
    TR_summer_2010 = timerange('2010-12-01', '2011-03-01'); 
    TR_summer_2011 = timerange('2011-12-01', '2012-03-01');
    TR_summer_2012 = timerange('2012-12-01', '2013-03-01'); 
    TR_summer_2013 = timerange('2013-12-01', '2014-03-01'); 
    TR_summer_2014 = timerange('2014-12-01', '2015-03-01');
    TR_summer_2015 = timerange('2015-12-01', '2016-03-01'); 
    TR_summer_2016 = timerange('2016-12-01', '2017-03-01'); 
    TR_summer_2017 = timerange('2017-12-01', '2018-03-01');
    TR_summer_2018 = timerange('2018-12-01', '2019-01-01'); 
    
    %     TR_fall_2006   = timerange('2006-03-01', '2006-06-01'); % there
    %     is no CALIPSO DATA in fall of 2006
    TR_fall_2007   = timerange('2007-03-01', '2007-06-01'); %% Fall: March 1st to June 1st
    TR_fall_2008   = timerange('2008-03-01', '2008-06-01');
    TR_fall_2009   = timerange('2009-03-01', '2009-06-01');
    TR_fall_2010   = timerange('2010-03-01', '2010-06-01');
    TR_fall_2011   = timerange('2011-03-01', '2011-06-01');
    TR_fall_2012   = timerange('2012-03-01', '2012-06-01');
    TR_fall_2013   = timerange('2013-03-01', '2013-06-01');
    TR_fall_2014   = timerange('2014-03-01', '2014-06-01');
    TR_fall_2015   = timerange('2015-03-01', '2015-06-01');
    TR_fall_2016   = timerange('2016-03-01', '2016-06-01');
    TR_fall_2017   = timerange('2017-03-01', '2017-06-01');
    TR_fall_2018   = timerange('2018-03-01', '2018-06-01');

%%
% Initializing of my variables for the loop below. That way on the first
% iteration there is something to fill in (see end of loop to understand).

% MAOD:
Total_MAOD_winter = [];
Total_Lat_MAOD_winter = [];
Total_Lon_MAOD_winter = [];

Total_MAOD_spring = [];
Total_Lat_MAOD_spring = [];
Total_Lon_MAOD_spring = [];

Total_MAOD_summer = [];
Total_Lat_MAOD_summer = [];
Total_Lon_MAOD_summer = [];

Total_MAOD_fall = [];
Total_Lat_MAOD_fall = [];
Total_Lon_MAOD_fall = [];

% Depolarization Ratio:
Total_Ice_winter = [];
Total_Lat_Ice_winter = [];
Total_Lon_Ice_winter = [];

Total_Ice_spring = [];
Total_Lat_Ice_spring = [];
Total_Lon_Ice_spring = [];

Total_Ice_summer = [];
Total_Lat_Ice_summer = [];
Total_Lon_Ice_summer = [];

Total_Ice_fall = [];
Total_Lat_Ice_fall = [];
Total_Lon_Ice_fall = [];

%Wind:
Total_Wind_winter = [];
Total_Lat_Wind_winter = [];
Total_Lon_Wind_winter = [];

Total_Wind_spring = [];
Total_Lat_Wind_spring = [];
Total_Lon_Wind_spring = [];

Total_Wind_summer = [];
Total_Lat_Wind_summer = [];
Total_Lon_Wind_summer = [];

Total_Wind_fall = [];
Total_Lat_Wind_fall = [];
Total_Lon_Wind_fall = [];


for i = 2007:2018
    
    disp(i)
    
    SEASON = {'winter', 'spring', 'summer', 'fall'};
    
    for j = 1:length(SEASON) % 1:4
        
        disp(j)
        
        if exist(sprintf('TR_%s_%d', SEASON{j}, i),'var') % only keep looping if the variable actually exists
            
            eval(sprintf('Lat_MAOD = TT_MAOD_CR(TR_%s_%d,:).Total_Latitude_Surface;', SEASON{j},i))
            eval(sprintf('Lon_MAOD = TT_MAOD_CR(TR_%s_%d, :).Total_Longitude_Surface;',SEASON{j}, i))
            eval(sprintf('MAOD     = TT_MAOD_CR(TR_%s_%d, :).CMOD_Surface;', SEASON{j}, i))
            
            eval(sprintf('Lat_Ice = TT_Ice(TR_%s_%d,:).Total_Latitude_Ice;', SEASON{j},i))
            eval(sprintf('Lon_Ice = TT_Ice(TR_%s_%d,:).Total_Longitude_Ice;', SEASON{j},i))
            eval(sprintf('Ice     = TT_Ice(TR_%s_%d,:).Total_Surface_532_Integrated_Depolarization_Ratio;', SEASON{j},i))
            
            bad_Ice_values = Ice <= -0.2 | Ice > 2;
            Ice(bad_Ice_values) = NaN;
            
            nan_ice        = isnan(Ice(:,1));
            Ice      = Ice(~nan_ice) ;
            Lat_Ice  = Lat_Ice(~nan_ice);
            Lon_Ice  = Lon_Ice(~nan_ice);
            
            
            eval(sprintf('Lat_Wind = TT_WINDS(TR_%s_%d,:).MASTER_Latitude;', SEASON{j}, i))
            eval(sprintf('Lon_Wind = TT_WINDS(TR_%s_%d,:).MASTER_Longitude;', SEASON{j}, i))
            eval(sprintf('Wind     = TT_WINDS(TR_%s_%d,:).MASTER_Winds;', SEASON{j}, i))
            
            
            if j == 1
                
                Total_MAOD_winter     = vertcat(Total_MAOD_winter,MAOD);
                Total_Lat_MAOD_winter = vertcat(Total_Lat_MAOD_winter, Lat_MAOD);
                Total_Lon_MAOD_winter = vertcat(Total_Lon_MAOD_winter, Lon_MAOD);
                                
                
                Total_Ice_winter     = vertcat(Total_Ice_winter, Ice);
                Total_Lat_Ice_winter = vertcat(Total_Lat_Ice_winter, Lat_Ice);
                Total_Lon_Ice_winter = vertcat(Total_Lon_Ice_winter, Lon_Ice);
                
                Total_Wind_winter     = vertcat(Total_Wind_winter, Wind);
                Total_Lat_Wind_winter = vertcat(Total_Lat_Wind_winter, Lat_Wind);
                Total_Lon_Wind_winter = vertcat(Total_Lon_Wind_winter, Lon_Wind);
                
            elseif j == 2
                
                Total_MAOD_spring     = vertcat(Total_MAOD_spring,MAOD);
                Total_Lat_MAOD_spring = vertcat(Total_Lat_MAOD_spring, Lat_MAOD);
                Total_Lon_MAOD_spring = vertcat(Total_Lon_MAOD_spring, Lon_MAOD);
                
                
                Total_Ice_spring     = vertcat(Total_Ice_spring, Ice);
                Total_Lat_Ice_spring = vertcat(Total_Lat_Ice_spring, Lat_Ice);
                Total_Lon_Ice_spring = vertcat(Total_Lon_Ice_spring, Lon_Ice);
                
                Total_Wind_spring = vertcat(Total_Wind_spring, Wind);
                Total_Lat_Wind_spring = vertcat(Total_Lat_Wind_spring, Lat_Wind);
                Total_Lon_Wind_spring = vertcat(Total_Lon_Wind_spring, Lon_Wind);
                
            elseif j == 3
                
                Total_MAOD_summer = vertcat(Total_MAOD_summer,MAOD);
                Total_Lat_MAOD_summer = vertcat(Total_Lat_MAOD_summer, Lat_MAOD);
                Total_Lon_MAOD_summer = vertcat(Total_Lon_MAOD_summer, Lon_MAOD);
                                 
                Total_Ice_summer = vertcat(Total_Ice_summer, Ice);
                Total_Lat_Ice_summer = vertcat(Total_Lat_Ice_summer, Lat_Ice);
                Total_Lon_Ice_summer = vertcat(Total_Lon_Ice_summer, Lon_Ice);
                
                Total_Wind_summer = vertcat(Total_Wind_summer, Wind);
                Total_Lat_Wind_summer = vertcat(Total_Lat_Wind_summer, Lat_Wind);
                Total_Lon_Wind_summer = vertcat(Total_Lon_Wind_summer, Lon_Wind);
                
            elseif j == 4
                
                Total_MAOD_fall = vertcat(Total_MAOD_fall,MAOD);
                Total_Lat_MAOD_fall = vertcat(Total_Lat_MAOD_fall, Lat_MAOD);
                Total_Lon_MAOD_fall = vertcat(Total_Lon_MAOD_fall, Lon_MAOD);
                
                
                Total_Ice_fall = vertcat(Total_Ice_fall, Ice);
                Total_Lat_Ice_fall = vertcat(Total_Lat_Ice_fall, Lat_Ice);
                Total_Lon_Ice_fall = vertcat(Total_Lon_Ice_fall, Lon_Ice);
                
                Total_Wind_fall = vertcat(Total_Wind_fall, Wind);
                Total_Lat_Wind_fall = vertcat(Total_Lat_Wind_fall, Lat_Wind);
                Total_Lon_Wind_fall = vertcat(Total_Lon_Wind_fall, Lon_Wind);
                
            end
            
            clear MAOD CR Wind Ice Lat_MAOD Lon_MAOD Lat_CR Lon_CR Lat_Wind Lon_Wind Lat_Ice Lon_Ice
            
        end
        
        
    end
end

%%

[MAOD_one_degree_winter, MAOD_OCC_winter, MAOD_STD_winter, MAOD_ERR_winter]  = hist_wt_occ_tot(Total_Lat_MAOD_winter, Total_Lon_MAOD_winter, Total_MAOD_winter, Aqua_Lat', Aqua_Lon');
[MAOD_one_degree_spring, MAOD_OCC_spring, MAOD_STD_spring, MAOD_ERR_spring]  = hist_wt_occ_tot(Total_Lat_MAOD_spring, Total_Lon_MAOD_spring, Total_MAOD_spring, Aqua_Lat', Aqua_Lon');
[MAOD_one_degree_summer, MAOD_OCC_summer, MAOD_STD_summer, MAOD_ERR_summer]  = hist_wt_occ_tot(Total_Lat_MAOD_summer, Total_Lon_MAOD_summer, Total_MAOD_summer, Aqua_Lat', Aqua_Lon');
[MAOD_one_degree_fall, MAOD_OCC_fall, MAOD_STD_fall, MAOD_ERR_fall]          = hist_wt_occ_tot(Total_Lat_MAOD_fall, Total_Lon_MAOD_fall, Total_MAOD_fall, Aqua_Lat', Aqua_Lon');

[Ice_one_degree_winter, Ice_OCC_winter, Ice_STD_winter, Ice_ERR_winter]  = hist_wt_occ_tot(Total_Lat_Ice_winter, Total_Lon_Ice_winter, Total_Ice_winter, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_spring, Ice_OCC_spring, Ice_STD_spring, Ice_ERR_spring]  = hist_wt_occ_tot(Total_Lat_Ice_spring, Total_Lon_Ice_spring, Total_Ice_spring, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_summer, Ice_OCC_summer, Ice_STD_summer, Ice_ERR_summer]  = hist_wt_occ_tot(Total_Lat_Ice_summer, Total_Lon_Ice_summer, Total_Ice_summer, Aqua_Lat', Aqua_Lon');
[Ice_one_degree_fall, Ice_OCC_fall, Ice_STD_fall, Ice_ERR_fall]          = hist_wt_occ_tot(Total_Lat_Ice_fall, Total_Lon_Ice_fall, Total_Ice_fall, Aqua_Lat', Aqua_Lon');

[Wind_one_degree_winter, Wind_OCC_winter, Wind_STD_winter, Wind_ERR_winter]  = hist_wt_occ_tot(Total_Lat_Wind_winter, Total_Lon_Wind_winter, Total_Wind_winter, Aqua_Lat', Aqua_Lon');
[Wind_one_degree_spring, Wind_OCC_spring, Wind_STD_spring, Wind_ERR_spring]  = hist_wt_occ_tot(Total_Lat_Wind_spring, Total_Lon_Wind_spring, Total_Wind_spring, Aqua_Lat', Aqua_Lon');
[Wind_one_degree_summer, Wind_OCC_summer, Wind_STD_summer, Wind_ERR_summer]  = hist_wt_occ_tot(Total_Lat_Wind_summer, Total_Lon_Wind_summer, Total_Wind_summer, Aqua_Lat', Aqua_Lon');
[Wind_one_degree_fall, Wind_OCC_fall, Wind_STD_fall, Wind_ERR_fall]          = hist_wt_occ_tot(Total_Lat_Wind_fall, Total_Lon_Wind_fall, Total_Wind_fall, Aqua_Lat', Aqua_Lon');



%%
% Only select pixels with greater than or equal to X number of observations

pixel_threshold_2007_2018 = Master_AOD_total_pixelcounts_2007_2018 >= 2;
% threshold_total_AOD_2007_2018 = Master_AOD_total_2007_2018;
threshold_coarse_AOD_2007_2018 = Master_AOD_coarse_2007_2018;
% threshold_fine_AOD_2007_2018 = Master_AOD_fine_2007_2018; 

% threshold_total_AOD_2007_2018(~pixel_threshold_2007_2018) = NaN;
threshold_coarse_AOD_2007_2018(~pixel_threshold_2007_2018) = NaN;
% threshold_fine_AOD_2007_2018(~pixel_threshold_2007_2018) = NaN;

% Any negative values should also get thrown out (as NAN)

threshold_coarse_AOD_2007_2018(threshold_coarse_AOD_2007_2018<0) = NaN;
% threshold_fine_AOD_2007_2018(threshold_fine_AOD_2007_2018<0) = NaN;
% threshold_total_AOD_2007_2018(threshold_total_AOD_2007_2018<0) = NaN;



t1 = datetime(2007,01,01);
t2 = datetime(2018,12,31);
times_2007_2018 = t1:caldays(1):t2; 
times_2007_2018 = times_2007_2018';
clear t1 t2




[winter_x, winter_y] = find(times_2007_2018.Month >= 6 & times_2007_2018.Month <= 8);
% Can check if this worked by typing 'times(winter_x)' in command window
times_2007_2018(winter_x)
% 
% AOD_total_winter = threshold_total_AOD_2007_2018; 
% AOD_total_winter = AOD_total_winter(:,:,winter_x);
% AOD_fine_winter = threshold_fine_AOD_2007_2018;
% AOD_fine_winter = AOD_fine_winter(:,:,winter_x);
AOD_coarse_winter = threshold_coarse_AOD_2007_2018;
AOD_coarse_winter = AOD_coarse_winter(:,:,winter_x);

Pixel_Counts_winter = Master_AOD_total_pixelcounts_2007_2018(:,:,winter_x);

% For Spring: 
[spring_x, spring_y] = find(times_2007_2018.Month >= 9 & times_2007_2018.Month <= 11);
times_2007_2018(spring_x)
% AOD_total_spring = threshold_total_AOD_2007_2018; 
% AOD_total_spring = AOD_total_spring(:,:, spring_x); 
% AOD_fine_spring = threshold_fine_AOD_2007_2018;
% AOD_fine_spring = AOD_fine_spring(:,:,spring_x);
AOD_coarse_spring = threshold_coarse_AOD_2007_2018;
AOD_coarse_spring = AOD_coarse_spring(:,:,spring_x);

Pixel_Counts_spring = Master_AOD_total_pixelcounts_2007_2018(:,:,spring_x);

% For Summer: 

[summer_x, summer_y] = find(times_2007_2018.Month >= 1 & times_2007_2018.Month <=2 | times_2007_2018.Month == 12);

times_2007_2018(summer_x)
% AOD_total_summer = threshold_total_AOD_2007_2018; 
% AOD_total_summer = AOD_total_summer(:,:, summer_x); 
% AOD_fine_summer = threshold_fine_AOD_2007_2018;
% AOD_fine_summer = AOD_fine_summer(:,:,summer_x);
AOD_coarse_summer = threshold_coarse_AOD_2007_2018;
AOD_coarse_summer = AOD_coarse_summer(:,:,summer_x);


Pixel_Counts_summer = Master_AOD_total_pixelcounts_2007_2018(:,:,summer_x);

% For Fall: 

[fall_x, fall_y] = find(times_2007_2018.Month >= 3 & times_2007_2018.Month <= 5); 
times_2007_2018(fall_x) 
% AOD_total_fall = threshold_total_AOD_2007_2018; 
% AOD_total_fall = AOD_total_fall(:,:, fall_x); 
% AOD_fine_fall = threshold_fine_AOD_2007_2018;
% AOD_fine_fall = AOD_fine_fall(:,:,fall_x);
AOD_coarse_fall = threshold_coarse_AOD_2007_2018;
AOD_coarse_fall = AOD_coarse_fall(:,:,fall_x);

Pixel_Counts_fall = Master_AOD_total_pixelcounts_2007_2018(:,:, fall_x);

% So now for the raw and smoothn climatologies, you can average all of
% these and plot... (for raw) or average all of them and then smoothn on
% the 2D. 
% 
% AOD_total_winter_mean = mean(AOD_total_winter,3 ,'omitnan');
% AOD_fine_winter_mean = mean(AOD_fine_winter, 3, 'omitnan');
AOD_coarse_winter_mean = mean(AOD_coarse_winter, 3, 'omitnan');
Pixel_Counts_winter_mean = mean(Pixel_Counts_winter, 3, 'omitnan');

% AOD_total_spring_mean = mean(AOD_total_spring, 3 , 'omitnan'); 
% AOD_fine_spring_mean = mean(AOD_fine_spring, 3, 'omitnan');
AOD_coarse_spring_mean = mean(AOD_coarse_spring, 3, 'omitnan');
Pixel_Counts_spring_mean = mean(Pixel_Counts_spring, 3, 'omitnan');

% AOD_total_summer_mean = mean(AOD_total_summer, 3, 'omitnan'); 
% AOD_fine_summer_mean = mean(AOD_fine_summer, 3, 'omitnan');
AOD_coarse_summer_mean = mean(AOD_coarse_summer, 3, 'omitnan');
Pixel_Counts_summer_mean = mean(Pixel_Counts_summer, 3, 'omitnan');

% AOD_total_fall_mean = mean(AOD_total_fall, 3, 'omitnan'); 
% AOD_fine_fall_mean = mean(AOD_fine_fall, 3, 'omitnan');
AOD_coarse_fall_mean = mean(AOD_coarse_fall, 3, 'omitnan');
Pixel_Counts_fall_mean = mean(Pixel_Counts_fall, 3, 'omitnan');



%%
make_it_tight = true;

subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], [0.02 0.02], [0.05 0.3]);

if ~make_it_tight,  clear subplot;  end

fig = figure;clf;


ax(1) = subplot(6,4,1);
plot_BSea_figure_subplot_allvars(Chl_a_winter_mean_rot, ...
    Aqua_Lon,...
    flip(Aqua_Lat),...
    [0 1]); 
colormap(ax(1),cmocean('algae'))

title({'Winter'}, 'FontSize',13);


% Regional boxes:
 lat_BSea = [-60*ones(1,21) -63*ones(1,21) -60 -63]; % Outline of a box
 lon_BSea = [-100:-80 -80:-1:-100 -100 -100];

 m_line(lon_BSea,lat_BSea,'linewi',2,'color','k');     % Area outline ...

 
 lat_BSea_chl = [-66*ones(1,21) -69*ones(1,21) -66 -69]; % Outline of a box
 lon_BSea_chl = [-87:-67 -67:-1:-87 -87 -87];

 m_line(lon_BSea_chl,lat_BSea_chl,'linewi',2,'color','k');     % Area outline ...

 lat_BSea_extra = [-73*ones(1,21) -70*ones(1,21) -70 -73]; % Outline of a box
 lon_BSea_extra = [-95:-75 -75:-1:-95 -95 -95];
 m_line(lon_BSea_extra,lat_BSea_extra,'linewi',2,'color','k');     % Area outline ...

 
m_text(-92,-61.5,'A','FontSize',15,'color','k','Rotation', 10,'fontweight','normal');
m_text(-78,-67.5,'B','FontSize',15,'color','k','fontweight','normal');
m_text(-87,-71.6,'C','FontSize',15,'color','k','Rotation', 7,'fontweight','normal');







ax(2) = subplot(6,4,2);
plot_BSea_figure_subplot_allvars(Chl_a_spring_mean_rot, ...
    Aqua_Lon,...
    flip(Aqua_Lat),...
    [0 1]); 
colormap(ax(2),cmocean('algae'))




title({'Spring'}, 'FontSize',13);


ax(3) = subplot(6,4,3);
plot_BSea_figure_subplot_allvars(Chl_a_summer_mean_rot, ...
    Aqua_Lon,...
    flip(Aqua_Lat),...
    [0 1]); 
colormap(ax(3),cmocean('algae'))

title({'Summer'}, 'FontSize',13);


ax(4) = subplot(6,4,4);
plot_BSea_figure_subplot_allvars(Chl_a_fall_mean_rot, ...
    Aqua_Lon,...
    flip(Aqua_Lat),...
    [0 1]); 
colormap(ax(4),cmocean('algae'))
title({'Fall'}, 'FontSize',13);

    
ax(5) = subplot(6,4,5);
plot_BSea_figure_subplot_allvars(Ice_one_degree_winter, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 0.8]); 
colormap(ax(5),cmocean('Ice'))





ax(6) = subplot(6,4,6);
plot_BSea_figure_subplot_allvars(Ice_one_degree_spring, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 0.8]); 
colormap(ax(6),cmocean('Ice'))



ax(7) = subplot(6,4,7);
plot_BSea_figure_subplot_allvars(Ice_one_degree_summer, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 0.8]); 
colormap(ax(7),cmocean('Ice'))




ax(8) = subplot(6,4,8);
plot_BSea_figure_subplot_allvars(Ice_one_degree_fall, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 0.8]); 
colormap(ax(8),cmocean('Ice'))






ax(9) = subplot(6,4,9);
plot_BSea_figure_subplot_allvars(Wind_one_degree_winter, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 10]); 
colormap(ax(9),cmocean('amp'))


ax(10) = subplot(6,4,10);
plot_BSea_figure_subplot_allvars(Wind_one_degree_spring, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 10]); 
colormap(ax(10),cmocean('amp'))

ax(11) = subplot(6,4,11);
plot_BSea_figure_subplot_allvars(Wind_one_degree_summer, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 10]); 
colormap(ax(11),cmocean('amp'))

ax(12) = subplot(6,4,12);
plot_BSea_figure_subplot_allvars(Wind_one_degree_fall, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 10]); 
colormap(ax(12),cmocean('amp'))





ax(13) = subplot(6,4,13);
plot_BSea_figure_subplot_allvars(SST_winter_mean_rot, ...
    Aqua_Lon,...
    flip(Aqua_Lat),...
    [-5 5]); 
colormap(ax(13),cmocean('balance'))


ax(14) = subplot(6,4,14);
plot_BSea_figure_subplot_allvars(SST_spring_mean_rot, ...
    Aqua_Lon,...
    flip(Aqua_Lat),...
    [-5 5]); 
colormap(ax(14),cmocean('balance'))

ax(15) = subplot(6,4,15);
plot_BSea_figure_subplot_allvars(SST_summer_mean_rot, ...
    Aqua_Lon,...
    flip(Aqua_Lat),...
    [-5 5]); 
colormap(ax(15),cmocean('balance'))

ax(16) = subplot(6,4,16);
plot_BSea_figure_subplot_allvars(SST_fall_mean_rot, ...
    Aqua_Lon,...
    flip(Aqua_Lat),...
    [-5 5]); 
colormap(ax(16),cmocean('balance'))




ax(17) = subplot(6,4,17);
plot_BSea_figure_subplot_allvars(MAOD_one_degree_winter, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 0.08]); 
colormap(ax(17),cmocean('tempo'))


ax(18) = subplot(6,4,18);
plot_BSea_figure_subplot_allvars(MAOD_one_degree_spring, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 0.08]); 
colormap(ax(18),cmocean('tempo'))

ax(19) = subplot(6,4,19);
plot_BSea_figure_subplot_allvars(MAOD_one_degree_summer, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 0.08]); 
colormap(ax(19),cmocean('tempo'))

    
ax(20) = subplot(6,4,20);
plot_BSea_figure_subplot_allvars(MAOD_one_degree_fall, ...
    Aqua_Lon,...
    Aqua_Lat,...
    [0 0.08]); 
colormap(ax(20),cmocean('tempo'))




ax(21) = subplot(6,4,21);
plot_BSea_figure_subplot_allvars(AOD_coarse_winter_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.12]); 
colormap(ax(21),cmocean('tempo'))


ax(22) = subplot(6,4,22);
plot_BSea_figure_subplot_allvars(AOD_coarse_spring_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.12]); 
colormap(ax(22),cmocean('tempo'))


ax(23) = subplot(6,4,23);
plot_BSea_figure_subplot_allvars(AOD_coarse_summer_mean, ...
   AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.12]); 
colormap(ax(23),cmocean('tempo'))

    
ax(24) = subplot(6,4,24);
plot_BSea_figure_subplot_allvars(AOD_coarse_fall_mean, ...
    AOD_Longitude_subset',...
    AOD_Latitude_subset',...
    [0 0.12]); 
colormap(ax(24),cmocean('tempo'))



hp4 = get(subplot(6,4,4),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)  hp4(2) 0.015  0.13]);
h.FontWeight = 'bold';
h.FontSize = 12;
hy = ylabel(h, 'Chl-{\it a} (mg m^{-3})', 'FontSize', 12);
% hy.FontSize = 18;

hp4 = get(subplot(6,4,8),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)  hp4(2) 0.015  0.13]);
h.FontWeight = 'bold';
h.FontSize = 12;
hy = ylabel(h, 'Sea Ice (\delta)', 'FontSize', 12);
% hy.FontSize = 18;

hp4 = get(subplot(6,4,12),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)  hp4(2) 0.015  0.13]);
h.FontWeight = 'bold';
h.FontSize = 12;
hy = ylabel(h, sprintf('Wind Speed (m s^{-1})'), 'FontSize', 12);
% hy.FontSize = 18;


hp4 = get(subplot(6,4,16),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)  hp4(2) 0.015  0.13]);
h.FontWeight = 'bold';
h.FontSize = 12;
hy = ylabel(h, 'SST (\circC)', 'FontSize', 12);
% hy.FontSize = 18;

hp4 = get(subplot(6,4,20),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)  hp4(2) 0.015  0.13]);
h.FontWeight = 'bold';
h.FontSize = 12;
hy = ylabel(h, sprintf('MAOD'), 'FontSize', 12);
% hy.FontSize = 18;


hp4 = get(subplot(6,4,24),'Position');
h = colorbar('Position', [hp4(1)+hp4(3)  hp4(2) 0.015  0.13]);
h.FontWeight = 'bold';
h.FontSize = 12;
hy = ylabel(h, sprintf('AOD_C'), 'FontSize', 12);
% hy.FontSize = 18;
%%

cd /Users/srishtidasarathy/Documents/'Documents - Srishti’s MacBook Pro '/Bowman/SouthernOcean_Causal_Analysis/Figures/NEW_WINDS_Figs

set(gcf,'PaperPositionMode','auto')
% set(gcf,'PaperPosition','fillpage') 

print(gcf,'V3_SST_included_Reorganized_withboundaries_Seasonal_Climatology_all_variables','-dpng','-r300');       %  *// 300 dpi



