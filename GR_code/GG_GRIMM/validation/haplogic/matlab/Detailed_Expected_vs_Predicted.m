% ##############################################################################
% # SCRIPT NAME:	Detailed_Expected_vs_Predicted.m
% # DESCRIPTION:	
% # PARAMETERS:	none
% # OUTPUT:	Generate expected vs. predicted plot paramters for haplogic
% # output has more percentile bins and percentage counts
% # TABLES:       
% #
% # DATE WRITTEN: 2011-10-26
% # DATE MODIFIED:
% # WRITTEN BY:   Abeer Madbouly
% #
% # REVISION HISTORY: 
% # REVISION DATE		REVISED BY	DESCRIPTION 
% # ------- ----------	--------------	-------------------------------------
% #
% #       COPYRIGHT (C) 2008 NATIONAL MARROW DONOR PROGRAM.  
% #               ALL RIGHTS RESERVED        
% ##############################################################################
function results_struct = Detailed_Expected_vs_Predicted(likelihoods, matches, match_value)

pop_data=[likelihoods, matches];
results_struct=struct('range_median',{},'Total_in_range',{},'Total_matches', {}, 'Predicted_matches',{} );

%create percentil bucket data
sorted_pop_data=sortrows(pop_data);
total_pop=size(sorted_pop_data,1);
%0%-10%
temp=sorted_pop_data(sorted_pop_data(:,1)>=0 & sorted_pop_data(:,1)<=0.05);
index_0_5 = size(temp,1);
sorted_pop_data_0_5=sorted_pop_data(1:index_0_5, :);
results_struct(1).range_median = 0.025;
results_struct(1).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_0_5(sorted_pop_data_0_5(:,2)== match_value);
results_struct(1).Total_matches = size(temp,1)/total_pop;
results_struct(1).Predicted_matches = results_struct(1).Total_matches/results_struct(1).Total_in_range;

temp=sorted_pop_data(sorted_pop_data(:,1)>0.05 & sorted_pop_data(:,1)<=0.1);
index_0_10 = index_0_5 + size(temp,1);
sorted_pop_data_5_10=sorted_pop_data(index_0_5+1:index_0_10, :);
results_struct(2).range_median = 0.075;
results_struct(2).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_5_10(sorted_pop_data_5_10(:,2)== match_value);
results_struct(2).Total_matches = size(temp,1)/total_pop;
results_struct(2).Predicted_matches = results_struct(2).Total_matches/results_struct(2).Total_in_range;

%10%-20%
temp=sorted_pop_data(sorted_pop_data(:,1)>0.1 & sorted_pop_data(:,1)<=0.15);
index_0_15 = index_0_10 + size(temp,1);
sorted_pop_data_10_15=sorted_pop_data((index_0_10+1):index_0_15, :);
results_struct(3).range_median = 0.125;
results_struct(3).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_10_15(sorted_pop_data_10_15(:,2)== match_value);
results_struct(3).Total_matches = size(temp,1)/total_pop;
results_struct(3).Predicted_matches = results_struct(3).Total_matches/results_struct(3).Total_in_range;

temp=sorted_pop_data(sorted_pop_data(:,1)>0.15 & sorted_pop_data(:,1)<=0.2);
index_0_20 = index_0_15 + size(temp,1);
sorted_pop_data_15_20=sorted_pop_data((index_0_15+1):index_0_20, :);
results_struct(4).range_median = 0.175;
results_struct(4).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_15_20(sorted_pop_data_15_20(:,2)== match_value);
results_struct(4).Total_matches = size(temp,1)/total_pop;
results_struct(4).Predicted_matches = results_struct(4).Total_matches/results_struct(4).Total_in_range;

%20%-30%
temp=sorted_pop_data(sorted_pop_data(:,1)>0.2 & sorted_pop_data(:,1)<=0.25);
index_0_25 = index_0_20 + size(temp,1);
sorted_pop_data_20_25=sorted_pop_data((index_0_20+1):index_0_25, :);
results_struct(5).range_median = 0.225;
results_struct(5).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_20_25(sorted_pop_data_20_25(:,2)== match_value);
results_struct(5).Total_matches = size(temp,1)/total_pop;
results_struct(5).Predicted_matches = results_struct(5).Total_matches/results_struct(5).Total_in_range;

temp=sorted_pop_data(sorted_pop_data(:,1)>0.25 & sorted_pop_data(:,1)<=0.3);
index_0_30 = index_0_25 + size(temp,1);
sorted_pop_data_25_30=sorted_pop_data((index_0_25+1):index_0_30, :);
results_struct(6).range_median = 0.275;
results_struct(6).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_25_30(sorted_pop_data_25_30(:,2)== match_value);
results_struct(6).Total_matches = size(temp,1)/total_pop;
results_struct(6).Predicted_matches = results_struct(6).Total_matches/results_struct(6).Total_in_range;

%30%-40%
temp=sorted_pop_data(sorted_pop_data(:,1)>0.3 & sorted_pop_data(:,1)<=0.35);
index_0_35 = index_0_30 + size(temp,1);
sorted_pop_data_30_35=sorted_pop_data((index_0_30+1):index_0_35, :);
results_struct(7).range_median = 0.325;
results_struct(7).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_30_35(sorted_pop_data_30_35(:,2)== match_value);
results_struct(7).Total_matches = size(temp,1)/total_pop;
results_struct(7).Predicted_matches = results_struct(7).Total_matches/results_struct(7).Total_in_range;

temp=sorted_pop_data(sorted_pop_data(:,1)>0.35 & sorted_pop_data(:,1)<=0.4);
index_0_40 = index_0_35 + size(temp,1);
sorted_pop_data_35_40=sorted_pop_data((index_0_35+1):index_0_40, :);
results_struct(8).range_median = 0.375;
results_struct(8).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_35_40(sorted_pop_data_35_40(:,2)== match_value);
results_struct(8).Total_matches = size(temp,1)/total_pop;
results_struct(8).Predicted_matches = results_struct(8).Total_matches/results_struct(8).Total_in_range;

%40%-50%
temp=sorted_pop_data(sorted_pop_data(:,1)>0.4 & sorted_pop_data(:,1)<=0.45);
index_0_45 = index_0_40 + size(temp,1);
sorted_pop_data_40_45=sorted_pop_data((index_0_40+1):index_0_45, :);
results_struct(9).range_median = 0.425;
results_struct(9).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_40_45(sorted_pop_data_40_45(:,2)== match_value);
results_struct(9).Total_matches = size(temp,1)/total_pop;
results_struct(9).Predicted_matches = results_struct(9).Total_matches/results_struct(9).Total_in_range;

temp=sorted_pop_data(sorted_pop_data(:,1)>0.45 & sorted_pop_data(:,1)<=0.5);
index_0_50 = index_0_45 + size(temp,1);
sorted_pop_data_45_50=sorted_pop_data((index_0_45+1):index_0_50, :);
results_struct(10).range_median = 0.475;
results_struct(10).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_45_50(sorted_pop_data_45_50(:,2)== match_value);
results_struct(10).Total_matches = size(temp,1)/total_pop;
results_struct(10).Predicted_matches = results_struct(10).Total_matches/results_struct(10).Total_in_range;

%50%-60%
temp=sorted_pop_data(sorted_pop_data(:,1)>0.5 & sorted_pop_data(:,1)<=0.55);
index_0_55 = index_0_50 + size(temp,1);
sorted_pop_data_50_55=sorted_pop_data((index_0_50+1):index_0_55, :);
results_struct(11).range_median = 0.525;
results_struct(11).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_50_55(sorted_pop_data_50_55(:,2)== match_value);
results_struct(11).Total_matches = size(temp,1)/total_pop;
results_struct(11).Predicted_matches = results_struct(11).Total_matches/results_struct(11).Total_in_range;

temp=sorted_pop_data(sorted_pop_data(:,1)>0.55 & sorted_pop_data(:,1)<=0.6);
index_0_60 = index_0_55 + size(temp,1);
sorted_pop_data_55_60=sorted_pop_data((index_0_55+1):index_0_60, :);
results_struct(12).range_median = 0.575;
results_struct(12).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_55_60(sorted_pop_data_55_60(:,2)== match_value);
results_struct(12).Total_matches = size(temp,1)/total_pop;
results_struct(12).Predicted_matches = results_struct(12).Total_matches/results_struct(12).Total_in_range;


%60%-70%
temp=sorted_pop_data(sorted_pop_data(:,1)>0.6 & sorted_pop_data(:,1)<=0.65);
index_0_65 = index_0_60 + size(temp,1);
sorted_pop_data_60_65=sorted_pop_data((index_0_60+1):index_0_65, :);
results_struct(13).range_median = 0.625;
results_struct(13).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_60_65(sorted_pop_data_60_65(:,2)== match_value);
results_struct(13).Total_matches = size(temp,1)/total_pop;
results_struct(13).Predicted_matches = results_struct(13).Total_matches/results_struct(13).Total_in_range;

temp=sorted_pop_data(sorted_pop_data(:,1)>0.65 & sorted_pop_data(:,1)<=0.7);
index_0_70 = index_0_65 + size(temp,1);
sorted_pop_data_65_70=sorted_pop_data((index_0_65+1):index_0_70, :);
results_struct(14).range_median = 0.675;
results_struct(14).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_65_70(sorted_pop_data_65_70(:,2)== match_value);
results_struct(14).Total_matches = size(temp,1)/total_pop;
results_struct(14).Predicted_matches = results_struct(14).Total_matches/results_struct(14).Total_in_range;

%70%-80%
temp=sorted_pop_data(sorted_pop_data(:,1)>0.7 & sorted_pop_data(:,1)<=0.75);
index_0_75 = index_0_70 + size(temp,1);
sorted_pop_data_70_75=sorted_pop_data((index_0_70+1):index_0_75, :);
results_struct(15).range_median = 0.725;
results_struct(15).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_70_75(sorted_pop_data_70_75(:,2)== match_value);
results_struct(15).Total_matches = size(temp,1)/total_pop;
results_struct(15).Predicted_matches = results_struct(15).Total_matches/results_struct(15).Total_in_range;

temp=sorted_pop_data(sorted_pop_data(:,1)>0.75 & sorted_pop_data(:,1)<=0.8);
index_0_80 = index_0_75 + size(temp,1);
sorted_pop_data_75_80=sorted_pop_data((index_0_75+1):index_0_80, :);
results_struct(16).range_median = 0.775;
results_struct(16).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_75_80(sorted_pop_data_75_80(:,2)== match_value);
results_struct(16).Total_matches = size(temp,1)/total_pop;
results_struct(16).Predicted_matches = results_struct(16).Total_matches/results_struct(16).Total_in_range;

%80-90%
temp=sorted_pop_data(sorted_pop_data(:,1)>0.8 & sorted_pop_data(:,1)<=0.85);
index_0_85 = index_0_80 + size(temp,1);
sorted_pop_data_80_85=sorted_pop_data((index_0_80+1):index_0_85, :);
results_struct(17).range_median = 0.825;
results_struct(17).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_80_85(sorted_pop_data_80_85(:,2)== match_value);
results_struct(17).Total_matches = size(temp,1)/total_pop;
results_struct(17).Predicted_matches = results_struct(17).Total_matches/results_struct(17).Total_in_range;

temp=sorted_pop_data(sorted_pop_data(:,1)>0.85 & sorted_pop_data(:,1)<=0.9);
index_0_90 = index_0_85 + size(temp,1);
sorted_pop_data_85_90=sorted_pop_data((index_0_85+1):index_0_90, :);
results_struct(18).range_median = 0.875;
results_struct(18).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_85_90(sorted_pop_data_85_90(:,2)== match_value);
results_struct(18).Total_matches = size(temp,1)/total_pop;
results_struct(18).Predicted_matches = results_struct(18).Total_matches/results_struct(18).Total_in_range;

%90-100%
temp=sorted_pop_data(sorted_pop_data(:,1)>0.9 & sorted_pop_data(:,1)<=0.95);
index_0_95 = index_0_90 + size(temp,1);
sorted_pop_data_90_95=sorted_pop_data((index_0_90+1):index_0_95, :);
results_struct(19).range_median = 0.925;
results_struct(19).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_90_95(sorted_pop_data_90_95(:,2)== match_value);
results_struct(19).Total_matches = size(temp,1)/total_pop;
results_struct(19).Predicted_matches = results_struct(19).Total_matches/results_struct(19).Total_in_range;

temp=sorted_pop_data(sorted_pop_data(:,1)>0.95 & sorted_pop_data(:,1)<=1);
index_0_100 = index_0_95 + size(temp,1);
sorted_pop_data_95_100=sorted_pop_data((index_0_95+1):index_0_100, :);
results_struct(20).range_median = 0.975;
results_struct(20).Total_in_range = size(temp,1)/total_pop;
temp=sorted_pop_data_95_100(sorted_pop_data_95_100(:,2)== match_value);
results_struct(20).Total_matches = size(temp,1)/total_pop;
results_struct(20).Predicted_matches = results_struct(20).Total_matches/results_struct(20).Total_in_range;
