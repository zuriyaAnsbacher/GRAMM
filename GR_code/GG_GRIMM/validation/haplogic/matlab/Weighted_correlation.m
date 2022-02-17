% ##############################################################################
% # SCRIPT NAME:	Weighted_correlation.m
% # DESCRIPTION:	
% # PARAMETERS:	The two vectors to be correlated y1 and y2, a vector of
% weights 
% # OUTPUT:	Weighted Pearson Correlation Coefficient
% # TABLES:       
% #
% # DATE WRITTEN: 2011-08-03
% # WRITTEN BY:   Abeer Madbouly
% #
% # REVISION HISTORY: 
% # REVISION DATE		REVISED BY	DESCRIPTION 
% # ------- ----------	--------------	-------------------------------------
% #
% #       COPYRIGHT (C) 2008 NATIONAL MARROW DONOR PROGRAM.  
% #               ALL RIGHTS RESERVED        
% ##############################################################################

function wt_corr=Weighted_correlation(y1,y2,weights)


% Calculate Weighted Average for both values
wt_sum1=0;
wt_sum2=0;
for i=1:10
    wt_sum1=wt_sum1 + weights(i)*y1(i);
    wt_sum2=wt_sum2 + weights(i)*y2(i);
end
Weighted_avg1=wt_sum1/sum(weights);
Weighted_avg2=wt_sum2/sum(weights);

%Calcualte Weighted Covariance for both
wt_cov=0;
for i=1:10
    wt_cov=wt_cov + (weights(i)*(y1(i)-Weighted_avg1)*(y2(i)-Weighted_avg2));
end
wt_cov=wt_cov/sum(weights); %sum(weights)=1

%Calcualte Weighted Correlation
    %Calcualte Weighted Covariance for first variable
    wt_cov1=0;
    for i=1:10
        wt_cov1=wt_cov1 + (weights(i)*(y1(i)-Weighted_avg1)*(y1(i)-Weighted_avg1));
    end
    %Calcualte Weighted Covariance for second variable
    wt_cov2=0;
    for i=1:10
        wt_cov2=wt_cov2 + (weights(i)*(y2(i)-Weighted_avg2)*(y2(i)-Weighted_avg2));
    end
    
    wt_corr= wt_cov/sqrt(wt_cov1*wt_cov2);
    
    