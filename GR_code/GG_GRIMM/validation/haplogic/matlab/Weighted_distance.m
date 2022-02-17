% ##############################################################################
% # SCRIPT NAME:	Weighted_distance.m
% # DESCRIPTION:	
% # PARAMETERS:	The two vectors to be correlated y1 and y2, a vector of weights
% # OUTPUT:	Weighted city block distance
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

function wt_dist=Weighted_distance(y1,y2,weights)

dist_vector=abs(y1-y2);
wt_dist_vector=[];
for i=1:10
    wt_dist_vector=[wt_dist_vector;weights(i)*dist_vector(i)];
end
wt_dist=sum(wt_dist_vector);
