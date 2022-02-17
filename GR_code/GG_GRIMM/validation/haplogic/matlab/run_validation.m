% ##############################################################################
% # SCRIPT NAME:	run_validation.m
% # DESCRIPTION:	
% # PARAMETERS:	none
% # OUTPUT:	Automate validation experiments for haplogic
% validation
% # TABLES:       
% #
% # DATE WRITTEN: 2011-5-18
% # WRITTEN BY:   Abeer Madbouly
% #
% # REVISION HISTORY: 
% # REVISION DATE		REVISED BY      DESCRIPTION 
% # 1       2011-8-4    Abeer Madbouly  Added Weighted regressions calcualtions
% # 2       2011-8-22   Abeer Madbouly  Modified output plot and files to be in one .csv file and one figure window
% # 3       2011-8-31   Abeer Madbouly  Added Weighted distance measure 
% # 4       2019-10-10   Abeer Madbouly Added changes to accomodate psudo dids and rids and output results for ONLY 9of10 and 10of10 
% # ------- ----------	--------------	-------------------------------------
% #
% #       COPYRIGHT (C) 2008 NATIONAL MARROW DONOR PROGRAM.  
% #               ALL RIGHTS RESERVED        
% ##############################################################################

function run_validation(filename, filter_drace, roc_thresh, roc_alpha)
    
    [rid,phen,rrace,rethn,hap3_rrace, hap2_rrace,did,drace,dethn,hap3_drace,hap2_drace,ss,rqst_dte, ...
     r_a,dpre_a,dpost_a,a_status,a_p0,a_p1,a_p2,r_b,dpre_b,dpost_b, b_status, ...
     b_p0,b_p1,b_p2,r_c,dpre_c,dpost_c,c_status,c_p0,c_p1,c_p2,r_drb1,dpre_drb1, ...
     dpost_drb1,drb1_status,drb1_p0,drb1_p1,drb1_p2,r_dqb1,dpre_dqb1,dpost_dqb1, ...
     dqb1_status,dqb1_p0,dqb1_p1,dqb1_p2,xof6,p_0of6,p_1of6,p_2of6,p_3of6,p_4of6, ...
     p_5of6,p_6of6,xof8,p_0of8,p_1of8,p_2of8,p_3of8,p_4of8,p_5of8,p_6of8,p_7of8,p_8of8, ...
     xof10,p_0of10,p_1of10,p_2of10,p_3of10,p_4of10,p_5of10,p_6of10,p_7of10,p_8of10,p_9of10, ...
     p_10of10]=textread(filename,'%s %d %s %s %s %s %s %s %s %s %s %s %s %s %s %s %c %f %f %f %s %s %s %c %f %f %f %s %s %s %c %f %f %f %s %s %s %c %f %f %f %s %s %s %c %f %f %f %c %f %f %f %f %f %f %f %c %f %f %f %f %f %f %f %f %f %s %f %f %f %f %f %f %f %f %f %f %f','delimiter','\t');
 
 %Create maximized figure for exporting plots
  h=figure;
  set(gcf, 'Position', get(0, 'Screensize'));
%   maximize(h);
 
 %Create output strcuture
 Validation_struct = struct('locus',{},'AUC',{},'SE',{},'Sensitivity',{},'Specificity',{},'PPV',{},'NPV',{},'wt_corr',{}, 'wt_dist',{});
 Validation_index=1;
 %%%%%%%%%%%%%%%%%%%%%%%
 %Race filter: Possible Values: AFA, API, CAU, DEC, HIS, MLT, NAM, OTH,
 %UNK, ALL
 %%%%%%%%%%%%%%%%%%%%%%%
 if strcmp('ALL', filter_drace) ~= 1
     race_index=[];  
     for i=1:size(hap2_drace,1)
          if strcmp(hap2_drace(i,:), filter_drace) == 1  
              race_index=[race_index;i];
          end       
      end
 end
 
 %p9of10
%  elseif strcmp(locus,'p9of10') == 1
    Validation_struct(Validation_index).locus='9of10'; 
     %Stratify by race
    if strcmp('ALL', filter_drace) ~= 1
        x9of10=xof10(race_index);
        p_9of10=p_9of10(race_index);
    else
        x9of10=xof10;
    end
    %Eliminate Indefinite and discrepant values
    char_xof10=char(x9of10);
    I_D_indeces=[];
    for i=1:size(char_xof10,1)
       if strcmp(char_xof10(i,1), 'I') == 0 & strcmp(char_xof10(i,1), 'D') == 0 
          I_D_indeces=[I_D_indeces;i]; 
       end
    end
    p_9of10_match=p_9of10(I_D_indeces);
    p_9of10_status=char_xof10(I_D_indeces, :);
    
    %Generate ROC matrix
    match_indeces=[];
    mismatch_indeces=[];
    for i=1:size(p_9of10_status,1)
       if strcmp(p_9of10_status(i,1), '9') == 1  
          match_indeces=[match_indeces;i];
       else
          mismatch_indeces=[mismatch_indeces;i]; 
       end
    end       
        
%     roc_input=ones(size(p_9of10_status),1);
    roc_input=ones(size(p_9of10_status));
    roc_input(mismatch_indeces)=0; 
    roc_input=[p_9of10_match roc_input];
    ROCdata=roc_AM(roc_input,roc_thresh, roc_alpha)
   
    %Expected vs. Actual
    Match_values=roc_input(:,2);
    Match_values(match_indeces)=2;
    results_struct = Detailed_Expected_vs_Predicted(p_9of10_match, Match_values, 2);
    x1=[results_struct.range_median]';
    y1=[results_struct.Predicted_matches]';
    y2=[0.025 0.075 0.125 0.175 0.225 0.275 0.325 0.375 0.425 0.475 0.525 0.575 0.625 0.675 0.725 0.775 0.825 0.875 0.925 0.975]';
    y3=[[results_struct.Total_matches]' , [results_struct.Total_in_range]'];
    %corr coeff
%     R = corr(y1,y2);
    
    %weighted corr coeff
    overall_in_range=sum([results_struct.Total_in_range]);
    weights=[results_struct.Total_in_range]'/overall_in_range;
    wt_R=Weighted_correlation(y1,y2,weights);
    %weighted city block distance
    wt_dist=Weighted_distance(y1,y2,weights);
    
%     figure
    [AX,H1,H2] = plotyy(x1,y1,x1,y3,'plot', 'bar');
    set(H1,'LineStyle','-.');
    set(AX(1),'YTick',[0. 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
    set(AX(2),'YTick',[0. 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
    set(AX(2),'XTick',[0. 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1])
    hold on
    axis([0 1 0 1]);
    plot([results_struct.range_median]',y2, 'r');
    title(strcat('Prob 9 of 10 Match Validation - Race = ', filter_drace, ', Corr Coeff = ',  num2str(wt_R), ', dist = ',  num2str(wt_dist)));
    hleg1 = legend('Actual Matches','Expected Matches', 'Total Matches', 'Total in Range');
    set(hleg1,'Location','SouthEastOutside');
    hold off   

    %plot file
    temp=strcat(filename,filter_drace);
    temp=strcat(temp,'_det_wt_dist_P9of10-');
    ImageFileName=strcat(temp, '.png');
    X=getframe(gcf);
    imwrite(X.cdata,ImageFileName,'png');
    %Results

    Validation_struct(Validation_index).AUC=ROCdata.AUC;
    Validation_struct(Validation_index).SE=ROCdata.SE;
    Validation_struct(Validation_index).Sensitivity=ROCdata.Sensitivity;
    Validation_struct(Validation_index).Specificity=ROCdata.Specificity;
    Validation_struct(Validation_index).PPV=ROCdata.PPV;
    Validation_struct(Validation_index).NPV=ROCdata.NPV;
    Validation_struct(Validation_index).wt_corr=wt_R;
    Validation_struct(Validation_index).wt_dist=wt_dist;
    Validation_index=Validation_index+1;
    
 %p10of10
%  elseif strcmp(locus,'p10of10') == 1
    Validation_struct(Validation_index).locus='10of10';
      %Stratify by race
    if strcmp('ALL', filter_drace) ~= 1
        xof10=xof10(race_index);
        p_10of10=p_10of10(race_index);
    end
    %Eliminate Indefinite and discrepant values
    char_xof10=char(xof10);
    I_D_indeces=[];
    for i=1:size(char_xof10,1)
       if strcmp(char_xof10(i,1), 'I') == 0 & strcmp(char_xof10(i,1), 'D') == 0 
          I_D_indeces=[I_D_indeces;i]; 
       end
    end
    p_10of10_match=p_10of10(I_D_indeces);
    p_10of10_status=char_xof10(I_D_indeces, :);
    
    %Generate ROC matrix
    match_indeces=[];
    mismatch_indeces=[];
    for i=1:size(p_10of10_status,1)
       if strcmp(p_10of10_status(i,:), '10') == 1  
          match_indeces=[match_indeces;i];
       else
          mismatch_indeces=[mismatch_indeces;i]; 
       end
    end       
        
    roc_input=ones(size(p_10of10_status));
    roc_input(mismatch_indeces)=0; 
    roc_input=[p_10of10_match roc_input];
    ROCdata=roc_AM(roc_input,roc_thresh, roc_alpha)
    
    %Expected vs. Actual
    Match_values=roc_input(:,2);
    Match_values(match_indeces)=2;
    results_struct = Detailed_Expected_vs_Predicted(p_10of10_match, Match_values, 2);
    x1=[results_struct.range_median]';
    y1=[results_struct.Predicted_matches]';
    y2=[0.025 0.075 0.125 0.175 0.225 0.275 0.325 0.375 0.425 0.475 0.525 0.575 0.625 0.675 0.725 0.775 0.825 0.875 0.925 0.975]';
    y3=[[results_struct.Total_matches]' , [results_struct.Total_in_range]'];
    %corr coeff
%     R = corr(y1,y2);
    
    %weighted corr coeff
    overall_in_range=sum([results_struct.Total_in_range]);
    weights=[results_struct.Total_in_range]'/overall_in_range;
    wt_R=Weighted_correlation(y1,y2,weights);
     %weighted city block distance
    wt_dist=Weighted_distance(y1,y2,weights);
    
%     figure
 set(gcf, 'Position', get(0, 'Screensize'));
    [AX,H1,H2] = plotyy(x1,y1,x1,y3,'plot', 'bar');
    set(H1,'LineStyle','-.');
    set(AX(1),'YTick',[0. 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
    set(AX(2),'YTick',[0. 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
    set(AX(2),'XTick',[0. 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1])
    hold on
    axis([0 1 0 1]);
    plot([results_struct.range_median]',y2, 'r');
    title(strcat('Prob 10 of 10 Match Validation - Race = ', filter_drace, ', Corr Coeff = ',  num2str(wt_R), ', dist = ',  num2str(wt_dist)));
    hleg1 = legend('Actual Matches','Expected Matches', 'Total Matches', 'Total in Range');
    set(hleg1,'Location','SouthEastOutside');
    hold off
    %plot file
    temp=strcat(filename,filter_drace);
    temp=strcat(temp,'_det_wt_dist_P10of10-');
    ImageFileName=strcat(temp, '.png');
    X=getframe(gcf);
    imwrite(X.cdata,ImageFileName,'png');

    %Results
        Validation_struct(Validation_index).AUC=ROCdata.AUC;
        Validation_struct(Validation_index).SE=ROCdata.SE;
        Validation_struct(Validation_index).Sensitivity=ROCdata.Sensitivity;
        Validation_struct(Validation_index).Specificity=ROCdata.Specificity;
        Validation_struct(Validation_index).PPV=ROCdata.PPV;
        Validation_struct(Validation_index).NPV=ROCdata.NPV;
        Validation_struct(Validation_index).wt_corr=wt_R;
        Validation_struct(Validation_index).wt_dist=wt_dist;

%Write result
        temp=strcat(filename,filter_drace,'_det');
        ResultsFileName=strcat(temp, '.csv');
        fileID = fopen(ResultsFileName, 'w');
        fprintf(fileID,'locus,AUC,SE,Sensitivity,Specificity,PPV,NPV,wt_corr, wt_dist\n');
        for k=1:Validation_index
            fprintf(fileID,'%s,',Validation_struct(k).locus);
            fprintf(fileID,'%f,',Validation_struct(k).AUC);
            fprintf(fileID,'%f,',Validation_struct(k).SE);
            fprintf(fileID,'%f,',Validation_struct(k).Sensitivity);
            fprintf(fileID,'%f,',Validation_struct(k).Specificity);
            fprintf(fileID,'%f,',Validation_struct(k).PPV);
            fprintf(fileID,'%f,',Validation_struct(k).NPV);
            fprintf(fileID,'%f,',Validation_struct(k).wt_corr);
            fprintf(fileID,'%f\n',Validation_struct(k).wt_dist);
        end

        fclose(fileID); 