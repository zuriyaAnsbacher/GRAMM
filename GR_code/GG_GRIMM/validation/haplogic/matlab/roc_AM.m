function ROCdata=roc(x,co, alpha)
% %
% Input: x - This is the data matrix. The first column is the column of the data value;
%            The second column is the column of the tag: unhealthy (1) and
%            healthy (0).
%          alpha - significance level (default 0.05)
%           co - threshold value
% Output: The ROC plot;
%         The Area under the curve with Standard error and Confidence
%         interval and comment.
%         Cut-off point for best sensitivity and specificity.
%         (Optional) the test performances at cut-off point.
%

clear args default nu
%Abeer
verbose=1;

tr=repmat('-',1,80);
lu=length(x(x(:,2)==1)); %number of unhealthy subjects
lh=length(x(x(:,2)==0)); %number of healthy subjects
z=sortrows(x,1);
%find unique values in z
labels=unique(z(:,1));
ll=length(labels); %count unique value
a=zeros(ll,2); %array preallocation
ubar=mean(x(x(:,2)==1),1); %unhealthy mean value
hbar=mean(x(x(:,2)==0),1); %healthy mean value
for K=1:ll
    if hbar<ubar
        TP=length(x(x(:,2)==1 & x(:,1)>labels(K)));
        FP=length(x(x(:,2)==0 & x(:,1)>labels(K)));
        FN=length(x(x(:,2)==1 & x(:,1)<=labels(K)));
        TN=length(x(x(:,2)==0 & x(:,1)<=labels(K)));
    else
        TP=length(x(x(:,2)==1 & x(:,1)<labels(K)));
        FP=length(x(x(:,2)==0 & x(:,1)<labels(K)));
        FN=length(x(x(:,2)==1 & x(:,1)>=labels(K)));
        TN=length(x(x(:,2)==0 & x(:,1)>=labels(K)));
    end
    a(K,:)=[TP/(TP+FN) TN/(TN+FP)]; %Sensitivity and Specificity
end

if hbar<ubar
    xroc=flipud([1; 1-a(:,2); 0]); yroc=flipud([1; a(:,1); 0]); %ROC points
    labels=flipud(labels);
else
    xroc=[0; 1-a(:,2); 1]; yroc=[0; a(:,1); 1]; %ROC points
end

Area=trapz(xroc,yroc); %estimate the area under the curve
%standard error of area
Area2=Area^2; Q1=Area/(2-Area); Q2=2*Area2/(1+Area);
V=(Area*(1-Area)+(lu-1)*(Q1-Area2)+(lh-1)*(Q2-Area2))/(lu*lh);
Serror=realsqrt(V);
if verbose
    %confidence interval
    cv=realsqrt(2)*erfcinv(alpha);
    ci=Area+[-1 1].*(cv*Serror);
    if ci(2)>1
        ci(2)=1;
    end
    %z-test
    SAUC=(Area-0.5)/Serror; %standardized area
    p=1-0.5*erfc(-SAUC/realsqrt(2)); %p-value
    %Performance of the classifier
    if Area==1
        str='Perfect test';
    elseif Area>=0.90 && Area<1
        str='Excellent test';
    elseif Area>=0.80 && Area<0.90
        str='Good test';
    elseif Area>=0.70 && Area<0.80
        str='Fair test';
    elseif Area>=0.60 && Area<0.70
        str='Poor test';
    elseif Area>=0.50 && Area<0.60
        str='Fail test';
    else
        str='Failed test - less than chance';
    end
    table=[labels'; yroc(2:end-1)'; 1-xroc(2:end-1)';]';
    
    %if partest.m was downloaded
    if p<=alpha
        %the best cut-off point is the closest point to (0,1)
        d=realsqrt(xroc.^2+(1-yroc).^2); %apply the Pitagora's theorem
        [~,J]=min(d); %find the least distance
%##Abeer commented out to set threshold as an argument         co=labels(J-1); %Set the cut-off point
               
        %table at cut-off point
        if hbar<ubar
            TP=length(x(x(:,2)==1 & x(:,1)>co));
            FP=length(x(x(:,2)==0 & x(:,1)>co));
            FN=length(x(x(:,2)==1 & x(:,1)<=co));
            TN=length(x(x(:,2)==0 & x(:,1)<=co));
        else
            TP=length(x(x(:,2)==1 & x(:,1)<co));
            FP=length(x(x(:,2)==0 & x(:,1)<co));
            FN=length(x(x(:,2)==1 & x(:,1)>=co));
            TN=length(x(x(:,2)==0 & x(:,1)>=co));
        end
        cotable=[TP FP; FN TN];
        [SS, PNp]=partest(cotable);
    end
end
if nargout
    ROCdata.AUC=Area;
    ROCdata.Cutoff=co;
    ROCdata.SE=Serror;
    ROCdata.p_value=p;
    if exist('SS','var')
        ROCdata.Sensitivity=SS(1);
        ROCdata.Specificity=SS(2);
    else
        ROCdata.Sensitivity=0;
        ROCdata.Specificity=0;
    end
    if exist('PNp','var')
        ROCdata.PPV=PNp(1);
        ROCdata.NPV=PNp(2);
    else
        ROCdata.PPV=0;
        ROCdata.NPV=0;
    end
end
