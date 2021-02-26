% This code estimate the Distance to default for the Duan's method replication

%% dtd estimation
% Elapsed time is 5999.028740 seconds.
clear; clc;
tic
load('extended_firmspecific_pre_dtd_fillmissing.mat');
load macro_subsample.mat;
load extended_firmlist_datadate.mat;

% create space for future variables to add 
%(DtD, Market_value_Asset,delta,mu,sigma)

firmspecific(:,15:19,:) = NaN(size(firmspecific,1),5,size(firmspecific,3));

%solve for firms individually to have individual values of delta
for i = 1:size(firmspecific,3)

    
starttime=firmlist(i,2);
endtime=firmlist(i,3);
exittype=firmlist(i,4);
len=endtime-starttime+1;

data=[firmspecific(starttime:endtime,:,i) macro(starttime:endtime,2)];

nbnan = size(data(any(isnan(data(:,[3 8 9 10 13])), 2), :),1);
data(any(isnan(data(:,[3 8 9 10 13])), 2), :) = [];

LongTermDebt = data(:,8);  %dlttq
CurrentLiabilities = data(:,9);  %lctq
LiabilitiesTotal = data(:,10); %ltq
LiabilitiesOther = max(LiabilitiesTotal-CurrentLiabilities-LongTermDebt,0);  
AssetsTotal = data(:,3);  %atq
r = mean(data(:,20))/100;
r = data(:,20)/100;
F = CurrentLiabilities+0.5*LongTermDebt;
market_value_equity = data(:,13);

if size(data,1)>0
    options = optimset('Display','off','MaxFunEvals',50000,'MaxIter',50000);
    LB = [-1 0 0];
    UB = [Inf Inf 1];
    x0 = [0 0.25 0.5];
    [x,fval] = fminsearchbnd(@log_likelihood_fun,x0,LB,UB,options,F, LiabilitiesOther, r, market_value_equity, AssetsTotal); % Call solver

    xL = 0 * ones(1,1);
    xU =  9999999 * ones(1,1);
    D = F +  x(3)*LiabilitiesOther;

    % warning off
    V =solve_dichotomy(@afunction,xL,xU,0.0001,D,r,market_value_equity,x(2));
    mu = x(1) ;
    sigma = x(2);
    delta = x(3) ;

    DTD = [NaN(nbnan,1) ;(log(V./D))./sigma];
    V = [NaN(nbnan,1); V];
    
    firmspecific(starttime:endtime,15:19,i) = [V delta*ones(len,1) DTD mu*ones(len,1)...
    sigma*ones(len,1)];
end


clear data datas x V DTD
end


% change Inf into NaN

c = find(firmspecific>9999999999999);
[a,b,d] = ind2sub(size(firmspecific),c);

for i = 1:size(a,1)
    firmspecific(a(i),b(i),d(i)) = NaN;
end

toc

if(overwrite == 1)
    save('extended_firmspecific_with_dtd_fillmissing.mat','firmspecific');
    % save('extended_firmspecific_with_dtd.mat','firmspecific');
end


%% graphs showing DTD
% 
% %FIRMspecific
% 
% [defindex,~] = find(firmlist(:,4)== 1) ;   %only def firms
% 
% for j = 1:size(defindex,1)
%     i = defindex(j);
%     label = firmlist(i,4);
%     if label == 0
%         color = 'red';
%     elseif label == 1
%         color = 'blue';
%     elseif label == 2
%         color = 'green';
%     end
%     
%     if any(firmspecific(:,17,i)>100)
%         bug = [bug i]
%     end
%     
%     plot(linspace(1,size(firmspecific(:,17,i),1),size(firmspecific(:,17,i),1)),firmspecific(:,17,i),color)
%     hold on
% end
% hold off
% 
% 
% % CONCLUSION : SOMETIMES EXTREME VALUES... DUE TO VERY LOW DEBT AND HIGH MARKET VALUE ASSET!
% % -> MAYBE GIVE DEF POINT INSTEAD?
% 
% 
% 
% % DATATENSOR
% 
% [defindex,~] = find(firmlist(:,4)== 1) ;   %only def firms
% 
% for i = 1:25
%     %i = defindex(j);
%     label = firmlist(i,4);
%     if label == 0
%         color = 'red';
%     elseif label == 1
%         color = 'blue';
%     elseif label == 2
%         color = 'green';
%     end
%     
%     if any(datatensor_deltas(:,4,i)>100)
%         bug = [bug i]
%     end
%     
%     plot(linspace(1,size(datatensor_deltas(:,4,i),1),size(datatensor_deltas(:,4,i),1)),datatensor_deltas(:,4,i),color)
%     hold on
% end
% hold off
% 
% 
% % CONCLUSION : SOMETIMES EXTREME VALUES... DUE TO VERY LOW DEBT AND HIGH MARKET VALUE ASSET!
% % -> MAYBE GIVE DEF POINT INSTEAD?
