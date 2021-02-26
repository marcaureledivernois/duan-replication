%%%%%%%%%%%%%%%%%%%%%% Duan replication %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is the main file for the intensity model estimation
%   macro        : A two dimension matrix containing the macroeconomics
%                  factors. The first dimension is time, second
%                  dimension is variable.
%   firmspecific: A three dimension matrix containing the firm specific
%                 data. The first dimension of this matrix is time, the
%                 second dimension is variable for the firm and the
%                 third dimension is company.
%   firmlist:     A m by 3 matrix. m: number of firms. First column
%                 is the starting date for the firm to appear in the
%                 dataset, second columun is the exit date. Third
%                 column is the exit type 0: surviving, 1: default, 
%                 2: other exit.


% Marc-Aurèle Divernois
% Comments : 
% 1. Need to add DtD
% 2. need winsorize, not truncate ! replace extremes by the quantile, not
%    delete row, otherwise will mess up firmlist startdate etc
% 4. if firm has a variable with only NaN, maybe remove firm bcz it mess up
%    fisher matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;

% deltas interpol
% load extended_datatensor_deltas.mat;
% firmspecific = datatensor_deltas;

% leveltrend interpol
load extended_datatensor_leveltrend.mat;
firmspecific = datatensor_leveltrend;

load datatensor_lagged_pd.mat   %test to compare with weird algo of predictdrop
firmspecific = datatensor_lagged;



% leveltrend fillmissing
% load extended_datatensor_leveltrend_fillmissing.mat
% firmspecific = datatensor_leveltrend;

load macro_subsample.mat;
load macrodata_pd.mat;    %test to compare with weird algo of predictdrop

load extended_firmlist_datadate.mat;
load firmlist_pd.mat;   %test to compare with weird algo of predictdrop
 
%%



N_param=size(firmspecific,2);
N_horizon=12;
fprintf('Total number of taus: %d (tau = 0 to %d)\n',N_horizon,N_horizon-1');
estimation_all=zeros(2*(N_param+3)+4,N_horizon);
std_all=zeros(2*(N_param+3)+4,N_horizon);
pvalue_all=zeros(2*(N_param+3)+4,N_horizon);

% indexreturn = macro(:,3); put this back
% interestrate = macro(:,2); put this back

macro = macrodata; % remove this when test predictdrop finished
indexreturn = macro(:,1); % remove this when test predictdrop finished
interestrate = macro(:,2); % remove this when test predictdrop finished

%estimation
%to be consistent with the paper, we only add bailout term to default
%intensity function for tau between 0 and 36

testdata = cell(N_horizon,3);
testdata2 = cell(N_horizon,3);
for tau=0:(N_horizon-1)
    fprintf('tau=%d, estimation in progress...',tau);
    tic
    %if tau<12
    %   [estimation1 std1 pvalue1 loglik]=mle_intensity([indexreturn interestrate], firmspecific, firmlist, 1, tau, 1); 
    %   [estimation2 std2 pvalue2 loglik]=mle_intensity([indexreturn interestrate], firmspecific, firmlist, 2, tau, 0); 
    %   estimation_all(:,tau+1)=[estimation1;NaN;NaN;estimation2];
    %   std_all(:,tau+1)=[std1;NaN;NaN;std2];
    %   pvalue_all(:,tau+1)=[pvalue1;NaN;NaN;pvalue2];
    %else

    
       [estimation1 std1 pvalue1 loglik testdata]=mle_intensity([indexreturn interestrate], firmspecific, firmlist, 1, tau, 0,testdata); 
       [estimation2 std2 pvalue2 loglik testdata2]=mle_intensity([indexreturn interestrate], firmspecific, firmlist, 2, tau, 0,testdata2); 
       estimation_all(:,tau+1)=[NaN;NaN;estimation1;NaN;NaN;estimation2];
       std_all(:,tau+1)=[NaN;NaN;std1;NaN;NaN;std2];
       pvalue_all(:,tau+1)=[NaN;NaN;pvalue1;NaN;NaN;pvalue2];
    %end
    toc
end
%output results to a text file
fid=fopen('results_replication_1991.txt','w+');
for i=1:size(estimation_all,1)
    for j=1:N_horizon
        if isnan(estimation_all(i,j))
            fprintf(fid,',');
            if j==N_horizon
                fprintf(fid,'\n');
            end
        else
            fprintf(fid,'%1.3f',estimation_all(i,j));
            if pvalue_all(i,j)<0.01
                fprintf(fid,'***,');
            elseif pvalue_all(i,j)<0.05
                fprintf(fid,'**,');
            elseif pvalue_all(i,j)<0.1
                fprintf(fid,'*,');
            else
                fprintf(fid,',');
            end
            if j==N_horizon
                fprintf(fid,'\n');
            end
        end
    end
    
    for j=1:N_horizon
        if isnan(estimation_all(i,j))
            fprintf(fid,',');
            if j==N_horizon
                fprintf(fid,'\n');
            end
        else
            fprintf(fid,'(%1.3f),',std_all(i,j));
            if j==N_horizon
                fprintf(fid,'\n');
            end
        end
    end
end
fclose(fid);

%% Plots

% Delta
% pars = ["intercept";
% "sp500";
% "treasury";
% "Cash/TA";
% "NI/TA";
% "Size";
% "DtD";
% "MB Ratio";
% "\Delta Cash/TA";
% "\Delta NI/TA";
% "\Delta Size";
% "\Delta DtD";
% "\Delta MB Ratio"];

% LevelTrend
pars = ["intercept";
"sp500";
"treasury";
"Cash/TA (level)";
"NI/TA (level)";
"Size (level)";
"DtD (level)";
"MB Ratio (level)";
"Cash/TA (trend)";
"NI/TA (trend)";
"Size (trend)";
"DtD (trend)";
"MB Ratio (trend)"];

 num = 15;
 x = linspace(0,35,36);
 ii = ceil(num/3);
 jj = ceil(num/ii);

 % alphas -> parameters for default forward intensity
for alphas = 3:15

    
    y1 = estimation_all(alphas,:);
    y2 = estimation_all(alphas,:) + 1.645*std_all(alphas,:);
    y3 = estimation_all(alphas,:) - 1.645*std_all(alphas,:);
    
    subplot(ii,jj,alphas-2)
    plot(x,y1,'b',x,y2,'g-.',x,y3,'r-.', 'LineWidth', 2);
    grid on
    title(sprintf('%s',pars(alphas-2)))
end

 % betas -> parameters for combined exit forward intensity
for betas = 18:30

    
    y1 = estimation_all(betas,:);
    y2 = estimation_all(betas,:) + 1.645*std_all(betas,:);
    y3 = estimation_all(betas,:) - 1.645*std_all(betas,:);
    
    subplot(ii,jj,betas-17)
    plot(x,y1,'b',x,y2,'g-.',x,y3,'r-.', 'LineWidth', 2);
    grid on
    title(sprintf('%s',pars(betas-17)))
end



