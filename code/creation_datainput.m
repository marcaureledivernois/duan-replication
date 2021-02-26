%%%%%%%%%%%%%%%%%%%%%% Duan replication %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code creates the input data matrices that will be used in the code
%   macro        : A two dimension matrix containing the macroeconomics
%                  factors. The first dimension is time, second
%                  dimension is variable.
%   firmspecific: A three dimension matrix containing the firm specific
%                 data. The first dimension of this matrix is time, the
%                 second dimension is variable for the firm and the
%                 third dimension is company.
%   firmlist:     A m by 4 matrix. m: number of firms. First column
%                 is the CUSIP. 2nd col is starting date for the firm to 
%                 appear in the dataset, third columun is the exit date.
%                 fourth column is the exit type 0: surviving, 1: default, 
%                 2: other exit.
% 
% First date : 2009Q4 (datacqtr) // 201001 (datadate)
% First date :                   // 196103 (datadate)  -- EXTENDED
% Last date  : 2018Q3 (datacqtr) // 201809 (datadate)
% Last date  :                   // 201812 (datadate)


% Marc-Aurèle Divernois
% Comments : 
% 1. need to know with which date work (datacqtr, datadate,..?)
%    for this, need to figure out out is dldate defined. also, for better
%    results, "cheat" a bit and take the smallest interval between dldate
%    and publication... otherwise will never be able to predict default
%    ---> need datadate i think bcz monthly frequency, and looking at duan
%    code, he adjust data for interval and only take macro var inbetween
% 2. took into account liquid&bankruptcy as default
% 3. need to check delisting date (closer to datadate or datacqtr?)
% 4. check comments in code below!
% 5. try estimation with fillmissing instead of interpol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% macro.mat                     : macro variables from 19610331 to 201812
%
% firmlist.mat                  : firmlist with a quarterly base (datacqtr)
% firmlist_datadate.mat         : firmlist with a monthly base (datadate)
% extended_firmlist_datadate.mat: firmlist with a monthly base (datadate)
%                                 for the extended dataset
%
% extended_firmspecific_raw.mat : raw extended dataset (before fillmissing)
% extended_firmspecific_raw.mat : raw extended dataset (with NaN)
% extended_datatensor.mat       : tensor before creation of trend/level
% extended_datatensor_deltas.mat: tensor with diff
% extended_database.mat         : raw csv database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
%% parameters

overwrite = 0;
SubsampleStart = 19910101;         %19610331
%%

%extendeddatabase : database from 19610331 to 2018
%database : database from 2010 to 2018

load('extended_database.mat');
database = extendeddatabase;

startyear = num2str(SubsampleStart); % 1961
startyear = str2double(startyear(1:4));
monthadjust = num2str(SubsampleStart);  %2        % because month of mindatadate is march
monthadjust = str2double(monthadjust(5:6))-1;

database = database(database{:,'datadate'}>=SubsampleStart,:);

CUSIP_list = unique(database(:,1));
m = size(CUSIP_list,1);

mindatadate = min(database{:,'datadate'}); %19610331
maxdatadate = max(database{:,'datadate'}); %20181231


%% FIRMLIST

CreateFirmlist;

%% FIRMSPECIFIC RAW
%Elapsed time is 550.161725 seconds.

tic
Nperiods = stddate(201812,startyear,monthadjust);
vars = {'gvkey','datadate','atq','cheq','niy','prccq','cshoq','dlttq','lctq','ltq'};
Nvar = size(vars,2);
Nfirms = m;

firmspecific = NaN(Nperiods,Nvar,Nfirms);


for i = 1:m
    tmpdata = database(index_of_first_appearance(i):index_of_last_appearance(i),vars);
        for j = 1:size(tmpdata,1)
            tmp = num2str(tmpdata{j,'datadate'});
            tmp2 = tmp(1:6);
            tmp3 = str2double(tmp2);
            t = stddate(tmp3,startyear,monthadjust);
            firmspecific(t,:,i) = tmpdata{j,vars};
        end
    clear tmpdata
end
toc


if(overwrite == 1)
    save('extended_firmspecific_raw.mat','firmspecific');
    load('extended_firmspecific_raw.mat');
end

%% CREATE VARIABLES
load('extended_firmspecific_raw.mat');

% Comments : 
% 2. maybe very old data are not good...  -> ok because work with ratios
% 3. for ML later : maybe need to do some standardization year by year
% 4. this variable is probably significative : Sign of net income
% 5. need to do trend and level -> i did diff_1
% 7. maybe lag variables
% 9. check variables (and also summary stats) to see if there is nothing
%    impossible... (currently i have some mb ratio negative ..) + dont
%    forget also to check delta sigma etc


% fill gaps
for i = 1:size(firmspecific,3)
    % replace 0 by NaN for variables that can't take 0 (at,cshoq,...)
    % so later it will be fillmissing by last value
    firmspecific(firmspecific(:,3,i)==0,3,i) = NaN;   %atq
    firmspecific(firmspecific(:,6,i)==0,6,i) = NaN;   %prccq
    firmspecific(firmspecific(:,7,i)==0,7,i) = NaN;   %cshoq
    
    % replace by NaN variables that can't be negative (lctq,ltq,...)
    firmspecific(firmspecific(:,9,i)<0,9,i) = NaN;   %lctq
    firmspecific(firmspecific(:,10,i)<0,10,i) = NaN;   %ltq
    firmspecific(firmspecific(:,8,i)<0,8,i) = NaN;   %dlttq
    
    % replace all NaN by last available value
     firmspecific(firmlist(i,2):firmlist(i,3),:,i) = ... 
       fillmissing(firmspecific(firmlist(i,2):firmlist(i,3),:,i),'previous',1);
    
    % replace all NaN by interpolated value
    % for j = 3:size(firmspecific,2)
    %    t = linspace(0.1, 10, numel(firmspecific(:,j,i)));
    %    nans = isnan(firmspecific(:,j,i));
    %    % replace all NaNs in x with linearly interpolated values
    %    if sum(nans) < size(firmspecific,1)-1
    %        firmspecific(nans,j,i) = interp1(t(~nans), firmspecific(~nans,j,i), t(nans));
    %    end
    %    clear t nans
    % end
end


% create variables

firmspecific(:,11,:) = firmspecific(:,4,:)./firmspecific(:,3,:);  %cash/ta
firmspecific(:,12,:) = firmspecific(:,5,:)./firmspecific(:,3,:);  %ni/ta
firmspecific(:,13,:) = firmspecific(:,6,:).*firmspecific(:,7,:);  %marketEquity

firmspecific(:,14,:) = log(firmspecific(:,13,:)./nanmean(firmspecific(:,13,:),3));  %size

% firmspecific(:,1,:) = gvkey
% firmspecific(:,2,:) = datadate
% firmspecific(:,3,:) = atq
% firmspecific(:,4,:) = cheq
% firmspecific(:,5,:) = niy
% firmspecific(:,6,:) = prccq
% firmspecific(:,7,:) = cshoq
% firmspecific(:,8,:) = dlttq
% firmspecific(:,9,:) = lctq
% firmspecific(:,10,:) = ltq
% firmspecific(:,11,:) = cash/ta
% firmspecific(:,12,:) = ni/ta
% firmspecific(:,13,:) = marketEquity
% firmspecific(:,14,:) = size
% firmspecific(:,15,:) = market_value_asset;
% firmspecific(:,16,:) = delta;
% firmspecific(:,17,:) = DtD;
% firmspecific(:,18,:) = mu;
% firmspecific(:,19,:) = sigma;
% firmspecific(:,20,:) = MBratio;


if(overwrite == 1)
%     save('extended_firmspecific_pre_dtd.mat','firmspecific');
%     load('extended_firmspecific_pre_dtd.mat');
    save('extended_firmspecific_pre_dtd_fillmissing.mat','firmspecific');
    load('extended_firmspecific_pre_dtd_fillmissing.mat');
end

%% DtD Estiation

DtD_estimation;

%% Extract only variables needed for Duan

%load('extended_firmspecific_with_dtd.mat');
load('extended_firmspecific_with_dtd_fillmissing.mat');

% create mb ratio
firmspecific(:,20,:) = firmspecific(:,15,:)./firmspecific(:,3,:);

datatensor = firmspecific(:,[11 12 14 17 20],:);

if(overwrite == 1)
    save('extended_datatensor_fillmissing.mat','datatensor');
    load('extended_datatensor_fillmissing.mat');
    
    %save('extended_datatensor.mat','datatensor');
    %load('extended_datatensor.mat');
end

% 5. need to do trend and level

% maybe better do compute t - t-1 , t - t-2, ...
% problem with Ret_1 : vector takes the form [0 0 r1 0 0 r2 0 0 r3...]
% maybe fillmissing? or first compute mean year by year of all variables...
% this way i get firm-month data that changes trough time

% datatensor(:,1,:) = CASH/TA
% datatensor(:,2,:) = NI/TA
% datatensor(:,3,:) = Size
% datatensor(:,4,:) = DtD
% datatensor(:,5,:) = MBratio
% datatensor(:,6,:) = diff_1 CASH/TA
% datatensor(:,7,:) = diff_1 NI/TA
% datatensor(:,8,:) = diff_1 Size
% datatensor(:,9,:) = diff_1 DtD
% datatensor(:,10,:) = diff_1 MBratio

% Winsorize

% first need to transform 3d into 2d

nbline = size(datatensor,1);
datatensor_trans = NaN(size(datatensor,1)*size(datatensor,3),size(datatensor,2));
for i = 1:size(datatensor,3)
    datatensor_trans((i-1)*nbline+1:i*nbline,:) = datatensor(:,:,i);
end


for i=1:size(datatensor_trans,2)
    [datatensor_trans(:,i),indx] = winsor(datatensor_trans(:,i),[5 95]);
end

for i = 1:size(datatensor,3)
    datatensor(:,:,i) = datatensor_trans((i-1)*nbline+1:...
    i*nbline,:);
end


m = size(datatensor,3);

if(deltas == 1)
    datatensor(:,6,:) = [NaN(1,1,m) ; diff(datatensor(:,1,:))];
    datatensor(:,7,:) = [NaN(1,1,m) ; diff(datatensor(:,2,:))];
    datatensor(:,8,:) = [NaN(1,1,m) ; diff(datatensor(:,3,:))];
    datatensor(:,9,:) = [NaN(1,1,m) ; diff(datatensor(:,4,:))];
    datatensor(:,10,:) = [NaN(1,1,m) ; diff(datatensor(:,5,:))];

    datatensor_deltas = datatensor;

    save('extended_datatensor_deltas.mat','datatensor_deltas');
    load('extended_datatensor_deltas.mat');
end

if(leveltrend == 1)
    windowsize = 11;  % 12 months
    datatensor_leveltrend = [movmean(datatensor,[windowsize 0]) datatensor-movmean(datatensor,[windowsize 0])];

    save('extended_datatensor_leveltrend_fillmissing.mat','datatensor_leveltrend');
    load('extended_datatensor_leveltrend_fillmissing.mat');
    
    %save('extended_datatensor_leveltrend.mat','datatensor_leveltrend');
    %load('extended_datatensor_leveltrend.mat');
end



%% MACRO
if(overwrite == 1)
    load('macro.mat','macro');
    macro = macro(macro(:,1)>=SubsampleStart,:);
    
    save('macro_subsample.mat','macro');
end












