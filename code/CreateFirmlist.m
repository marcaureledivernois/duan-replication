%% FIRMLIST
CUSIP_list = unique(database(:,1));
m = size(CUSIP_list,1);

index_of_first_appearance = [1; find(diff(database{:,1}) ~= 0)+1 ];


start_date = NaN(m,1);
for i = 1:m
    tmp = num2str(database{index_of_first_appearance(i),'datadate'});
    tmp2 = tmp(1:6);
    tmp3 = str2double(tmp2);
    start_date(i,1) = tmp3;
end

std_start_date = NaN(length(start_date),1);
for i = 1:length(start_date)
    std_start_date(i) =  stddate(start_date(i),startyear,monthadjust);
end

std_start_date;

% date of last entry of the firm in the sample

index_of_last_appearance = [find(diff(database{:,1}) ~= 0) ; size(database,1)];


exit_date = NaN(m,1);
for i = 1:m
    tmp = num2str(database{index_of_last_appearance(i),'datadate'});
    tmp2 = tmp(1:6);
    tmp3 = str2double(tmp2);
    exit_date(i,1) = tmp3;
end
    
std_exit_date = NaN(length(exit_date),1);
for i = 1:length(exit_date)
    std_exit_date(i) =  stddate(exit_date(i),startyear,monthadjust);
end

std_exit_date;

% delisting code  : 02 - Bankruptcy 
%                   03 - Liquidation

dlrsn = database{index_of_last_appearance,'dlrsn'};

deletion_reason = NaN(size(dlrsn));

for i = 1:length(deletion_reason)
    if isnan(dlrsn(i))
        deletion_reason(i) = 0;
    elseif dlrsn(i) == 2     
        deletion_reason(i) = 1;
    elseif dlrsn(i) == 3     
        deletion_reason(i) = 1;
    else
        deletion_reason(i) = 2;
    end
end


deletion_reason;
% merge

firmlist = [CUSIP_list{:,1} std_start_date std_exit_date deletion_reason];
save('extended_firmlist_datadate.mat','firmlist');
