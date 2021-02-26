%%% function to map unstandardized datadate (ex : 201008) to 
%%% a standardized date (ex: 8).

%%% potential problem : what should i take as datadate? fiscal year?
%%% calendar year ?... need to match corretcly firms ! ... maybe some work
%%% to do here :-)

%monthadj : takes into account month of starting 

function std_date = stddate(unstand_date,startyear,monthadjust)
str = num2str(unstand_date);
yearstr = str(1:4);
monthstr = str(5:6);
yeardbl = str2double(yearstr);
monthdbl = str2double(monthstr);
std_date = (yeardbl-startyear)*12+monthdbl-monthadjust;
end

