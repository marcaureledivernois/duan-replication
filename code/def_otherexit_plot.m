%survdef each year
clear all;clc;

load extended_datatensor_leveltrend.mat;
firmspecific = datatensor_leveltrend;

load macro_subsample.mat;
load extended_firmlist_datadate.mat;


%% plot of default, other exits and surviving

%def
def = firmlist(firmlist(:,4)==1,:);
periodicity = 12;

ndef = [];
ndef(1) = sum(def(:,3) < periodicity);

for i=2:336/periodicity
    ndef(i) = sum(def(:,3) < periodicity*(i) & def(:,3) >= periodicity*(i-1));
end


%other exit
oe = firmlist(firmlist(:,4)==2,:);
noe = []; 
noe(1) = sum(oe(:,3) < periodicity);

for i=2:336/periodicity
    noe(i) = sum(oe(:,3) < periodicity*(i) & oe(:,3) >= periodicity*(i-1));
end


%surviving
nentering =  [];
nentering(1) = sum(firmlist(:,2) <= 12);

for i=2:336/periodicity
    nentering(i) = sum(firmlist(:,2) <= 12*(i) & firmlist(:,2) > 12*(i-1));
end

nsurv = nentering - ndef - noe;
cumsum_nsurv = cumsum(nsurv);

figure
subplot(3,1,1)
stairs((1:336/periodicity) + 1991,ndef)
grid on
title('Total number of defaulted firms per year')

subplot(3,1,2)
stairs((1:336/periodicity) + 1991  ,noe)
grid on
title('Total number of firms exited for other reasons per year')

subplot(3,1,3)
stairs((1:336/periodicity) + 1991  ,cumsum_nsurv)
grid on
title('Total number of firms surviving each year')


