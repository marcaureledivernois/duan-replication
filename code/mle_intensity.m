
%Purpose: estimate the parameters for the intensity function and generate
%the standard error for the estimation.

%Input: macro       : A two dimension matrix containing the macroeconomics
%                     factors. The first dimension is time, second
%                     dimension is variable.
%       firmspecific: A three dimension matrix containing the firm specific
%                     data. The first dimension of this matrix is time, the
%                     second dimension is variable for the firm and the
%                     third dimension is company.
%       firmlist:     A m by 3 matrix. m: number of firms. First column
%                     is the starting date for the firm to appear in the
%                     dataset, second columun is the exit date. Third
%                     column is the exit type 0: surviving, 1: default, 2: other exit.
%       intensitytype:type of intensity to be estimated 1: default, 2:
%                     other exit
%       tau:          horizon, same as the tau in the paper
%       bailout:      whether to include the bailout term, 1: yes, 0: no

%Output: estimation:  estimation of the parameters in the intensity
%                     function
%        stder:       standard error 
%        pvalue:      p-value
%        loglik:      pseudo log-likelihood

function [estimation stder pvalue loglik testdata]=mle_intensity(macro, firmspecific, firmlist, intensitytype, tau, bailout,testdata)

%initialization
tol=1.0e-12;
max_iter=50;
t_bailout=212;
traintest=true;

N_macro=size(macro,2);
N_firmspecific=size(firmspecific,2);
NPAR=N_macro+N_firmspecific+1;

%adjust firmspecific and firmlist
firmlist(:,3)=firmlist(:,3)-tau;
firmspecific=firmspecific(:,:,firmlist(:,3)>=firmlist(:,2));  %keep only firms where t_0 =< t_c
firmlist=firmlist(firmlist(:,3)>=firmlist(:,2),:);            %keep only firms where t_0 =< t_c
nobs=size(firmlist,1);                                        %nobs = number of firms
%generate the input for likelihood function
k=0;
firmdata=zeros(sum(firmlist(:,3)-firmlist(:,2)+1),NPAR+3);
%firmdata=zeros(sum(firmlist(:,3)-firmlist(:,2)+1),NPAR+2);    % 3 more parameters : time, intensity type and constant

%firmdata(:,1) = bailout
%firmdata(:,2) = constant
%firmdata(:,3) = SPreturn
%firmdata(:,4) = TBills
%firmdata(:,5,:) = CASH/TA
%firmdata(:,6,:) = NI/TA
%firmdata(:,7,:) = Size
%firmdata(:,8,:) = DtD
%firmdata(:,9,:) = MBratio
%firmdata(:,10,:) = Ret_1 CASH/TA
%firmdata(:,11,:) = Ret_1 NI/TA
%firmdata(:,12,:) = Ret_1 Size
%firmdata(:,13,:) = Ret_1 DtD
%firmdata(:,14,:) = Ret_1 MBratio
%firmdata(:,15,:) = periodsleft
%firmdata(:,16,:) = exittype


for i=1:nobs
    starttime=firmlist(i,2);
    endtime=firmlist(i,3);
    exittype=firmlist(i,4);
    len=endtime-starttime+1;
    if size(macro,2) > 0 && size(firmspecific,1) >0 
        firmdata(k+1:k+len,:)=[max((starttime:endtime)'-t_bailout,0) ones(len,1) macro(starttime:endtime,:) firmspecific(starttime:endtime,:,i) endtime-(starttime:endtime)' exittype*ones(len,1)];
    elseif size(macro,2) > 0 
        firmdata(k+1:k+len,:)=[max((starttime:endtime)'-t_bailout,0) ones(len,1) macro(starttime:endtime,:) endtime-(starttime:endtime)' exittype*ones(len,1)];
    else
        firmdata(k+1:k+len,:)=[max((starttime:endtime)'-t_bailout,0) ones(len,1) firmspecific(starttime:endtime,:,i) endtime-(starttime:endtime)' exittype*ones(len,1)];
    end
    k=k+len;
end
firmdata=firmdata(1:k,:);


index0=find(1-(1-(firmdata(:,end)==0)).*(1-firmdata(:,end-1)>0)); %surviving
firm0=firmdata(index0,1:end-2);
index1=find((firmdata(:,end)==1).*(firmdata(:,end-1)==0));%default
firm1=firmdata(index1,1:end-2);
index2=find((firmdata(:,end)==2).*(firmdata(:,end-1)==0));%other exit
firm2=firmdata(index2,1:end-2);

% removing all lines that contain at least one NaN (missing value)
firm0(any(isnan(firm0), 2), :) = [];
firm1(any(isnan(firm1), 2), :) = [];
firm2(any(isnan(firm2), 2), :) = [];

firm0(any(isinf(firm0), 2), :) = [];
firm1(any(isinf(firm1), 2), :) = [];
firm2(any(isinf(firm2), 2), :) = [];

out_firm0 = firm0(:,2:end);
out_firm1 = firm1(:,2:end);
out_firm2 = firm2(:,2:end);

traintest = true;
% train/test split
if(traintest == true)
    p = 0.8;      % proportion of rows to select for training
    N = size(out_firm0,1);  % total number of rows 
    tf = false(N,1);    % create logical index vector
    tf(1:round(p*N)) = true ;    
    tf = tf(randperm(N));   % randomise order
    firm0_train = out_firm0(tf,:); 
    firm0_test = out_firm0(~tf,:);

    N = size(out_firm1,1);  % total number of rows 
    tf = false(N,1);    % create logical index vector
    tf(1:round(p*N)) = true ;    
    tf = tf(randperm(N));   % randomise order
    firm1_train = out_firm1(tf,:); 
    firm1_test = out_firm1(~tf,:);

    N = size(out_firm2,1);  % total number of rows 
    tf = false(N,1);    % create logical index vector
    tf(1:round(p*N)) = true ;    
    tf = tf(randperm(N));   % randomise order
    firm2_train = out_firm2(tf,:); 
    firm2_test = out_firm2(~tf,:);

    testdata{tau+1,1} = firm0_test;
    testdata{tau+1,2} = firm1_test;
    testdata{tau+1,3} = firm2_test;
    
%Optimization and estimation using analytical gradient and hessian
    if bailout==1
        [param, f0, h, noerror, noiter]=opt_analytical(firm0(:,2:end), firm1(:,2:end), firm2(:,2:end), intensitytype, tol, max_iter);
        [param, f0, h, noerror, noiter]=opt_analytical_bailout(firm0, firm1, firm2, tol, max_iter, param);
    else
        [param, f0, h, noerror, noiter]=opt_analytical(firm0_train, firm1_train, firm2_train, intensitytype, tol, max_iter);
    end
    hess = -h;

else
    %Optimization and estimation using analytical gradient and hessian
    if bailout==1
       [param, f0, h, noerror, noiter]=opt_analytical(firm0(:,2:end), firm1(:,2:end), firm2(:,2:end), intensitytype, tol, max_iter);
       [param, f0, h, noerror, noiter]=opt_analytical_bailout(firm0, firm1, firm2, tol, max_iter, param);
    else
       [param, f0, h, noerror, noiter]=opt_analytical(firm0(:,2:end), firm1(:,2:end), firm2(:,2:end), intensitytype, tol, max_iter);
    end
    hess = -h;
end

%calculate standard error and p-value
if bailout==1
    fisher=fisher_analytical_bailout(param,macro,firmspecific,firmlist,t_bailout);   % need to check in this function. somethings fishy with my data... probably -inf again
else
    fisher=fisher_analytical(param,macro,firmspecific,firmlist,intensitytype);
end
SIG=hess\fisher/hess;
stder=sqrt(diag(SIG));
estimation=param;
loglik=-f0;
%adjust the bailout parameters to be consistent with the paper
if bailout==1
    stder(1)=2*abs(estimation(1))*stder(1);
    stder(2)=2*abs(estimation(2))*stder(2);
    estimation(1)=-estimation(1)^2;
    estimation(2)=estimation(2)^2;
end
pvalue=2*(1-normcdf(abs(estimation./stder)));
if noerror==1 || noiter>max_iter
    display('Warning: the program did not converge, please further check the data and codes!')
end
end