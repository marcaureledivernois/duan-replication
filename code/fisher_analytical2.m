%compute fisher information matrix using the analytical form of gradient
%for the case without bailout term
function fisher=fisher_analytical2(param,firmdata,firmlist,intensitytype)
delta_t=1/12;
fisher=zeros(size(param,1),size(param,1));
i=1;
k=1;
while i<=size(firmlist,1)
    starttime=firmlist(i,2);
    endtime=firmlist(i,3);
    exittype=firmlist(i,4);
    len=endtime-starttime+1;
    
    covariate= firmdata(k:k+len-1,2:14,i); % need remove some var to coincide with covariate of fisher_analytical.m
    k = k + len;

    % winsor
    for i=5:14
        [covariate(:,i),indx] = winsor(firmdata(:,i),[5 95]);
    end
    
    
    if intensitytype==1
       if len==1
           if exittype==0 || exittype==2
               g=-delta_t*exp(covariate*param)*covariate;
           else
               temp=delta_t*exp(covariate*param);
               g=exp(-temp)*temp/(1-exp(-temp))*covariate;
           end
       else
           if exittype==0 || exittype==2
               temp=sparse((1:len)', (1:len)', delta_t * exp(covariate*param), len, len);
               g=nansum(-temp*covariate);
           else
               temp0=sparse((1:(len-1))', (1:(len-1))', delta_t * exp(covariate(1:end-1,:)*param), (len-1), (len-1));
               temp1=delta_t*exp(covariate(end,:)*param);
               temp1=exp(-temp1)*temp1/(1-exp(-temp1));
               g=nansum(-temp0*covariate(1:end-1,:))+temp1*covariate(end,:);
           end
       end
    else
        if len==1
            if exittype==0
                g=-delta_t*exp(covariate*param)*covariate;
            elseif exittype==1
                g=zeros(1,size(param,1));
            else
                temp=delta_t*exp(covariate*param);
                g=exp(-temp)*temp/(1-exp(-temp))*covariate;
            end
            if isnan(g)
               g=0;
            end
        else
            if exittype==0
                temp=sparse((1:len)', (1:len)', delta_t * exp(covariate*param), len, len);
                g=nansum(-temp*covariate);
            elseif exittype==1
                temp=sparse((1:(len-1))', (1:(len-1))', delta_t * exp(covariate(1:end-1,:)*param), (len-1), (len-1));
                g=nansum(-temp*covariate(1:end-1,:));
            else
                temp0=sparse((1:(len-1))', (1:(len-1))', delta_t * exp(covariate(1:end-1,:)*param), (len-1), (len-1));
                temp1=delta_t*exp(covariate(end,:)*param);
                temp1=exp(-temp1)*temp1/(1-exp(-temp1));
                g=nansum(-temp0*covariate(1:end-1,:))+temp1*covariate(end,:);
            end
        end
    end
    
%     if sum(isnan(g))>0      %maybe remove firm in this case... may mess up fisher matrix
%         g=0;
%     end
    
    fisher=fisher+g'*g;
    if (sum(any(isnan(fisher))) > 0) || sum(any(isinf(fisher))) > 0
        fprintf("NaN or inf entered into fisher matrix at index %d", i);
        break;
    end
    i=i+1;
end

