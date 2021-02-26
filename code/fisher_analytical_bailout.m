%compute fisher information matrix using the analytical form of gradient
%for the case with bailout term
function fisher=fisher_analytical_bailout(param,macro,firmspecific,firmlist,t_bailout)
delta_t=1/12;
fisher=zeros(size(param,1),size(param,1));
g=zeros(1,size(param,1));
i=1;
while i<=size(firmlist,1)
    starttime=firmlist(i,1);
    endtime=firmlist(i,2);
    exittype=firmlist(i,3);
    len=endtime-starttime+1;
    if size(macro,2) > 0 && size(firmspecific,1) >0 
        covariate=[ones(len,1) macro(starttime:endtime,:) firmspecific(starttime:endtime,:,i)];
    elseif size(macro,2) > 0 
        covariate=[ones(len,1) macro(starttime:endtime,:) endtime-(starttime:endtime)'];
    else
        covariate=[ones(len,1) firmspecific(starttime:endtime,:,i) endtime-(starttime:endtime)'];
    end
    
    if len==1
        temp0=delta_t*exp(-param(1)^2*exp(-param(2)^2*max(starttime-t_bailout,0))*(starttime>t_bailout)+covariate*param(3:end));
        temp1=param(1)^2*exp(-param(2)^2*max(starttime-t_bailout,0))*(starttime>t_bailout);
        temp2=2*param(1)*exp(-param(2)^2*max(starttime-t_bailout,0))*(starttime>t_bailout);
        temp3=2*param(2)*max(starttime-t_bailout,0);
        if exittype==0 || exittype==2
            g(1)=temp0*temp2;
            g(2)=-temp0*temp1*temp3;
            g(3:end)=-temp0*covariate;
        else
            g(1)=-exp(-temp0)*temp0*temp2/(1-exp(-temp0));
            g(2)=exp(-temp0)*temp0*temp1*temp3/(1-exp(-temp0));
            g(3:end)=exp(-temp0)*temp0/(1-exp(-temp0))*covariate;
        end
    else
        temp0=delta_t*exp(-param(1)^2*exp(-param(2)^2*max((starttime:endtime)'-t_bailout,0)).*((starttime:endtime)'>t_bailout)+covariate*param(3:end));
        temp1=param(1)^2*exp(-param(2)^2*max((starttime:endtime)'-t_bailout,0)).*((starttime:endtime)'>t_bailout);
        temp2=2*param(1)*exp(-param(2)^2*max((starttime:endtime)'-t_bailout,0)).*((starttime:endtime)'>t_bailout);
        temp3=2*param(2)*max((starttime:endtime)'-t_bailout,0);
        if exittype==0 || exittype==2
            temp=sparse((1:len)', (1:len)', temp0, len, len);
            g(1)=sum(temp0.*temp2);
            g(2)=-sum(temp0.*temp1.*temp3);
            g(3:end)=sum(-temp*covariate);
        else
            g(1)=sum(temp0(1:end-1).*temp2(1:end-1))-exp(-temp0(end))*temp0(end)*temp2(end)/(1-exp(-temp0(end)));
            g(2)=-sum(temp0(1:end-1).*temp1(1:end-1).*temp3(1:end-1))+exp(-temp0(end))*temp0(end)*temp1(end)*temp3(end)/(1-exp(-temp0(end)));
            temp=sparse((1:(len-1))', (1:(len-1))', temp0(1:end-1), (len-1), (len-1));
            g(3:end)=-sum(temp*covariate(1:end-1,:))+exp(-temp0(end))*temp0(end)/(1-exp(-temp0(end)))*covariate(end,:);
        end
    end
    
    fisher=fisher+g'*g;
    i=i+1;
end

