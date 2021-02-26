%ML estimation of the model with bailout term using analytical gradient and hessian

function [param, f0, h, noerror, noiter]=opt_analytical_bailout(firm0, firm1, firm2, tol, max_iter,initial_value)

%initialization
delta_t=1/12;
size0 = size(firm0,1);
size1 = size(firm1,1);
size2 = size(firm2,1);
epsilon= 1.0e10;
mu = 2.0;
param=zeros(size(firm0,2)+1,1);
param(1:2)=[1;0.2];
param(3:end)=initial_value;
f0 = -sum(exp(-param(1)^2*exp(-param(2)^2*firm0(:,1)).*(firm0(:,1)>0)+firm0(:,2:end)*param(3:end))*delta_t)+sum(log(1-exp(-exp(-param(1)^2*exp(-param(2)^2*firm1(:,1)).*(firm1(:,1)>0)+firm1(:,2:end)*param(3:end))*delta_t)))-sum(exp(-param(1)^2*exp(-param(2)^2*firm2(:,1)).*(firm2(:,1)>0)+firm2(:,2:end)*param(3:end))*delta_t);

f0 = -f0;

noerror = 0;
noiter = 0;
while epsilon>tol && noerror < 1 && noiter < max_iter+1;
    noiter = noiter+1;
    %compute gradient and hessian according to their analytical form, and
    %optimize using the (damped) Newton method
    g=zeros(size(param));
    h=zeros(size(param,1),size(param,1));
    temp0_0=delta_t*exp(-param(1)^2*exp(-param(2)^2*firm0(:,1)).*(firm0(:,1)>0)+firm0(:,2:end)*param(3:end));
    temp0_1=delta_t*exp(-param(1)^2*exp(-param(2)^2*firm1(:,1)).*(firm1(:,1)>0)+firm1(:,2:end)*param(3:end));
    temp0_2=delta_t*exp(-param(1)^2*exp(-param(2)^2*firm2(:,1)).*(firm2(:,1)>0)+firm2(:,2:end)*param(3:end));
    temp1_0=param(1)^2*exp(-param(2)^2*firm0(:,1)).*(firm0(:,1)>0);
    temp1_1=param(1)^2*exp(-param(2)^2*firm1(:,1)).*(firm1(:,1)>0);
    temp1_2=param(1)^2*exp(-param(2)^2*firm2(:,1)).*(firm2(:,1)>0);
    temp2_0=2*param(1)*exp(-param(2)^2*firm0(:,1)).*(firm0(:,1)>0);
    temp2_1=2*param(1)*exp(-param(2)^2*firm1(:,1)).*(firm1(:,1)>0);
    temp2_2=2*param(1)*exp(-param(2)^2*firm2(:,1)).*(firm2(:,1)>0);
    temp3_0=2*param(2)*firm0(:,1);
    temp3_1=2*param(2)*firm1(:,1);
    temp3_2=2*param(2)*firm2(:,1);
    temp4_0=2*exp(-param(2)^2*firm0(:,1)).*(firm0(:,1)>0);
    temp4_1=2*exp(-param(2)^2*firm1(:,1)).*(firm1(:,1)>0);
    temp4_2=2*exp(-param(2)^2*firm2(:,1)).*(firm2(:,1)>0);
    temp5_0=2*firm0(:,1);
    temp5_1=2*firm1(:,1);
    temp5_2=2*firm2(:,1);
    
    %compute gradient(3:end) and hessian(3:end,3:end) first. (the first two
    %parameters are bailout parameters)
    temp0 = sparse((1:size0)', (1:size0)', temp0_0, size0, size0);
    temp1 = temp0_1;
    temp2 = sparse((1:size2)', (1:size2)', temp0_2, size2, size2);
    temp1_g = sparse((1:size1)', (1:size1)', exp(-temp1).*temp1./(1-exp(-temp1)), size1, size1);
    temp1_h = sparse((1:size1)', (1:size1)', (1-temp1-exp(-temp1))./(1-exp(-temp1)), size1, size1);
    g0 = -temp0*firm0(:,2:end);
    g1 = temp1_g*firm1(:,2:end);
    g2 = -temp2*firm2(:,2:end);
    h0 = firm0(:,2:end)'*g0;
    g0 = sum(g0)';
    h1 = firm1(:,2:end)'*temp1_h*g1;
    g1 = sum(g1)';
    h2 = firm2(:,2:end)'*g2;
    g2= sum(g2)';
    h(3:end,3:end) = -(h0+h1+h2);
    g(3:end) = -(g0+g1+g2);
    
    %gradient(1) and gradient(2)
    g_1=sum(temp0_0.*temp2_0)-sum(exp(-temp0_1).*temp0_1.*temp2_1./(1-exp(-temp0_1)))+sum(temp0_2.*temp2_2);
    g_2=-sum(temp0_0.*temp1_0.*temp3_0)+sum(exp(-temp0_1).*temp0_1.*temp1_1.*temp3_1./(1-exp(-temp0_1)))-sum(temp0_2.*temp1_2.*temp3_2);
    g(1)=-g_1;
    g(2)=-g_2;
    
    %hessian(1:2,1:2)
    h_11=sum(temp0_0.*(temp4_0-temp2_0.^2))-sum((exp(-temp0_1).*(temp0_1.^2).*(temp2_1.^2)+(1-exp(-temp0_1)).*exp(-temp0_1).*temp0_1.*(temp4_1-temp2_1.^2))./(1-exp(-temp0_1)).^2)+sum(temp0_2.*(temp4_2-temp2_2.^2));
    h_22=sum(temp0_0.*temp1_0.*(temp3_0.^2.*(1-temp1_0)-temp5_0))-sum((exp(-temp0_1).*(temp0_1.^2).*(temp1_1.^2).*(temp3_1.^2)-(1-exp(-temp0_1)).*(exp(-temp0_1)).*temp0_1.*temp1_1.*((temp1_1-1).*(temp3_1.^2)+temp5_1))./(1-exp(-temp0_1)).^2)+sum(temp0_2.*temp1_2.*(temp3_2.^2.*(1-temp1_2)-temp5_2));
    h_12=sum(temp0_0.*temp2_0.*temp3_0.*(temp1_0-1))+sum((exp(-temp0_1).*(temp0_1).^2.*temp1_1.*temp2_1.*temp3_1-(1-exp(-temp0_1)).*exp(-temp0_1).*temp0_1.*temp2_1.*temp3_1.*(temp1_1-1))./(1-exp(-temp0_1)).^2)+sum(temp0_2.*temp2_2.*temp3_2.*(temp1_2-1));
    h(1,1)=-h_11;
    h(2,2)=-h_22;
    h(1,2)=-h_12;
    h(2,1)=-h_12;
    
    %hessian(1,3:end) and hessian(3:end,1)
    temp0=sparse((1:size0)', (1:size0)', temp0_0.*temp2_0, size0, size0);
    temp1=sparse((1:size1)', (1:size1)', exp(-temp0_1).*temp0_1.*(1-temp0_1-exp(-temp0_1)).*temp2_1./(1-exp(-temp0_1)).^2, size1, size1);
    temp2=sparse((1:size2)', (1:size2)', temp0_2.*temp2_2, size2, size2);
    h_13=sum(temp0*firm0(:,2:end))-sum(temp1*firm1(:,2:end))+sum(temp2*firm2(:,2:end));
    h(1,3:end)=-h_13;
    h(3:end,1)=-h_13';
    
    %hessian(2,3:end) and hessian(3:end,2)
    temp0=sparse((1:size0)', (1:size0)', temp0_0.*temp1_0.*temp3_0, size0, size0);
    temp1=sparse((1:size1)', (1:size1)', exp(-temp0_1).*temp0_1.*(1-temp0_1-exp(-temp0_1)).*temp1_1.*temp3_1./(1-exp(-temp0_1)).^2, size1, size1);
    temp2=sparse((1:size2)', (1:size2)', temp0_2.*temp1_2.*temp3_2, size2, size2);
    h_23=-sum(temp0*firm0(:,2:end))+sum(temp1*firm1(:,2:end))-sum(temp2*firm2(:,2:end));
    h(2,3:end)=-h_23;
    h(3:end,2)=-h_23';
    
    A = h + mu*eye(size(h,1));
    
        
    m=0;
    while m==0;
        try
            chol(A);
            m=1;
        catch
            m=0;
            mu = 2.0*mu;
            A = h + mu*eye(size(h,1));
        end
    end
    
    if det(A) > 1.0e-15;
        d = -inv(A)*g;
    else
        d = -g/norm(g);
    end
    param1=param+d;
    

    f1 = -sum(exp(-param1(1)^2*exp(-param1(2)^2*firm0(:,1)).*(firm0(:,1)>0)+firm0(:,2:end)*param1(3:end))*delta_t)+sum(log(1-exp(-exp(-param1(1)^2*exp(-param1(2)^2*firm1(:,1)).*(firm1(:,1)>0)+firm1(:,2:end)*param1(3:end))*delta_t)))-sum(exp(-param1(1)^2*exp(-param1(2)^2*firm2(:,1)).*(firm2(:,1)>0)+firm2(:,2:end)*param1(3:end))*delta_t);

    f1 = -f1;
           
    if isnan(f1);
        noerror = 1;
    else
        param = param1;
        mu = mu / 2;
        epsilon = g'*g;
        f0 = f1;
    end
    
end