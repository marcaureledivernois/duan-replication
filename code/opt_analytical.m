%ML estimation of the model without bailout term using analytical gradient and hessian

function [param, f0, h, noerror, noiter]=opt_analytical(firm0, firm1, firm2, intensitytype,tol, max_iter)

%initialization
delta_t=1/12;
size0 = size(firm0,1);
size1 = size(firm1,1);
size2 = size(firm2,1);
epsilon= 1.0e10;
mu = 2.0;
param=zeros(size(firm0,2),1);
if intensitytype==1
    f0 = -sum(exp(firm0*param)*delta_t)+sum(log(1-exp(-exp(firm1*param)*delta_t)))-sum(exp(firm2*param)*delta_t);
else
    f0 = -sum(exp(firm0*param)*delta_t)+sum(log(1-exp(-exp(firm2*param)*delta_t)));
end
f0 = -f0;

noerror = 0;
noiter = 0;
while epsilon>tol && noerror < 1 && noiter < max_iter+1;
    noiter = noiter+1;
    %compute gradient and hessian according to their analytical form, and
    %optimize using the (damped) Newton method
    if intensitytype==1 ;
        
        temp0 = sparse((1:size0)', (1:size0)', delta_t * exp(firm0*param), size0, size0);
        temp1 = delta_t*exp(firm1*param);
        temp2 = sparse((1:size2)', (1:size2)', delta_t * exp(firm2*param), size2, size2);
        temp1_g = sparse((1:size1)', (1:size1)', exp(-temp1).*temp1./(1-exp(-temp1)), size1, size1);
        temp1_h = sparse((1:size1)', (1:size1)', (1-temp1-exp(-temp1))./(1-exp(-temp1)), size1, size1);
        g0 = -temp0*firm0;
        g1 = temp1_g*firm1;
        g2 = -temp2*firm2;
        h0 = firm0'*g0;
        g0 = sum(g0)';
        h1 = firm1'*temp1_h*g1;
        g1 = sum(g1)';
        h2 = firm2'*g2;
        g2= sum(g2)';
        h = -(h0+h1+h2);
        g = -(g0+g1+g2);
        
        A = h + mu*eye(size(h,1));
    elseif intensitytype ==2 ;
        
        temp0 = sparse((1:size0)', (1:size0)', delta_t * exp(firm0*param), size0, size0);
        temp2 = delta_t*exp(firm2*param);
        temp2_g = sparse((1:size2)', (1:size2)', exp(-temp2).*temp2./(1-exp(-temp2)), size2, size2);
        temp2_h = sparse((1:size2)', (1:size2)', (1-temp2-exp(-temp2))./(1-exp(-temp2)), size2, size2);
        g0 = -temp0*firm0;
        g2 = temp2_g*firm2;
        h0 = firm0'*g0;
        g0 = sum(g0)';
        h2 = firm2'*temp2_h*g2;
        g2 = sum(g2)';
        h = -(h0+h2);
        g = -(g0+g2);
        
        A = h + mu*eye(size(h,1));
    end
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
    
    if intensitytype==1
        f1 = -sum(exp(firm0*param1)*delta_t)+sum(log(1-exp(-exp(firm1*param1)*delta_t)))-sum(exp(firm2*param1)*delta_t);
    else
        f1 = -sum(exp(firm0*param1)*delta_t)+sum(log(1-exp(-exp(firm2*param1)*delta_t)));
    end
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