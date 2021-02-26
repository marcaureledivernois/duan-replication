function L = log_likelihood_fun(x, F, LiabilitiesOther, r, market_value_equity, AssetsTotal)

% mu = x(1) ;
% sigma = x(2);
% delta = x(3) ;

n = size(market_value_equity,1);
V = zeros(n,1);
D = F +  x(3)*LiabilitiesOther;

xL = 0 * ones(1,1);
xU =  9999999 * ones(1,1);
warning off
V =solve_dichotomy(@afunction,xL,xU,0.0001,D,r,market_value_equity,x(2));

W(:,1) = log(V(2:end,1)./V(1:end-1,1).*AssetsTotal(1:end-1,1)./AssetsTotal(2:end,1)) ...
    - 3/12*(x(1) - 0.5*x(2)^2);

d_t(:,1) = (x(2))^(-1)*(log(V(2:end,1)./D(2:end,1)) + (r(2:end) + 0.5*x(2)^2));

L = -(n-1)/2*log(2*pi) -0.5*(n-1)*(log(3/12*x(2)^2)) -sum(W.^2./(2*3/12*x(2)^2),1) ...
    -sum(log(V(2:end,1)./AssetsTotal(2:end,1)),1) -sum(log(normcdf(d_t(2:end))),1);

L = -L;

end