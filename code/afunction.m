function y=afunction(x,D,r,E,sigma)

d1 = (log(x./D)+(r+sigma^2))./sigma; 
d2 = d1 - sigma;
y  = x.*normcdf(d1)-exp(-r).*D.*normcdf(d2)-E;
end