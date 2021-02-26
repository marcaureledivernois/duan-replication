function xsol=solve_dichotomy(fun,xL,xU,tol,varargin)
% solve_dichotomy: computes the zeros of a function using the dichotomy algorithm. Returns an
% error if one of the bounds for x is not verified. Expects all elements to
% be column vectors
fL = feval(fun,xL,varargin{:}); % evaluate function at lower boundary
fU = feval(fun,xU,varargin{:}); % and at upper boundary

% [fL fU]

% check if problem makes sens. Sides of function must be of opposite signs
x=fL.*fU;
e=x>0;
idx=find(e==1);
se=sum(e);
if se>0
    [idx fL(idx) fU(idx)]
    error('For the indices indicated, lower or upper bounds do not delimit boundary, problem in solve_dichotomy')
end
d=xU-xL;
dmax=max(d);
NIter=floor((log(dmax/tol))/log(2))+1; % this is the required number of iterations to reach required precision
warning(['Dichotomy: Required number of iterations to reach desired accuracy ' num2str(NIter)]);

for j=1:NIter
    xmi = (xU+xL)/2; % get middle of upper and lower bound
    fmi = feval(fun,xmi,varargin{:}); % evaluate function at midpoint
    e   = fL.*fmi>0; % where lower bound needs to be moved
    xL  = e.*xmi +(1-e).*xL; % move lower and upper bounds where required
    xU  = (1-e).*xmi + e.*xU;
    fL  = e.*fmi +(1-e).*fL; % update values of funcsiotns where required
    fU  = (1-e).*fmi + e.*fU; 
end
xsol=(xU+xL)/2; % since ingredients already computed, cost of geting one last improvement is cheap
