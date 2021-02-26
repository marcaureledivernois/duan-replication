function x = TestDichotomy()
% tests if program that computes zero of a vector function does so
% correctly
clc

xL = -7 * ones(2,1);
xU =  7 * ones(2,1);
s  = (1:1:500)';
s=s/100;
tol=0.0001;
A=1;
warning on
D = [50;60];
r = 0.05
E = [50;30];
sigma = 0.05;
tic
x=solve_dichotomy(@afunction,xL,xU,tol,D,r,E,sigma);
toc
end
%=================================================================


