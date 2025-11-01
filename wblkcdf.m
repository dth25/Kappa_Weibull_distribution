function y = wblkcdf(x, xs, m, k)
%   WBLKCDF is the CDF of the kappa Weibull probability distribution
%   
%   INPUT
%
%   X:     X values
%   XS:     Weibull scale parameter (scalar)
%   M:      Weibull Modulus (scalar)
%   K:      Kappa values (scalar)
%
%   OUTPUT
%
%   Y:      Y=P(X <=x; xs,m,k)

xn = x/xs; 
xnm = xn.^m; 
y = 1 -  expk(-xnm,k) ;

