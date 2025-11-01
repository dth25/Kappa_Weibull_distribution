function y  = expk(x, k)
%   EXPK is the kappa exponential: y = expk(x,k)

if k~=0
y = ( sqrt(1 + (k*x).^2) + k*x ).^ (1/k);
else
    y=exp(x);
end;

