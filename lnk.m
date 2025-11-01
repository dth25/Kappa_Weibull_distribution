function [ y ] = lnk( x,k )
%   LNK is the kappa logarithm Y= ( X^K - X^(-K) )/2K
%
%   The classical loGarithm follows at the limit K -->0
%
%   INPUT
%
%   X:      Vector of real values
%   K:      Kappa number in [0, 1)
%
%   OUTPUT
%
%   Y:      Kappa logarithm of X
if k~=0
    y = (x.^k - x.^(-k))./( 2*k);
elseif k==0
    y=log(x);
end

end

