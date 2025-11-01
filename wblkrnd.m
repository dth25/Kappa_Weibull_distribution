function r = wblkrnd(xs, m, k, varargin)
%WBLKRND Random arrays from the kappa-Weibull distribution.
%   R = WBLKRND(XS,M,K) returns a random number chosen from the
%   kappa-Weibull distribution with scale parameter XS and shape parameter M.
%
%   R = WBLKRND(XS,M,K,A,B,...) or R = WBLKRND(XS,M,K,[A,B,...]) returns an
%   A-by-B-by-... array.
%
%   See also WBLKCDF, WBLRND.

%   WBLKRND is based on WBLRND and uses the inversion method.

%   References:
%     [1] Clementi F., Gallegati M. and Kaniadakis G. A kappa-generalized
%         statistical mechanics approach to income analysis, J. Stat. Mech.
%         (2009) P02037

%   $Revision: 1.1 $  $Date: 2013/03/13 $

% Check input arguments
if nargin < 3
    error('stats:wblkrnd:TooFewInputs','Requires at least three input arguments.');
end

% k>1 is valid. In [1] they restrict k to [0,1) 
if k<0
    error('stats:wblkrnd:ArgumentsRange','Kappa should be positive.');
end

if xs<0 || m<0
    error('stats:wblkrnd:ArgumentsRange','Arguments should be positive.');
end

% Generate uniform random values, and apply the kappa-Weibull inverse CDF.
r =  xs * (-lnk(1-rand(varargin{:}), k)).^(1/m);
