function [y, g]=lnL_wbk(x,data,cens, freq)
%   LNL_WBK calculates the negative log-likelihood of the kappa Weibull
%   distribution for a given data set
%
%   INPUT
%
%   X:      Kappa Weibull parameters (scale, Weibull modulus, kappa)
%   DATA:   Nx1 vector of measurements
%
%   OUTPUT
%
%   Y:      Negative log-likelihood scaled by the number of points
%   G:      Gradient of NLL with respect to the parameters

%   Use a negative log-likelihood per sample point

N=length(data);

t0=x(1); m=x(2); kk=x(3);


%   If the parameter k is very small, use the Weibull log-likelihood

%   Correct zero data before taking the logarithm

i0 = data==0; 
datac = data; 
datac(i0) = eps; 

lndatc = log(datac); 
lnt0 = log(t0); 

if kk > 10^-6
    tt =  (data/t0).^(m);
    dsq =  1 + (kk*tt).^2 ;
    
    yv1= - N *log(m) + N* m * lnt0 ;
    
    yv2= - (m-1)* lndatc;
    
    yv3= 0.5* log( dsq );
    
    temp = expk (-tt,kk);
    iz = find(temp ==0);
    %   For values of tt >>1  expk(-tt) may give zero; to avoid this note that
    %   log(expk(-tt)) = log(exp(-tt)) = -tt, for tt >> 1.
    yv4= -  log( temp );
    if ~isempty(iz)
        yv4(iz) = - tt(iz);
    end
    yv= yv1 + sum(yv2 + yv3  + yv4);
    y=yv/N;
    
    %   The term g_ij is the derivative of yvj with respect to parameter i
    %   where p_1 = t0, p_2 = m, p_3 = kappa.
    
    sdsq = sqrt(dsq);
    
    g11 = N*m/t0;
    g13 = - (m * kk^2/t0) *(tt.^2) ./ dsq ;
    g14 = -(m /t0) * tt ./ sdsq  ;
    g1 = g11 + sum(g13+g14);
    
    g21 = -N/m + N * lnt0;
    g22 = - lndatc ;
    g23 =  ((kk*tt).^2 ) .*  (lndatc - lnt0)  ./ dsq ;
    g24 =  (lndatc - lnt0) .* tt ./sdsq ;
    g2 = g21 + sum(g22+g23+g24);
    
    
    g32= kk * (tt.^2) ./ dsq ;
    
    %    f3 = - ((data/t0).^(m)) ./ sqrt(1 + (kk^2)*(data/t0).^(2*m)) ;
    g33_num = -( -kk^2 *(tt.^2) + kk*tt .* sdsq + log(sdsq -kk*tt) + ...
        kk^2 * log(sdsq -kk*tt) .* (tt.^2) - ...
        kk* log(sdsq -kk*tt).*sdsq .* tt );
    g33_den = kk^2 * sdsq .* (kk*tt - sdsq);
    
%     if any(~isfinite(g33_num)) || any(~isreal(g33_num))
%         disp('g33_num');
%     elseif any(~isfinite(g33_den)) || any(~isreal(g33_den))
%         disp('g33_den');
%     end
    
    
    g33 = g33_num./g33_den;
    
    g3 = sum(g32+g33);
    
    g = [ g1 g2 g3]/N;
 
    %   For large t try a different rescaling of the terms that does not
    %   involve infinities
    if any(~isfinite(g)) || any(~isreal(g))
        disp('alt'); 
        ttidsq =  1./( kk^2 + (tt).^(-2) ) ;
        
%        min(ttidsq), max(ttidsq),
        sttidsq = sqrt(ttidsq);
        
        g13 = - (m * kk^2/t0) * ttidsq;
        g14 = -(m /t0) * sttidsq ;
        g1 = g11 + sum(g13+g14);
        
        g23 =  (kk.^2 ) .* (lndatc - lnt0)  .* ttidsq ;
        g24 =  (lndatc - lnt0) .* sttidsq ;
        g2 = g21 + sum(g22+g23+g24);
        
        g32= kk * ttidsq ;
        
        
        %    f3 = - ((data/t0).^(m)) ./ sqrt(1 + (kk^2)*(data/t0).^(2*m)) ;
        g33_num = -( -kk^2 * ttidsq + kk* sttidsq + log(sdsq -kk*tt) .* ...
           (ttidsq./(tt.^2) + kk^2 .* ttidsq - kk .* sttidsq ) );
        
        g33_den = kk^2 *  (kk*sttidsq - 1);
        
        g33 = g33_num./g33_den;
        
%         if any(~isfinite(g33_num)) || any(~isreal(g33_num))
%             disp('g33_num second');
%         elseif any(~isfinite(g33_den)) || any(~isreal(g33_den))
%             disp('g33_den second');
%         end

%   The logarithm in g33_num can become -Inf for large values of tt. This
%   divergence, however, is canceled by the faster divergence of the term
%   that multiplies the logarithm, which for large tt behaves as tt^-2. 

         if any(~isfinite(g33_num)) || any(~isreal(g33_num))
             g33=0;
         end
       
        g3 = sum(g32+g33);
        
        g = [ g1 g2 g3]/N;
            
    end
else
    
    y = wbllike([t0 m], data)/N; % Negative Log_likelihood for Weibull per sample point
    
end

if y==Inf
    y=realmax;
elseif y== -Inf
    y=realmin;
end




