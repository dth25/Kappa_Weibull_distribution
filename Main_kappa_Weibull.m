% The function MAIN_KAPPA_WEIBULL calculates the fit of the data set DATA to 
% the WEIBULL and the KAPPA Weibull distributions.  The empirical CDF of the 
% data is first constructed. 
%  
% The kappa-Weibull CDF fit is performed using the MLE method by minimization
% of the Negative Log-Likelihood (NLL) of the KAPPA Weibull distribution.
% The optimization uses a two-step approach: 
% First, it employs the function FMINCON using the ACTIVE-SET method 
% and is then followed by FMINSEARCH (Simplex search method), which uses 
% the output of FMINCON as initial condition.  Gradient information is 
% supplied to FMINCON via the objective function lnL_wbk 
%  
% Plots of the CDF for the empirical distribution as well as the best-fit
% Weibull and kappa-Weibull distributions are generated. The functions   
% Phi(tau) (Weibull plots) for the three CDFs are also plotted. Finally,
% the quantile-quantile plot (empirical to kappa-Weibull) is generated.  
%
%   The code has been tested with Matlab R2015b
% 
%==========================================================================
%   INPUT VARIABLES
%==========================================================================
%   DATA:           Array of data values
%   DATANAM:        String array containing name of data (used for
%   labeling)
%   
%==========================================================================
%   OUTPUT VARIABLES
%==========================================================================
%
%   BETA_W:         Parameters of Weibull
%   BETA_KW:        Parameters of Kappa-Weibull
%   NLL_KW:         Negative Log-Likelihood of Fitted Kappa-Weibull
%   NLL_W:          Negative Log-Likelihood of  Fitted Weibull
%==========================================================================
%   EXAMPLE:
%==========================================================================
%   load Cairo_wind_speed; 
%   [beta_w, beta_kw, NLL_W, NLL_KW] = Main_kappa_Weibull(wind, 'Cairo wind');
%   close all; clearvars;
%   load Carbon_fiber
%   [beta_w, beta_kw, NLL_W, NLL_KW] = Main_kappa_Weibull(data_carbon, 'Tensile strength')
%
%==========================================================================
%   FUNCTIONS USED: expk, lnk, lnL_wbk, wblkcdf 
%==========================================================================
%
%	Copyright (C) 2022 Dionisis Hristopulos
%   email: dchristopoulos@tuc.gr
% 
%	Last Modified: July 2022
% 
%==========================================================================
%	REFERENCES (If you use this code please cite the following)
%==========================================================================
%
% [1] D.T. Hristopulos and A. Baxevani, "Kaniadakis functions beyond 
% statistical mechanics: weakest-link scaling, power-law tails, and 
% modified lognormal distribution," Entropy, 2022. 
%  
% [2] D. T. Hristopulos, M. P. Petrakis, and G. Kaniadakis, "Finite-size 
% effects on return interval distributions for weakest-link-scaling 
% systems," Physical Review E 89, 052142, 28 May 2014.
% https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.052142
% 
% [3] D. T. Hristopulos, M. P. Petrakis, and G. Kaniadakis, 
% "Weakest-link scaling and extreme events in finite-sized systems," 
% Entropy, 17(3):1103-1122, 2015. https://doi.org/10.3390/e17031103
% 
% For more information on the datasets used to test this code see 
% references in [1].

function [beta_w, beta_kw, NLL_W, NLL_KW] = Main_kappa_Weibull(data, datanam)

data(data==0)=[];

L  = length(data);

%==========================================================================
%   Calculate empirical probability distribution
%==========================================================================

[Fwt,xwt,flo,fup]=ecdf(data);

if sum(data <0) > 0
    disp('Negative data - I will stop here'); 
    return;
end

%==========================================================================
%   Kappa-Weibull distribution

optionsk = optimoptions('fmincon', 'Algorithm','active-set','GradObj','on',...
    'MaxFunEvals',20000,'MaxIter',20000,'TolFun',1.E-5);
%   The optimization is performed using Gradient information 

ts1=prctile(data,30); ts2=prctile(data,70);
mmin=0.3; mmax=5;
kmin= 0.00001; kmax=2;
x0 = [ max(eps, mean(data));  0.7;  0.5];

error = 100; ic = 1; 

LB = [ max(eps, ts1); mmin; kmin];
UB = [ max(eps, ts2); mmax; kmax];


while error> 0.0001

[x,fval,~,output,lambda,grad] = fmincon(@(x) lnL_wbk(x,data), ...
    x0,[],[],[],[], LB, UB, [], optionsk);

x01 = x;
[x,fval,exitflag,output] = fminsearch(@(x) lnL_wbk(x,data),x01);
x(x<0)=0;   % This takes care of very small negative values for kappa

error = sqrt(norm (( x - x0)./x0) ); 

x0 = x; ic  = ic+1;

end
ic = ic -1; 
beta_kw = x;
NLL_KW = L * fval; % Multiply with sample size because lnL_wbk maximizes likelihood per sampling point

cdf_wk=wblkcdf(xwt,beta_kw(1),beta_kw(2),beta_kw(3));

%==========================================================================
%   Fit data to Weibull distribution
%==========================================================================

%   Correct zero data before taking the logarithm

i0 = data==0; 
datac = data; 
datac(i0) = eps; 
[beta_w] = mle(datac,'Distribution','wbl');

cdf_w = wblcdf(xwt,beta_w(1),beta_w(2));

NLL_W = wbllike(beta_w,datac);
%==========================================================================
%   Plots of cumulative probability distribution
%==========================================================================
figure(1);
np=length(cdf_wk);

semilogx(xwt, Fwt, 'b-','LineWidth', 1.5); hold on;
semilogx(xwt, cdf_w, 'g:','LineWidth', 1.5); hold on;
semilogx(xwt, cdf_wk, 'r--','LineWidth', 1.5); hold on;
%semilogx(xwt,cdf_wk2,'m.','LineWidth',2.5); hold on;
axis tight
legend('Empirical','Weibull','\kappa Weibull','Location','Best');
grid on;
xlabel(datanam);
ylabel('CDF');
% title(['M_{c}=0',blanks(1),'#Events=',num2str(np)]);

title([' Scale=',num2str(beta_kw(1),'% 10.2E'),',', ...
    ' m=',num2str(beta_kw(2),'% 10.2f'),',', ...
    ' \kappa=',num2str(beta_kw(3),'% 10.2f')]);
set(gca, 'FontSize', 14);

beta_kw(2:3),

%==========================================================================
%   Weibull Plots
%==========================================================================

xx=log(xwt);
yy_e= log(log(1./(1-Fwt)));
yy_w= log(log(1./(1-cdf_w)));
yy_wk = log(log(1./(1-cdf_wk)));
%yy_wk2 = log(log(1./(1-cdf_wk2)));

%==========================================================================
figure(2);
plot(xx, yy_e, 'b-', 'LineWidth', 1.5); hold on;
plot(xx, yy_wk, 'r--', 'LineWidth', 1.5); hold on;
plot(xx, yy_w, 'g-.', 'LineWidth', 1.5);

xlabel(['ln(' datanam ')']);
ylabel(['\Phi(' datanam ')']);
legend('Empirical', '\kappa Weibull', 'Weibull', 'Location', 'Best');
grid on;
title('Weibull plot')
% title([' Scale=',num2str(beta_kw(1),'% 10.2E'),',', ...
%     ' m=',num2str(beta_kw(2),'% 10.2f'),',', ...
%     ' \kappa=',num2str(beta_kw(3),'% 10.2f')]);
set(gca, 'FontSize', 14);
axis tight

%==========================================================================
%   CONSTRUCTION OF Q-Q PLOTS
%==========================================================================

pp=1:1:99;

P_data=prctile(data,pp);
xs=beta_kw(1); m=beta_kw(2); k=beta_kw(3);

pval=1-pp/100;
P_model=xs*((lnk(1./pval,k)).^(1/m));

%==========================================================================
figure(3);  % Q-Q plot 
loglog(P_data, P_model, 'ro', 'MarkerFace', 'r', 'MarkerSize', 6); 
hold on;
xx0=min(min(P_data), min(P_model));
yy0=max(max(P_data), max(P_model));
line([xx0 yy0], [xx0 yy0],'LineWidth',2);
set(gca,'FontSize',14);

xlabel(datanam,'FontSize',14);
ylabel('\kappa-Weibull','FontSize',14);
title('Q-Q plot');
% title([' \tau_0=', num2str(beta_kw(1), '% 10.2E'),',', ...
%     ' m=', num2str(beta_kw(2), '% 10.2f'),',', ...
%     ' \kappa=', num2str(beta_kw(3), '% 10.2f')]);
axis tight;
set(gca, 'FontSize', 14);
axis equal
grid on
