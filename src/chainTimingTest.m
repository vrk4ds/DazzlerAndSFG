%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Split Step Nonlinear Propagation Working file %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Constants %%%%

clear
% Length
m = 10^0; mm = 10^-3*m; um = 10^-6*m; nm = 10^-9*m;
% Time
s = 10^0; ps = 10^-12*s; fs = 10^-15*s;
% Energy
J = (m^2)/(s^2); mJ = 10^-3 * J; uJ = 10^-6 * J;
% Physics Constants
c0 = 299792458*(m/s); eps0 = 8.854187817*10^-12/m;

%% Chain Timing Testing

if ~exist('this','var')
    load('chain.mat');
end


timing = zeros(4,100);

for ii = 1:100
    a = tic;

    b = tic;
    ssnlObj = ssnl(this,'copy');
    timing(2,ii) = toc(b);

    % figure(1); prev.plotField({1,101},3,[-20 20],0)
    % title('Crys 1, Field 3 (515 nm), Time')
    % figure(2); prev.plotField({2,101},3,[514.5 515.5],0)
    % title('Crys 1, Field 3 (515 nm), Spectrum')

    c = tic;
    ssnlObj = ssnl('chain',ssnlObj,rulesSHG);
    timing(3,ii) = toc(c);

    d = tic;
    ssnlObj.propagate(0)
    timing(4,ii) = toc(d);

    % figure(3); prev.plotField({1,101},2,[-20 20],0)
    % title('Crys 2, Field 2 (257.5 nm), Time')
    % figure(4); prev.plotField({2,101},2,[257 258],0)
    % title('Crys 2, Field 2 (257.5 nm), Spectrum')

    timing(1,ii) = toc(a);
    ii
end

mean(timing,2)
sum(timing,2)/60
