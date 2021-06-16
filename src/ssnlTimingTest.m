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

%% Timing test for SSNL class progamatically %%%%

N = 1;
timing = zeros(4,N);

for ii = 1:N
    
    a = tic;
    
    b = tic;
    ssnlInit.props_crys = 'BBO';
    ssnlInit.props_len = 2*mm;
    ssnlInit.props_theta = 23.29;
    ssnlInit.props_mixType = 'SFG';
    ssnlInit.props_lams = [1030*nm,1030*nm,515*nm];
    ssnlInit.props_ks = (2*pi)./ssnlInit.props_lams;
    ssnlInit.props_omegas = c0 .* ssnlInit.props_ks;
    ssnlInit.props_taus = [330*fs,330*fs,20*fs];
    ssnlInit.props_energies = [25*uJ,25*uJ,0];
    ssnlInit.props_spotRad = 400*um;
    
%     tayScale = 9.0 + ((rand-1)*0.25);
    tayScale = 1;
    
    ssnlInit.props_specPhases =...
        [-3.27*(ps^2),(0.42*tayScale)*(ps^3),0,0;...
        3.27*(ps^2),(-0.42*tayScale)*(ps^3),0,0;...
        0,0,0,0];
    
    npts = 2^14;
    
    
    ssnlObj = ssnl(ssnlInit,'copy');
    timing(2,ii) = toc(b);
    
    c = tic;
    ssnlObj.genEqns;
    ssnlObj.genGrids(npts);
    ssnlObj.genField;
    timing(3,ii) = toc(c);
    
    d = tic;
    ssnlObj.propagate(0);
    timing(4,ii) = toc(d);
    
    timing(1,ii) = toc(a);
    ii
end

mean(timing,2)
sum(timing,2)/60;

ssnlObj.energyF(ssnlObj.eField(1,101,3,:),ssnlObj.eqns_Index{3}(ssnlObj.props_lams(3),ssnlObj.props_theta),ssnlObj.props_spotRad)/uJ

%%%% Plotting %%%%

figure(1);clf;
ssnlObj.plotField({'t',1:101},3,[-25 25],0)

figure(2);clf;
ssnlObj.plotField({'w',1:101},3,[514 516],0)







