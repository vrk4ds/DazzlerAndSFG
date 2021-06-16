function genEqns(this)

syms w theta l lCtr field1 field2 field3 dOmega kk2 kk3 kk4 kk5 phiTay

%% index of refraction equations %%

[nO,nE_Theta,this.props_dNL] = this.indexEqns(this.props_crys);

if strcmpi(this.props_mixType,'SFG')
    
    this.eqns_Index{1} = matlabFunction(nO); % f( lambda )
    this.eqns_Index{2} = matlabFunction(nO); % f( lambda )
    this.eqns_Index{3} = matlabFunction(nE_Theta); % f( lambda, theta )
    
elseif strcmpi(this.props_mixType,'SHG')
    
    this.eqns_Index{1} = matlabFunction(nO); % f( lambda )
    this.eqns_Index{2} = matlabFunction(nE_Theta); % f( lambda, theta )
    
end

%% k-vector equations and derivatives %%

if strcmpi(this.props_mixType,'SFG')
    
    this.eqns_k{1,1} =...
        matlabFunction((w/this.const_c)*nO(2*pi*this.const_c/w)); % f( omega )
    this.eqns_k{1,2} =...
        matlabFunction(diff(this.eqns_k{1,1},w)); % f( omega )
    this.eqns_k{1,3} =...
        matlabFunction(diff(this.eqns_k{1,2},w)); % f( omega )
    
    this.eqns_k{2,1} =...
        matlabFunction((w/this.const_c)*nO(2*pi*this.const_c/w)); % f( omega )
    this.eqns_k{2,2} =...
        matlabFunction(diff(this.eqns_k{1,1},w)); % f( omega )
    this.eqns_k{2,3} =...
        matlabFunction(diff(this.eqns_k{1,2},w)); % f( omega )
    
    this.eqns_k{3,1} =...
        matlabFunction((w/this.const_c)*nE_Theta(2*pi*this.const_c/w,theta),...
        'Vars',{w,theta}); % f( omega, theta )
    this.eqns_k{3,2} =...
        matlabFunction(diff(this.eqns_k{2,1},w),'Vars',{w,theta}); % f( omega, theta )
    this.eqns_k{3,3} =...
        matlabFunction(diff(this.eqns_k{2,2},w),'Vars',{w,theta}); % f( omega, theta )
    
elseif strcmpi(this.props_mixType,'SHG')
    
    this.eqns_k{1,1} =...
        matlabFunction((w/this.const_c)*nO(2*pi*this.const_c/w)); % f( omega )
    this.eqns_k{1,2} =...
        matlabFunction(diff(this.eqns_k{1,1},w)); % f( omega )
    this.eqns_k{1,3} =...
        matlabFunction(diff(this.eqns_k{1,2},w)); % f( omega )
    
    this.eqns_k{2,1} =...
        matlabFunction((w/this.const_c)*nE_Theta(2*pi*this.const_c/w,theta),...
        'Vars',{w,theta}); % f( omega, theta )
    this.eqns_k{2,2} =...
        matlabFunction(diff(this.eqns_k{2,1},w),'Vars',{w,theta}); % f( omega, theta )
    this.eqns_k{2,3} =...
        matlabFunction(diff(this.eqns_k{2,2},w),'Vars',{w,theta}); % f( omega, theta )
    
end

%% non-linear mixing equations

if strcmpi(this.props_mixType,'SFG')
    
    % This is the nonlinear coefficient based on the crystal axis and lambda
    sNL = this.props_dNL * 1i * [...
        (2*this.props_ks(1)/this.eqns_Index{1}(this.props_lams(1))),...
        (2*this.props_ks(2)/this.eqns_Index{2}(this.props_lams(2))),...
        (2*this.props_ks(3)/this.eqns_Index{3}(this.props_lams(3),this.props_theta))...
        ];
    
    % Each of these equations is a function of two fields. I have labeled them
    % here for clarity on which ones are used in each equation. The choice of
    % which equation to use is how the field you are solving for is decided.
    this.eqns_NL{1} = matlabFunction( sNL(1) .* conj(field2) .* field3 );
    this.eqns_NL{2} = matlabFunction( sNL(2) .* conj(field1) .* field3 );
    this.eqns_NL{3} = matlabFunction( sNL(3) .* field1 .* field2 );
    
elseif strcmpi(this.props_mixType,'SHG')
    
    % This is the nonlinear coefficient based on the crystal axis and lambda
    sNL = this.props_dNL * 1i * [...
        (2*this.props_ks(1)/this.eqns_Index{1}(this.props_lams(1))),...
        (this.props_ks(2)/this.eqns_Index{2}(this.props_lams(2),this.props_theta))...
        ];
    
    % Each of these equations is a function of two fields. I have labeled them
    % here for clarity on which ones are used in each equation. The choice of
    % which equation to use is how the field you are solving for is decided.
    this.eqns_NL{1} = matlabFunction( sNL(1) .* conj(field2) .* field3 );
    this.eqns_NL{2} = matlabFunction( sNL(2) .* field1 .* field2 );
    
end

%% Taylor expansion of the phase of the pulse

dOmega(l,lCtr) = (2*pi*this.const_c) * ( (1/lCtr) - (1/l) );
this.eqns_phase = matlabFunction(...
    ((kk2/factorial(2)) .* (dOmega(l,lCtr)^2)) +...
    ((kk3/factorial(3)) .* (dOmega(l,lCtr)^3)) +...
    ((kk4/factorial(4)) .* (dOmega(l,lCtr)^4)) +...
    ((kk5/factorial(5)) .* (dOmega(l,lCtr)^5)),...
    'Vars',{l,lCtr,kk2,kk3,kk4,kk5});


end