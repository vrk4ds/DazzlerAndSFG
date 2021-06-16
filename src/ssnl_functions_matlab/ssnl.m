%% Beam Creation and Propagation Master Class


% This class holds all the creation and propagation functions and is
% designed to be fully extensible and give access to all neccessary
% functions and properties at any given time.



classdef ssnl < handle
    
    %% Properties %%
    properties
        
        props_crys
        props_len
        props_dNL
        props_theta
        props_mixType
        props_lams
        props_omegas
        props_ks
        props_taus
        props_energies
        props_spotRad
        props_specPhases
        
        grid_nPts
        grid_dt
        grid_dz
        grid_nZSteps
        grid_dw
           
        list_t
        list_lambda
        list_omega
        list_dOmega
        list_k
        
        eqns_Index
        eqns_k
        eqns_NL
        eqns_phase
        
        chain_Fields
                        
        eField % Numbering: Fourier Space(t=1), z, field, values (ie. [t,dz,lam1,:])
        
    end
    
    properties (Constant)
        
        const_c = 299792458;
        const_eps0 = 8.854187817 * 10^-12;
        const_fwhm = 4 * log(2);
        
    end
    
    
    %% Methods %%
    methods
        
        %%%% Initialization Function %%%%
        function this = ssnl(varargin)
            
            if nargin == 0
                return
            elseif nargin == 1
                if isstruct(varargin{1}) || isa(varargin{1},'ssnl')
                    
                    this.inputProperties(varargin{1});
                                        
                elseif strcmpi(varargin{1},'gui')
                    
                    this.setupGUI();
                    
                end
                                
            elseif nargin == 2
                if isstruct(varargin{1}) || isa(varargin{1},'ssnl')
                    
                    this.inputProperties(varargin{1});
                    
                    if strcmpi(varargin{2},'copy') || varargin{2} == 1
                        return
                    end
                    
                end
                
            elseif nargin == 3
                
                if strcmpi(varargin{1},'chain') || strcmpi(varargin{1},'chainprev')
                    
                    this.chainPrev(varargin{2},varargin{3});
                    
                end
                    
                
            end
            
            if isempty(this.eqns_Index)
                this.genEqns;
            end
            
            this.genGrids;
            this.genField;
            
        end
        
        
        setupGUI(this);
        chainPrev(this,prev,rules);
        
        plotField(this,fInd,fNums,range,norm);
        propagate(this,textFlag);
        
        genEqns(this);
        genGrids(this,varargin)
        genField(this);
                        
        aFieldPeak = fieldPK(this,ePulse,mRad,tau,n);
        energyField = energyF(this,field,n,mRad);
        intenField = intenF(this,field,n);
        
        inputProperties(this,vec)
        output = outputProperties(this,str)
        
        
    end
    
    methods (Static)
       
        varargout = indexEqns(crysName);
        
        fwhm = FWHM(lst);
        
        field = fft(field);
        field = ifft(field);
        
        intenPeak = intenPK(ePulse,mRad,tau);
        
        outStruct = chainStruct();
        
    end
    
end