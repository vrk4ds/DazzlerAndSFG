function genGrids(this,nPts,dt,nZ)

if exist('nPts','var') && ~isempty(nPts)
    this.grid_nPts = nPts;
else
    if isempty(this.grid_nPts)
        this.grid_nPts = 2^14;
    end
end

if exist('dt','var') && ~isempty(dt)
    this.grid_dt = dt;
else
    if isempty(this.grid_dt)
        this.grid_dt = this.props_taus(1)/10;
    end
end

if exist('nZ','var') && ~isempty(nZ)
    this.grid_nZSteps = nZ;
else
    if isempty(this.grid_nZSteps)
        this.grid_nZSteps = 100;
    end
end

nFields = length(this.props_lams);

this.grid_dz = this.props_len / (this.grid_nZSteps - 1);
this.grid_dw = (2*pi)/(this.grid_nPts * this.grid_dt);


this.list_t = this.grid_dt*((-((this.grid_nPts/2)-1)):(this.grid_nPts/2));
this.list_dOmega = this.grid_dw*((-((this.grid_nPts/2)-1)):(this.grid_nPts/2));

if strcmpi(this.props_mixType,'SFG')
    
    for ii = 1:nFields
        
        this.list_omega(ii,:) = this.list_dOmega + this.props_omegas(ii);
        this.list_lambda(ii,:) = (2*pi*this.const_c)./this.list_omega(ii,:);
        
        if ii ~= nFields
            this.list_k(ii,:) = ...
                ( (2*pi)./this.list_lambda(ii,:) )...
                .* this.eqns_Index{ii}(this.list_lambda(ii,:))...
                - (this.list_dOmega .* this.eqns_k{1,2}(this.props_omegas(1)));
        elseif ii == nFields
            this.list_k(ii,:) = ...
                ( (2*pi)./this.list_lambda(ii,:) )...
                .* this.eqns_Index{ii}(this.list_lambda(ii,:),this.props_theta)...
                - (this.list_dOmega .* this.eqns_k{1,2}(this.props_omegas(1)));
        end
        
    end
    
elseif strcmpi(this.props_mixType,'SHG')
    
    for ii = 1:nFields
        
        this.list_omega(ii,:) = this.list_dOmega + this.props_omegas(ii);
        this.list_lambda(ii,:) = (2*pi*this.const_c)./this.list_omega(ii,:);
        
        if ii ~= nFields
            this.list_k(ii,:) = ...
                ( (2*pi)./this.list_lambda(ii,:) )...
                .* this.eqns_Index{ii}(this.list_lambda(ii,:))...
                - (this.list_dOmega .* this.eqns_k{1,2}(this.props_omegas(1)));
        elseif ii == nFields
            this.list_k(ii,:) = ...
                ( (2*pi)./this.list_lambda(ii,:) )...
                .* this.eqns_Index{ii}(this.list_lambda(ii,:),this.props_theta)...
                - (this.list_dOmega .* this.eqns_k{1,2}(this.props_omegas(1)));
        end
        
    end
    
end

end