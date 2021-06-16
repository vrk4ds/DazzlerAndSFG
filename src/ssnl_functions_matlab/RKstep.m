function RKstep(this,zStep)
%%% Speecialized RK4 equation to deal with class structure

nVars = size(this.eField,3);
N = this.grid_nPts;

if nVars == 2
    ordF = [1,2;1,1];
elseif nVars == 3
    ordF = [2,3;1,3;1,2];
end

rk0 = zeros(1,1,nVars,N);
rk1 = zeros(1,1,nVars,N);
rk2 = zeros(1,1,nVars,N);
rk3 = zeros(1,1,nVars,N);

for ii = 1:nVars
    rk0(1,1,ii,:) = this.grid_dz * this.eqns_NL{ii}(...
        this.eField(1,zStep,ordF(ii,1),:),...
        this.eField(1,zStep,ordF(ii,2),:) );
end

for ii = 1:nVars
    tmpField = this.eField(1,zStep,:,:) + rk0/2;
    rk1(1,1,ii,:) = this.grid_dz * this.eqns_NL{ii}(...
        tmpField(1,1,ordF(ii,1),:),...
        tmpField(1,1,ordF(ii,2),:) );
end

for ii = 1:nVars
    tmpField = this.eField(1,zStep,:,:) + rk1/2;
    rk2(1,1,ii,:) = this.grid_dz * this.eqns_NL{ii}(...
        tmpField(1,1,ordF(ii,1),:),...
        tmpField(1,1,ordF(ii,2),:) );
end

for ii = 1:nVars
    tmpField = this.eField(1,zStep,:,:) + rk2;
    rk3(1,1,ii,:) = this.grid_dz * this.eqns_NL{ii}(...
        tmpField(1,1,ordF(ii,1),:),...
        tmpField(1,1,ordF(ii,2),:) );
end

this.eField(1,zStep,:,:) = this.eField(1,zStep,:,:) + rk0/6 + rk1/3 + rk2/3 + rk3/6;

end