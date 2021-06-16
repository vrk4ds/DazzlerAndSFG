function propagate(this,textFlag)

if ~exist('textFlag','var')
    textFlag = 1;
end

nFields = size(this.eField,3);

for iz = 2:(this.grid_nZSteps+1)
    
    for iF = 1:nFields
        this.eField(1,iz,iF,:) = this.ifft( this.eField(2,iz-1,iF,:) .*...
            permute(exp(1i.*this.list_k(iF,:).*dzStep(this,iz)),[4,3,1,2])...
            );
    end
    if iz <= this.grid_nZSteps + 1
        this.RKstep(iz);
        for iF = 1:nFields
            this.eField(2,iz,iF,:) = this.fft( this.eField(1,iz,iF,:) );
        end
    end
    
    if textFlag
        disp(['Z-Step: ',num2str(iz)]);
    end
    
end


end

function dz = dzStep(this,iz)

if iz == 2 || iz == this.grid_nZSteps + 1
    dz = this.grid_dz/2;
else
    dz = this.grid_dz;
end

end