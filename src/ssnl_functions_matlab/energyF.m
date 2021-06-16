function energyField = energyF(this,field,n,mRad)

energyField = ((pi/2)*mRad^2) * sum(this.intenF(field,n)) * this.grid_dt;

end