function aFieldPeak = fieldPK(this,ePulse,mRad,tau,n)

aFieldPeak = sqrt( this.intenPK(ePulse,mRad,tau) ./ (2*this.const_eps0*n*this.const_c) );

end