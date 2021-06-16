function intenPeak = intenPK(ePulse,mRad,tau)

intenPeak = ( 2*ePulse*sqrt(log(16)) )/( (pi^(3/2))*(mRad^2)*tau );

end