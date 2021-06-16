function intenField = intenF(this,field,n)

intenField = real( (2*this.const_eps0*n*this.const_c) .* conj(field).*field );

end