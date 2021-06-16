import numpy as np
from numpy.fft import fftshift, ifftshift, ifft, fft

def fwhm(lst):
    '''Calculate the fwhm given XXX?'''
    #not a requirement that lst is real?
    if isinstance(lst, complex):
        lst = abs(lst)**2

    halfmax = max(lst)/2
    fwhm    = halfmax[np.where(halfmax < lst)]
    return fwhm

def fft(field):
    '''fft with 0 shift?'''
    #https://numpy.org/doc/stable/reference/generated/numpy.fft.ifftshift.html
    return ifftshift(fft(fftshift(field)))

def ifft(field):
    '''inverse fft with 0 shift?'''
    #https://numpy.org/doc/stable/reference/generated/numpy.fft.fftshift.html
    return ifftshift(ifft(fftshift(field)))

def inten_peak(epulse, mrad, tau):
    '''Calculate peak intensity?'''
    #double check log type
    numerator   = 2*epulse*np.sqrt(np.log(16))
    denominator = (np.pi**(3/2))*(mrad**2)*tau 
    return numerator/denominator

## --- Put these into a class? ----
def intenF(self,field,n):
    ''' '''
    # not done
    #np.real( (2*self.const_eps0*n*self.const_c) .* conj(field).*field )
    return

def energyf(self, field,n,mRad):
    '''Calculate energy of XXX?'''
    #not done
    const = (np.pi/2)*mRad**2
    return #(const * sum(self.intenF(this,field,n)) * self.grid_dt;


