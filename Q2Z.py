##############################################################################
# Created by Giovanni Pireddu
# giovanni.pireddu(at)sorbonne-universite.fr
#
# Script to compute admittance / impedance from the time series of the charge
# Assuming a MetalWalls format as described in:
# https://gitlab.com/ampere2/metalwalls/-/wikis/output-files#total_charges.out
#
##############################################################################


import numpy as np
from scipy.interpolate import lagrange
from numba import njit, prange

import timeit

complextype= complex

@njit(parallel=True)
def FilonLagrange(DFT,freq,Time,coeffs):
    """
    Parameters
    ----------
    DFT : complex
        Array result of the Fourier-Laplace transform
    freq : float
        Array of frequencies (rad/s)
    Time : float
        Array of time points (s)
    coeffs : float
        Coefficients from the Lagrange interpolation

    Returns
    -------
    DFT : complex
        Array result of the Fourier-Laplace transform

    """
    for ifreq in prange(len(freq)):
        ip= 0
        for it in range(0,len(Time)-2,2): #Skip intermediate points
            t1= Time[it]
            t2= Time[it+2]
            a= coeffs[ip,0]
            b= coeffs[ip,1]
            c= coeffs[ip,2]
            ff= freq[ifreq]
            DFT[ifreq] += (1/ff**3) * (np.exp(-1j*ff*t2) * (a * (ff * t2 * (2 + 1j * ff * t2) - 2 * 1j) 
                             + ff * (1j * b * ff * t2 + b + 1j * c * ff)) - 
                             1j * np.exp(-1j * ff * t1) * (a * (-2 + ff * t1 * (ff * t1 - 2*1j)) +
                                                           ff * (c * ff + b * (ff * t1 - 1j))))
            ip += 1
            
    return DFT

def LagrangeInterpol(coeffs,Time,Signal):
    """   
    Parameters
    ----------
    coeffs : float
        Coefficients from the Lagrange interpolation
    Time : float
        Array of time points (s)
    Signal : TYPE
        DESCRIPTION.

    Returns
    -------
    coeffs : float
        Coefficients from the Lagrange interpolation

    """
    ipp= 0
    for ip in range(1,len(Time)-1,2):
        x= Time[ip-1:ip+2]
        y= Signal[ip-1:ip+2]
        p= lagrange(x,y)   
        cc= p.coef       
        coeffs[ipp,:] = cc
        ipp += 1
        
    return np.array(coeffs)
  
def AdmFromQ(Time,QACF,freq):
    """
    Parameters
    ----------
    Time : float
        Array of time points (s)
    QACF : float
        Total charge autocorrelation function
    freq : float
        Array of frequencies (rad/s)

    Returns
    -------
    Adm : complex
        Admittance array (1/ohm)

    """
    coeffs= np.zeros(((len(Time))//2,3))
    Coeffs= LagrangeInterpol(coeffs,Time,QACF)
    DFT= FilonLagrange(np.zeros(len(freq),dtype=complex),freq,Time,Coeffs)
    Adm= np.zeros(len(DFT),dtype=complex)
    for i in range(len(freq)):
        Adm[i]= beta * ((freq[i]**2)*DFT[i] + 1j * freq[i] * QACF[0] )
    return Adm

def WKACF(x,time):
    """
    Parameters
    ----------
    x : float
        Signal
    time : float
        Array of time points (s)
        
    Returns
    -------
    ACF : float
        Autocorrelation function
    time : float
        Array of time points (s) [Half the length of the original array]

    """
    Ctt = np.fft.fftn(x-np.mean(x))
    CC = Ctt[:] * np.conjugate(Ctt[:])
    CC[:] = np.fft.ifftn(CC[:])
    ACF = (CC[:len(CC)//2]).real / len(CC)
    time= time[:len(time)//2]
    return ACF, time

start = timeit.default_timer()

print('Start')

# Parameters
temperature= 298 #K
timeperstep= 1E-15 #s
nfreq= 100

# Constants
k = 1.3806485279e-23 #J.K-1
e = 1.602176620898e-19 #C
beta= 1/(k*temperature)

# Load data
Data= np.loadtxt('total_charges.out',skiprows=3)

print('Data loaded')

# Convert in SI units
Charges= Data[:,1] * e
Time= np.arange(0,len(Data),1) * 1E-15

# QACF
QACF, redTime = WKACF(Charges,Time)

np.savetxt('QACF.out',np.column_stack((redTime,QACF)),header='Time (s) / QACF (C^2)')

print('QACF computed')

# List of frequencies
hif= (2*np.pi)/(redTime[1]*3)
lof= (2*np.pi)/redTime[-1]
freq= np.logspace(np.log10(lof),np.log10(hif),nfreq)

# Admittance / Impedance
Adm= AdmFromQ(redTime,QACF, freq)
Imp= 1/Adm

print('Admittance / Impedance computed')

np.savetxt('Admittance.out',np.column_stack((freq,Adm.real,Adm.imag)),
           header='Frequency (rad/s) / Re[Y] (1/Ohm) / Im[Y] (1/Ohm)')
np.savetxt('Impedance.out',np.column_stack((freq,Adm.real,Adm.imag)),
           header='Frequency (rad/s) / Re[Z] (Ohm) / Im[Z] (Ohm)')

stop = timeit.default_timer()

print('Time: ', stop - start)  

print('Done.')





