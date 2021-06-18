import numpy as np
from math import factorial
from collections import Counter
from constructPatternCount import constructPatternCount
from calcS_fromPatternCount import calcS_fromPatternCount

def calcPESCcurves(data,n=5,max_delay=100):
    '''
    function calcPESCcurves - Returns PE(tau) and SC(tau)
    Input:
        data  - array
        n     - embedding dimension (default=5)
        max_delay - largest embedding delays to loop through (default=100)
    Output:
        C(tau) - Normalized Jensen-Shannon complexity as function of embedding delay tau
        H - Normalized Shannon Perumation Entropy as function of embedding delay tau
    '''
    nfac = factorial(n)
    delay_array = np.arange(1,max_delay)
    num_delays=len(delay_array)
    PEs=np.zeros([num_delays])
    SCs=np.zeros([num_delays])
    for loop_delay in np.arange(len(delay_array)):		
        if (loop_delay%100)==0: print( 'On Delay ',delay_array[loop_delay])
        permstore_counter = []
        permstore_counter = Counter(permstore_counter)
        tot_perms = 0
        arr,nperms,n0 = constructPatternCount(data,n=5,delay=delay_array[loop_delay])
        permstore_counter = permstore_counter+arr
        tot_perms = tot_perms+nperms
        PE_tot,PE_tot_Se = calcS_fromPatternCount(permstore_counter,tot_perms,n0)
        C =  -2.*((PE_tot_Se - 0.5*PE_tot - 0.5*np.log2(nfac))
                    /((1 + 1./nfac)*np.log2(nfac+1) - 2*np.log2(2*nfac) 
                    + np.log2(nfac))*(PE_tot/np.log2(nfac)))
        PEs[loop_delay]=PE_tot/np.log2(nfac)
        SCs[loop_delay]=C
    return PEs,SCs