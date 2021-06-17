import numpy as np
from math import factorial

def calcS_fromPatternCount(count,tot_perms,n):
    '''
    

    Parameters
    ----------
    count : A Counter object
        Count occurance result from constructPatternCount()
    tot_perms : integer
        total number of permutations from constructPatternCount()
    n : integer
        embedding dimension from constructPatternCount()

    Returns
    -------
    S - shannon permutation entropy
    Se - Shannon + Uniform permuation Entropy

    '''
    Ptot=tot_perms
    N=factorial(n)
    invPtot=1./Ptot     #Inverse for later calcuations
    S = 0.
    Se = 0.
    for q in iter(count.values()):
        q*=invPtot #convert to probability
        S += -q * np.log2(q)
        q+=1./N
        q/=2
        Se += -q * np.log2(q)
    for i in range(len(count),N):
        q=1./2/N
        Se += -q * np.log2(q)
    return S,Se