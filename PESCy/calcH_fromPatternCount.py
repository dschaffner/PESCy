import numpy as np
from collections import Counter

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
    S - 

    '''
    
    T=np.array(data)
    if len(T.shape)>1:
        raise TypeError('Data must be a 1-D array')
    t = len(T)
    Ptot = t - delay*(n - 1)    #Total number of order n permutations in T
    #print 'Number of permutations = ', Ptot
    A = []			 #Array to store each permutation
    
    for i in range(Ptot):	#Will run through all possible n segments of T
        A.append(''.join(T[i:i+(n-1)*delay+1:delay].argsort().astype(str)))
    #Count occurance of patterns
    count=Counter(A)
    return count,Ptot