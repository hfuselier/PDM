import numpy as np
from load import load_data

def convert(data):
    row = data.shape[0]

    p = data[:,3]
    q = data[:,4]


    d = {}

    d['E'] = data[data[:,5] == 60]
    d['C'] = data[data[:,5] == 0]
    d['o'] = data[np.logical_and(data[:,5] != 60 , data[:,5] != 0)]
    d['p'] = data[:,3]
    d['q'] = data[:,4]


    
    for G in ('E','C','o') :
        d['q'+ G] = d[G][:,4]
        d['p'+ G] = d[G][:,3]
        d['t'+ G] = d[G][:,5]
    for i in range(3):
        d['s'+str(i+1)+'o'] = d['o'][:,i]
        
    
    for i in range(1,4):
        d['sig'+str(i)] = data[:,i-1]
        
    d['conf'] = np.unique(d['sig3'])
        
    return d







