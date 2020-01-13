import numpy as np
from load import load_data
from convert import convert
from pmc import *

def error_computation(P1,P2):
    ## Error based on AsigI+BsigII+CsigIII=1
    P1_coeff = [P1.A,P1.B,P1.C]
    P2_coeff = [P2.A,P2.B,P2.C]
    sP1 = np.dot(P1.data[:,:3],P1_coeff) #yi P1
    sP2 = np.dot(P2.data[:,:3],P2_coeff) #yi P2
    sP = np.concatenate((sP1,sP2),axis=None) #yi
    nb_err = sP.size
    err_P1 = sP1-1
    err_P2 = sP2-1
    err = sP-1
    mean_err_P1 = np.mean(np.square(err_P1))
    mean_err_P2 = np.mean(np.square(err_P2))
    mean_err = np.mean(np.square(err))
    
    return mean_err, mean_err_P1, mean_err_P2

def standard_dev(P,data,criterion):
    #Standard deviation
    d = convert(data)
    m = np.unique(d['C'][:,1]).size + np.unique(d['E'][:,1]).size + np.unique(d['o'][:,1]).size
    conf = np.unique(data[:,1])
    mconf = conf.size
    nC = 0
    setC = 0
    nE = 0
    setE = 0
    no = 0
    seto = 0
    Si_int = 0
    for i in range(mconf):
        for j in range(data[:,0].size):
            # Standard deviation of COMPRESSION data
            if data[j,1]==conf[i] and data[j,5] == 0:
                nC = nC+1
                sig_test = data[j,0]
                sig_calc = (1/P.A)*(-P.B*data[j,1]-P.C*data[j,2]+1)
            # Standard deviation of EXTENSION data
            elif data[j,1]==conf[i] and data[j,5] == 60:
                nE = nE+1
                sig_test = data[j,0]
                sig_calc = (1/P.A)*(-P.B*P.data[j,1]-P.C*P.data[j,2]+1)
            # Standard deviation of MULTI AXIAL data
            elif P.data[j,1]==conf[i] and P.data[j,5] != 0 and P.data[j,5] != 60:
                no = no+1
                sig_test = P.data[j,0]
                sig_calc = (1/P.A)*(-P.B*P.data[j,1]-P.C*P.data[j,2]+1)
                sig_calc = P.data[j,2]+P.Co*np.sqrt((P.m/P.Co)*P.data[j,2]+1)
                setE = setE + np.square(sig_test-sig_calc)
        
        if nC ==0:
            sC=0 
        else:
            sC = np.sqrt((1/nC)*setC)
        if nE ==0:
            sE=0 
        else:
            sE = np.sqrt((1/nE)*setE)
        if no ==0:
            so=0 
        else:
            so = np.sqrt((1/no)*seto)
    
        Si_int = Si_int + sC+sE+so
        
    S = Si_int/m
    
    return S