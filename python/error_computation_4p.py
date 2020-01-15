import numpy as np
from load import load_data
from convert import convert
from pmc_4p import *

def error_computation_4p(P1,P2,data):
    ## Error based on AsigI+BsigII+CsigIII=1
    P1_coeff = [P1.A,P1.B,P1.C]
    P2_coeff = [P2.A,P2.B,P2.C]
    sP1 = np.dot(data[:3],P1_coeff) #yi P1
    sP2 = np.dot(data[:3],P2_coeff) #yi P2
    err_P1 = abs(sP1-1)/np.sqrt(P1.A**2+P1.B**2+P1.C**2)
    err_P2 = abs(sP2-1)/np.sqrt(P2.A**2+P2.B**2+P2.C**2)
    
    return err_P1, err_P2

def standard_dev_4p(P,data,criterion):
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
                # Mohr-Coulomb and Paul-Mohr-Coulomb criteria
                if criterion == 'MC' or criterion == 'PMC':
                    sig_calc = (1/P.A)*(-P.B*data[j,1]-P.C*data[j,2]+1)
                # Hoek-Brown criterion
                elif criterion == 'HB':
                    sig_calc = data[j,2]+P.Co*np.sqrt((P.m/P.Co)*data[j,2]+1)
                setC = setC + np.square(sig_test-sig_calc)
            # Standard deviation of EXTENSION data
            elif data[j,1]==conf[i] and data[j,5] == 60:
                nE = nE+1
                sig_test = data[j,0]
                # Mohr-Coulomb and Paul-Mohr-Coulomb criteria
                if criterion == 'MC' or criterion == 'PMC':
                    sig_calc = (1/P.A)*(-P.B*data[j,1]-P.C*data[j,2]+1)
                # Hoek-Brown criterion
                elif criterion == 'HB':
                    sig_calc = data[j,2]+P.Co*np.sqrt((P.m/P.Co)*data[j,2]+1)
                setE = setE + np.square(sig_test-sig_calc)
            # Standard deviation of MULTI AXIAL data
            elif data[j,1]==conf[i] and data[j,5] != 0 and data[j,5] != 60:
                no = no+1
                sig_test = data[j,0]
                # Mohr-Coulomb and Paul-Mohr-Coulomb criteria
                if criterion == 'MC' or criterion == 'PMC':
                    sig_calc = (1/P.A)*(-P.B*data[j,1]-P.C*data[j,2]+1)
                # Hoek-Brown criterion
                elif criterion == 'HB':
                    sig_calc = data[j,2]+P.Co*np.sqrt((P.m/P.Co)*data[j,2]+1)
                seto = seto + np.square(sig_test-sig_calc)
        
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