import numpy as np
from load import load_data
from convert import convert
from pmc2 import *

def error_computation(P1,P2):

    ## Error based on Asig1+Bsig2+Csig3=1
    P1_coeff = [P1.A,P1.B,P1.C]
    P2_coeff = [P2.A,P2.B,P2.C]
    sP1 = np.dot(P1.data[:,:3],P1_coeff) #yi P1
    sP2 = np.dot(P2.data[:,:3],P2_coeff) #yi P2
    sP = np.concatenate((sP1,sP2),axis=None) #yi
    nb_err = sP.size
    err_P1 = sP1-1
    err_P2 = sP2-1
    err = sP-1
    sum_err_P1 = np.sum(abs(err_P1))
    mean_err_P1 = np.mean(abs(err_P1))
    sum_err_P2 = np.sum(abs(err_P2))
    mean_err_P2 = np.mean(abs(err_P2))
    sum_err = np.sum(abs(err))
    mean_err = np.mean(abs(err))
    
    # R^2 = 1 - SSres/SStot
    #SSres = Sum(yi-fi)^2
    
    SSres = np.sum(np.square(err))
    SSres_P1 = np.sum(np.square(err_P1))
    SSres_P2 = np.sum(np.square(err_P2))
        
    # SStot = Sum(yi-ymean)^2
    ymean = np.mean(sP)
    SStot = np.sum(np.square(sP-ymean))
    ymean_P1 = np.mean(sP1)
    SStot_P1 = np.sum(np.square(sP1-ymean))
    ymean_P2 = np.mean(sP2)
    SStot_P2 = np.sum(np.square(sP2-ymean))
        
    # R^2 = 1 - SSres/SStot
    Rsqu = 1- (SSres/SStot)
    Rsqu_P1 = 1 - SSres_P1/SStot_P1
    Rsqu_P2 = 1 - SSres_P2/SStot_P2

    return mean_err, mean_err_P1, mean_err_P2
    #return Rsqu, Rsqu_P1, Rsqu_P2