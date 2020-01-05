import numpy as np
from load import load_data
from convert import convert
from pmc2 import *

def brute_force_two_sets(data_a,data_b,data_c):
    bp_set = []
    P1_c = data_c
    P2_c = data_c
    mean_error = []
    
    for i in tqdm(range(data_a.shape[0]-1,0,-1)):
        for j in range(data_b.shape[0]-1,0,-1):
            bp_set.append([i,j])
            P1_a = data_a[i:,:]
            P2_a = data_a[:i,:]
            P1_b = data_b[j:,:]
            P2_b = data_b[:j,:]
            P1, P2 = create_P1_and_P2_iteration(P1_a,P2_a,P1_b,P2_b,P1_c,P2_c)
            #if P1.phyC > P1.phyE or P2.phyC > P2.phyE:
            #    break
            P1_coeff = [P1.A,P1.B,P1.C]
            P2_coeff = [P2.A,P2.B,P2.C]
            sP1 = np.dot(np.transpose(P1.sig123),P1_coeff)
            sP2 = np.dot(np.transpose(P2.sig123),P2_coeff)
            err_P1 = abs(1-sP1)
            err_P2 = abs(1-sP2)
            sumP1 = np.sum(err_P1)
            sumP2 = np.sum(err_P2)
            mean_error.append((sumP1 + sumP2)/(len(err_P1)+len(err_P2)))

    mean_error = np.array(mean_error)
    bp_set = np.array(bp_set)
    min_error_index = np.nanargmin(mean_error)
    bpa, bpb=bp_set[min_error_index]
    P1_a=data_a[int(bpa):,:]
    P2_a=data_a[:int(bpa),:]
    P1_b=data_b[int(bpb):,:]
    P2_b=data_b[:int(bpb),:]
    
    return mean_error, P1_a, P2_a, P1_b, P2_b, P1_c, P2_c

def brute_force_three_sets(data_a,data_b,data_c):
    bp_set = [] 
    mean_error = []
    a=0
    
    for i in tqdm(range(data_a.shape[0]-1,0,-1)) :
        for j in range(data_b.shape[0]-1,0,-1) :
            for k in range(data_c.shape[0]-1,0,-1) :
                bp_set.append([i,j,k])
                P1_a=data_a[i:,:]
                P2_a=data_a[:i,:]
                P1_b=data_b[j:,:]
                P2_b=data_b[:j,:]
                P1_c=data_c[k:,:]
                P2_c=data_c[:k,:]
                P1, P2 = create_P1_and_P2_iteration(P1_a,P2_a,P1_b,P2_b,P1_c,P2_c) 
                #if P1.phyC > P1.phyE or P2.phyC > P2.phyE:
                #    break
                P1_coeff = [P1.A,P1.B,P1.C]
                P2_coeff = [P2.A,P2.B,P2.C]
                sP1 = np.dot(np.transpose(P1.sig123),P1_coeff)
                sP2 = np.dot(np.transpose(P2.sig123),P2_coeff)
                err_P1 = abs(1-sP1)
                err_P2 = abs(1-sP2)
                sumP1 = np.sum(err_P1)
                sumP2 = np.sum(err_P2)
                mean_error.append((sumP1 + sumP2)/(len(err_P1)+len(err_P2)))
                a=a+1
        
    mean_error = np.array(mean_error)
    min_error_index = np.nanargmin(mean_error)
    bpa, bpb, bpc = bp_set[min_error_index]
    P1_a=data_a[int(bpa):,:]
    P2_a=data_a[:int(bpa),:]
    P1_b=data_b[int(bpb):,:]
    P2_b=data_b[:int(bpb),:]
    P1_c=data_c[int(bpc):,:]
    P2_c=data_c[:int(bpc),:]
    
    return mean_error, P1_a, P2_a, P1_b, P2_b, P1_c, P2_c

def planes_def(data,d):
    
    lin_C = d['pC'].size
    lin_E = d['pE'].size
    lin_o = d['po'].size
    
    ## INITIALIZATION
    data_C = get_C(data)[np.argsort(get_C(data)[:, 3])] # sorted with p (mean stress)
    data_E = get_E(data)[np.argsort(get_E(data)[:, 3])] # sorted with p (mean stress)
    data_o = get_o(data)[np.argsort(get_o(data)[:, 5])] # sorted with teta (Lode angle)
    
    ## ITERATIONS
    if lin_C == 0 :
        mean_error, P1_E, P2_E, P1_o, P2_o, P1_C, P2_C = brute_force_two_sets(data_E,data_o,data_C)
    if lin_E == 0 :
        mean_error, P1_C, P2_C, P1_o, P2_o, P1_E, P2_E = brute_force_two_sets(data_C,data_o,data_E)
    if lin_o == 0 :
        mean_error, P1_C, P2_C, P1_E, P2_E, P1_o, P2_o = brute_force_two_sets(data_C,data_E,data_o)
    if lin_C != 0 and lin_E != 0 and lin_o != 0 :
        mean_error, P1_C, P2_C, P1_E, P2_E, P1_o, P2_o = brute_force_three_sets(data_C,data_E,data_o)
        
    
    ## RESULTS
    P1_init, P2_init = create_P1_and_P2_iteration(P1_C,P2_C,P1_E,P2_E,P1_o,P2_o)

    if P1_init.Vo < P2_init.Vo :
        P1 = P2_init
        P2 = P1_init
    elif P1_init.Vo > P2_init.Vo :
        P1 = P1_init
        P2 = P2_init
    
    print(P1.sol)
    print(P2.sol)
    
    ## ERROR COMPUTATION
    ## Error based on Asig1+Bsig2+Csig3=1
    P1_coeff = [P1.A,P1.B,P1.C]
    P2_coeff = [P2.A,P2.B,P2.C]
    sP1 = np.dot(np.transpose(P1.sig123),P1_coeff) #yi P1
    sP2 = np.dot(np.transpose(P2.sig123),P2_coeff) #yi P2
    sP = np.concatenate((sP1,sP2),axis=None) #yi
    nb_err = sP.size
    err_P1 = abs(1-sP1)
    err_P2 = abs(1-sP2)
    err = abs(1-sP)
    f_P1 = 1*np.ones(err_P1.shape)
    f_P2 = 1*np.ones(err_P1.shape)
    f = np.concatenate((f_P1,f_P2),axis=None)
    
    # R^2 = 1 - SSres/SStot
    #SSres = Sum(yi-fi)^2
    nb_err = err.size
    nb_err_P1 = err_P1.size
    nb_err_P2 = err_P2.size
    
    SSres = 0
    SSres_P1 = 0
    SSres_P2 = 0
    for i in range(0,nb_err-1,1) :
        SSres = SSres + (err[i])**2
    for j in range(0,nb_err_P1-1,1) :
        SSres_P1 = SSres_P1 + (err_P1[j])**2
    for k in range(0,nb_err_P2-1,1) :
        SSres_P2 = SSres_P2 + (err_P2[k])**2
        
    # SStot = Sum(yi-ymean)^2
    ymean = np.mean(abs(err))
    SStot = 0
    ymean_P1 = np.mean(abs(err_P1))
    SStot_P1 = 0
    ymean_P2 = np.mean(abs(err_P2))
    SStot_P2 = 0
    for i in range(0,nb_err-1,1) :
        SStot = SStot + (f[i]-ymean)**2
    for j in range(0,nb_err_P1-1,1) :
        SStot_P1 = SStot_P1 + (f_P1[j]-ymean_P1)**2
    for k in range(0,nb_err_P2-1,1) :
        SStot_P2 = SStot_P2 + (f_P2[k]-ymean_P2)**2
        
    # R^2 = 1 - SSres/SStot
    Rsqu = 1 - SSres/SStot
    Rsqu_P1 = 1 - SSres_P1/SStot_P1
    Rsqu_P2 = 1 - SSres_P2/SStot_P2

    return P1, P2, Rsqu, Rsqu_P1, Rsqu_P2