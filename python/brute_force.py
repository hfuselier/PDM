import numpy as np
from load import load_data
from convert import convert
from pmc import *
from error_computation import *

def brute_force_one_set(data_a,data_b,data_c):
    bp_set = []
    P1_b = data_b
    P2_b = data_b
    P1_c = data_c
    P2_c = data_c
    mean_error = []
    mean_error_P1 = []
    mean_error_P2 = []
    
    for i in tqdm(range(data_a.shape[0]-1,0,-1)):
        bp_set.append([i])
        P1_a = data_a[i:,:]
        P2_a = data_a[:i,:]
        P1, P2 = create_P1_and_P2(P1_a,P2_a,P1_b,P2_b,P1_c,P2_c)
        mean_err, mean_err_P1, mean_err_P2 = error_computation(P1,P2)
        mean_error.append(mean_err)
        mean_error_P1.append(mean_err_P1)
        mean_error_P2.append(mean_err_P2)

    mean_error = np.array(mean_error)
    mean_error_P1 = np.array(mean_error_P1)
    mean_error_P2 = np.array(mean_error_P2)
    bp_set = np.array(bp_set)
    min_error_index = np.nanargmin(mean_error)
    bpa=bp_set[min_error_index]
    P1_a=data_a[int(bpa):,:]
    P2_a=data_a[:int(bpa),:]
    #Err = mean_error[min_error_index]
    #Err_P1 = mean_error_P1[min_error_index]
    #Err_P2 = mean_error_P2[min_error_index]
    
    return P1_a, P2_a, P1_b, P2_b, P1_c, P2_c

def brute_force_two_sets(data_a,data_b,data_c):
    bp_set = []
    P1_c = data_c
    P2_c = data_c
    mean_error = []
    mean_error_P1 = []
    mean_error_P2 = []
    
    for i in tqdm(range(data_a.shape[0]-1,0,-1)):
        for j in range(data_b.shape[0]-1,0,-1):
            bp_set.append([i,j])
            P1_a = data_a[i:,:]
            P2_a = data_a[:i,:]
            P1_b = data_b[j:,:]
            P2_b = data_b[:j,:]
            P1, P2 = create_P1_and_P2(P1_a,P2_a,P1_b,P2_b,P1_c,P2_c)
            mean_err, mean_err_P1, mean_err_P2 = error_computation(P1,P2)
            mean_error.append(mean_err)
            mean_error_P1.append(mean_err_P1)
            mean_error_P2.append(mean_err_P2)

    mean_error = np.array(mean_error)
    mean_error_P1 = np.array(mean_error_P1)
    mean_error_P2 = np.array(mean_error_P2)
    bp_set = np.array(bp_set)
    min_error_index = np.nanargmin(mean_error)
    bpa, bpb=bp_set[min_error_index]
    P1_a=data_a[int(bpa):,:]
    P2_a=data_a[:int(bpa),:]
    P1_b=data_b[int(bpb):,:]
    P2_b=data_b[:int(bpb),:]
    #Err = mean_error[min_error_index]
    #Err_P1 = mean_error_P1[min_error_index]
    #Err_P2 = mean_error_P2[min_error_index]
    
    return P1_a, P2_a, P1_b, P2_b, P1_c, P2_c #Err,Err_P1, Err_P2, 

def brute_force_three_sets(data_a,data_b,data_c):
    bp_set = [] 
    mean_error = []
    mean_error_P1 = []
    mean_error_P2 = []
    
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
                P1, P2 = create_P1_and_P2(P1_a,P2_a,P1_b,P2_b,P1_c,P2_c) 
                mean_err, mean_err_P1, mean_err_P2 = error_computation(P1,P2)
                mean_error.append(mean_err)
                mean_error_P1.append(mean_err_P1)
                mean_error_P2.append(mean_err_P2)
        
    mean_error = np.array(mean_error)
    mean_error_P1 = np.array(mean_error_P1)
    mean_error_P2 = np.array(mean_error_P2)
    min_error_index = np.nanargmin(mean_error)
    bp_set = np.array(bp_set)
    bpa, bpb, bpc = bp_set[min_error_index]
    P1_a=data_a[int(bpa):,:]
    P2_a=data_a[:int(bpa),:]
    P1_b=data_b[int(bpb):,:]
    P2_b=data_b[:int(bpb),:]
    P1_c=data_c[int(bpc):,:]
    P2_c=data_c[:int(bpc),:]
    #Err = mean_error[min_error_index]
    #Err_P1 = mean_error_P1[min_error_index]
    #Err_P2 = mean_error_P2[min_error_index]
    
    return P1_a, P2_a, P1_b, P2_b, P1_c, P2_c #Err, Err_P1, Err_P2, 

def planes_def(data,d):
    
    lin_C = d['pC'].size
    lin_E = d['pE'].size
    lin_o = d['po'].size
    
    ## INITIALIZATION
    data_C = get_C(data)[np.argsort(get_C(data)[:, 3])] # sorted with p (mean stress)
    data_E = get_E(data)[np.argsort(get_E(data)[:, 3])] # sorted with p (mean stress)
    data_o = get_o(data)[np.argsort(get_o(data)[:, 5])] # sorted with teta (Lode angle)
    
    ## ITERATIONS
    if lin_C != 0 and lin_E == 0 and lin_o == 0:
        P1_C, P2_C, P1_E, P2_E, P1_o, P2_o = brute_force_one_set(data_C,data_E,data_o)
    if lin_C == 0 and lin_E != 0 and lin_o != 0:
        P1_E, P2_E, P1_o, P2_o, P1_C, P2_C = brute_force_two_sets(data_E,data_o,data_C)
    if lin_E == 0 and lin_C != 0 and lin_o != 0:
        P1_C, P2_C, P1_o, P2_o, P1_E, P2_E = brute_force_two_sets(data_C,data_o,data_E)
    if lin_o == 0 and lin_C != 0 and lin_E != 0:
        P1_C, P2_C, P1_E, P2_E, P1_o, P2_o = brute_force_two_sets(data_C,data_E,data_o)
    if lin_C != 0 and lin_E != 0 and lin_o != 0 :
        P1_C, P2_C, P1_E, P2_E, P1_o, P2_o = brute_force_three_sets(data_C,data_E,data_o)
        
    
    ## RESULTS
    P1_init, P2_init = create_P1_and_P2(P1_C,P2_C,P1_E,P2_E,P1_o,P2_o)

    if P1_init.Vo < P2_init.Vo :
        P1 = P2_init
        P2 = P1_init
    elif P1_init.Vo > P2_init.Vo :
        P1 = P1_init
        P2 = P2_init
    
    S_P1 = standard_dev(P1,'PMC')
    S_P2 = standard_dev(P2,'PMC')
    
    
    return P1, P2, S_P1, S_P2