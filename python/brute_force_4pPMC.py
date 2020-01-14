import numpy as np
from load import load_data
from convert import convert
from pmc_4p import *
from mc import *
from hb import *
from error_computation import error_computation
from error_computation_4pPMC import *

def brute_force_one_set(data_a):
    bp_set = []
    mean_error = []
    mean_error_P1 = []
    mean_error_P2 = []
    
    for i in tqdm(range(data_a.shape[0]-1,0,-1)):
        bp_set.append([i])
        P1_a = data_a[i:,:]
        P2_a = data_a[:i,:]
        P1, P2 = create_P1_and_P2(P1_a,P2_a)
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
    Err = mean_error[min_error_index]
    Err_P1 = mean_error_P1[min_error_index]
    Err_P2 = mean_error_P2[min_error_index]
    
    return P1_a, P2_a

def allocation(P1,P2,data_C,data_E,data_o):
    P1_data = []
    P2_data = []
    data = np.concatenate([data_C,data_E,data_o])
    for i in range(data.shape[0]):
        err_P1, err_P2 = error_computation_4p(P1,P2,data[i,:])
        if err_P1 < err_P2:
            P1_data.append(data[i,:])
        elif err_P1 > err_P2:
            P2_data.append(data[i,:])

    return np.array(P1_data),np.array(P2_data)

def planes_def_PMC(data,fit_pts):
    
    ## INITIALIZATION
    data_C = np.concatenate([data[0:1,:],data[4:5,:],data[8:24,:]])
    data_C = data_C[np.argsort(data_C[:, 3])] # sorted with p (mean stress)
    data_E = get_E(data)[np.argsort(get_E(data)[:, 3])] # sorted with p (mean stress)
    data_o = get_o(data)[np.argsort(get_o(data)[:, 5])] # sorted with teta (Lode angle)
    
    ## ITERATIONS
    P1_C, P2_C = brute_force_one_set(fit_pts) 
    ## RESULTS
    P1_init, P2_init = create_P1_and_P2(P1_C,P2_C)
    if P1_init.Vo < P2_init.Vo :
        P1 = P2_init
        P2 = P1_init
    elif P1_init.Vo > P2_init.Vo :
        P1 = P1_init
        P2 = P2_init 
    print(P1.sol)
    print(P2.sol)
    
    ## ALLOCATIONS OF E AND o DATA TO P1 OR P2
    P1_data, P2_data = allocation(P1,P2,data_C,data_E,data_o)
    print(P1_data)
    print('o')
    print(P2_data)
    S_P1 = standard_dev_4p(P1,np.array(P1_data),'PMC')
    S_P2 = standard_dev_4p(P2,np.array(P2_data),'PMC')
    
    
    return P1, P2, S_P1, S_P2

def planes_def_MC_HB(data,fit_pts,criterion):
    
    ## INITIALIZATION
    data_C = np.concatenate([data[0:1,:],data[4:5,:],data[8:24,:]])
    data_C = data_C[np.argsort(data_C[:, 3])] # sorted with p (mean stress)
    data_E = get_E(data)[np.argsort(get_E(data)[:, 3])] # sorted with p (mean stress)
    data_o = get_o(data)[np.argsort(get_o(data)[:, 5])] # sorted with teta (Lode angle)
    
    ## RESULTS
    if criterion == 'MC':
        P= Plane_MC(fit_pts)
    elif criterion =='HB':
        P= Plane_HB(fit_pts)
    print(P.sol)
    
    non_fit_data = np.concatenate((data_C,data_E,data_o))
    S = standard_dev_4p(P,non_fit_data,criterion)
    
    
    return P, S