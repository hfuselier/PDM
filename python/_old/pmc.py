import numpy as np
from numpy import sin, tan ,cos, pi, sqrt
import math

def get_C(data):
    return data[data[:,5] == 0]
def get_E(data):
    return data[data[:,5] == 60]
def get_o(data):
    return data[np.logical_and(data[:,5] != 60 , data[:,5] != 0)]

def get_p(data):
    return data[:,3]
def get_q(data):
    return data[:,4]
def get_t(data):
    return data[:,5]


def create_P1_and_P2(data,d,bpC,bpE):
    P1 = {}
    P2 = {}
    
    P1_C = get_C(data)[get_p(get_C(data)) >= bpC]
    P1_E = get_E(data)[get_p(get_E(data)) >= bpE]
    P1_o = get_o(data)[get_o(data)[:,0] < 2*get_o(data)[:,0]]

    P2_C = get_C(data)[get_p(get_C(data)) < bpC]
    P2_E = get_E(data)[get_p(get_E(data)) < bpE]
    P2_o = get_o(data)[get_o(data)[:,0] >= 2*get_o(data)[:,0]]

    P1 = {}
    P2 = {}
    P1['data'] = np.concatenate((P1_C,P1_E,P1_o),axis=0)
    P2['data'] = np.concatenate((P2_C,P2_E,P2_o),axis=0)
    
    # sperating datas into 2 planes the planes
    for i in ('p','q','t'):
        key = i+'C'
        P1[key] = d[key][d['pC'] >= bpC]
        P2[key] = d[key][d['pC'] <= bpC]

        key = i+'E'
        P1[key] = d[key][d['pE'] >= bpE]
        P2[key] = d[key][d['pE'] <= bpE]

        key = i+'o'
        P1[key] = d[key][d['s1o'] <  2*d['s2o']]
        P2[key] = d[key][d['s1o'] >= 2*d['s2o']]

        P1[i+'tot'] = np.concatenate((P1[i+'C'],P1[i+'E'],P1[i+'o']))
        P2[i+'tot'] = np.concatenate((P2[i+'C'],P2[i+'E'],P2[i+'o']))
    
    P1['ttot']*=pi/180
    P2['ttot']*=pi/180
    
    #solving plane equations
    for P in (P1,P2):
        a = P['qtot']*np.sin(P['ttot'])

        A = np.ndarray((P['ptot'].shape[0],3))
        A[:,0] = P['ptot']
        A[:,1] = a
        A[:,2] = 1

        b = P['qtot']*np.cos(P['ttot'])

        x = np.linalg.pinv(A)@b
        print(x)
        
        bc = x[2]
        k = x[1]
        P['Vo'] = bc/x[0]
        be = 2*bc/(1-sqrt(3)*k)
        P['phyC'] = math.asin(3*bc/(6*P['Vo']+bc))
        print(3*be/(6*P['Vo']-be))
        P['phyE'] = math.asin(3*be/(6*P['Vo']-be))
        P['Sol'] = np.array([bc, be, k, P['Vo'], P['phyC']*180/pi, P['phyE']*180/pi])
        
        P['Nc'] = (1-sin(P['phyC']))/(2*sin(P['phyC']))
        P['Ne'] = (1-sin(P['phyE']))/(2*sin(P['phyE']))
        
        P['A'] = ((1-sin(P['phyC']))/(2*sin(P['phyC'])))/P['Vo']
        P['B'] = ((sin(P['phyC'])-sin(P['phyE']))/(2*sin(P['phyC'])*sin(P['phyE'])))/P['Vo']
        P['C'] = -((1+sin(P['phyE']))/(2*sin(P['phyE'])))/P['Vo']
        
    #if P1['Vo'] < P2['Vo']:
     #   P1, P2 = P2, P1
        
    return P1, P2

def qCfit(P,x):
    phyCP = P['phyC']
    phyEP = P['phyE']
    VoP = P['Vo']
    return (6*sin(phyCP)/(3-sin(phyCP)))*x+(6*VoP*tan(phyCP)*cos(phyCP)/(3-sin(phyCP)));

def qEfit(P,x):
    phyCP = P['phyC']
    phyEP = P['phyE']
    VoP = P['Vo']
    return -(6*sin(phyEP)/(3+sin(phyEP)))*x-(6*VoP*tan(phyEP)*cos(phyEP)/(3+sin(phyEP)));

def sig1(P,x,y):
    return (1/P['Nc'])*(P['Vo']+P['Nc']*x+P['Ne']*(y-x)+y)

def get_z_pmc(P,sigma_num,configuration,x,y):
    Vo = P['Vo']
    A, B, C = P['A'], P['B'], P['C']
    if configuration == 'x':
        if sigma_num == 1:
            return (1/A)*(1-x*B-y*C)
        elif signma_num == 2:
            return (1/B)*(1-x*A-y*C)
        elif signma_num == 3:
            return (1/C)*(1-x*A-y*B)
    elif configuration == 'y':
        if sigma_num == 1:
            return (1/A)*(1-y*B-x*C)
        elif signma_num == 2:
            return (1/B)*(1-y*A-x*C)
        elif signma_num == 3:
            return (1/C)*(1-y*A-x*B)
        
        
conf_cycle = [1,2,3,3,2,1]
form_cycle = ['x','y']

def get_plane_n(P,num):
    return get_plane(P,conf_cycle[num%6],form_cycle[(num//3)%2])

def get_plane(P,sigma,form):
    P_coef = [P['A'],P['B'],P['C']]
    normal = np.zeros(3)
    normal[2] = P_coef[sigma-1]
    del P_coef[sigma-1]
    if form == 'y':
        P_coef = P_coef[::-1]
        # P_coef = [P_coef[1],P_coef[0]]
    normal[:2] = P_coef
    #origin = 1/(normal*3)
    
    return normal #origin