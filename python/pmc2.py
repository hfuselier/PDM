import numpy as np
from numpy import sin, tan ,cos, pi, sqrt, arcsin
from tqdm.auto import tqdm

def get_C(data): # Gives all points from conventional compression (CTC)
    return data[data[:,5] == 0]
def get_E(data): # Gives all points from conventional extension (CTE)
    return data[data[:,5] == 60]
def get_o(data): # Gives all points from True triaxial testing (TT)
    return data[np.logical_and(data[:,5] != 60 , data[:,5] != 0)]

def get_p(data): # Gives mean stress for all the points 
    return data[:,3]
def get_q(data): # Gives deviatoric stress of all the points
    return data[:,4]
def get_t(data): # Gives Lodge angle of all the points
    return data[:,5]

# Creates fitting Planes P2 and P1
def create_P1_and_P2_initialization(data,bpC,bpE): 
    P1_C = get_C(data)[get_p(get_C(data)) >= bpC]
    P2_C = get_C(data)[get_p(get_C(data)) < bpC]

    P1_E = get_E(data)[get_p(get_E(data)) >= bpE]
    P2_E = get_E(data)[get_p(get_E(data)) < bpE]
    
    # The factor 2 shouldn't be changed. Normally is the same for all the cases
    P1_o = get_o(data)[get_o(data)[:,0] > 2*get_o(data)[:,1]]
    P2_o = get_o(data)[get_o(data)[:,0] >= 2*get_o(data)[:,1]]

    P1 = Plane((P1_C,P1_E,P1_o))
    P2 = Plane((P2_C,P2_E,P2_o))
    
    return P1, P2

def create_P1_and_P2_iteration(P1_C,P2_C,P1_E,P2_E,P1_o,P2_o): 
    P1 = Plane((P1_C,P1_E,P1_o))
    P2 = Plane((P2_C,P2_E,P2_o))
    
    return P1, P2
            
def qCfit(P,x): # Compute q for plane fitting of compression data (p-q coordinate system)
    phyCP = P.phyC
    VoP = P.Vo
    return (6*sin(phyCP)/(3-sin(phyCP)))*x+(6*VoP*tan(phyCP)*cos(phyCP)/(3-sin(phyCP)));

def qEfit(P,x): # Compute q for plane fitting of extension data (p-q coordinate system)
    phyEP = P.phyE
    VoP = P.Vo
    return -(6*sin(phyEP)/(3+sin(phyEP)))*x-(6*VoP*tan(phyEP)*cos(phyEP)/(3+sin(phyEP)));

def sig1(P,x,y): # Compute sig for plane fitting (sig2-sig1 coordinate system)
    return (1/P.Nc)*(P.Vo+P.Nc*x+P.Ne*(y-x)+y)

conf_cycle = [1,2,2,3,3,1]
form_cycle = ['x','x','y','y','x','y']
conf_cycle = [1,3,3,2,2,1]
form_cycle = ['y','x','y','y','x','x']
#conf_cycle = [1,2,3,3,2,1]
#form_cycle = ['x','y','x','y','x','y']

def get_plane_normal_12_cycle(P1,P2,num,offset):
    num = (num + offset)%12
    p_cycle = [P1,P2,P2,P1] 
    P = p_cycle[num%4]
    return get_plane_normal_6_cycle(P,num//2)

def get_plane_normal_6_cycle(P,num):
    return get_plane_normal(P,conf_cycle[num%6],form_cycle[num%6])

def get_plane_normal(P,sigma,form):
    
    if sigma == 1 and form == 'x' : 
        normal = [P.A,P.B,P.C]
    if sigma == 1 and form == 'y' : 
        normal = [P.A,P.C,P.B]
    if sigma == 2 and form == 'x' : 
        normal = [P.B,P.A,P.C]
    if sigma == 2 and form == 'y' : 
        normal = [P.C,P.A,P.B]
    if sigma == 3 and form == 'x' : 
        normal = [P.B,P.C,P.A]
    if sigma == 3 and form == 'y' : 
        normal = [P.C,P.B,P.A]
    
    return np.array(normal) #origin

def p_plane_intersection_12(P1,P2,plane_dist,offset):
    p_plane_n = np.ones(3)/np.sqrt(3)
    p_plane_o = p_plane_n*plane_dist
    
    eq = np.zeros((12,3,3))
    b = np.array([[1,1,np.dot(p_plane_n,p_plane_o)]])
    
    for i in range(12):
        eq[i] = np.array([get_plane_normal_12_cycle(P1,P2,i,offset),get_plane_normal_12_cycle(P1,P2,i-1,offset),p_plane_n])
        
    pts = np.linalg.solve(eq,b).transpose()
    return pts

def p_plane_intersection_6(P,plane_dist):
    p_plane_n = np.ones(3)/np.sqrt(3)
    p_plane_o = p_plane_n*plane_dist

    eq = np.zeros((6,3,3))
    b = np.array([[1,1,np.dot(p_plane_n,p_plane_o)]])
    for i in range(6):
        eq[i] = np.array([get_plane_normal_6_cycle(P,i),get_plane_normal_6_cycle(P,i-1),p_plane_n])
    
    pts = np.linalg.solve(eq,b).transpose()
    return pts

def p_planes(pconst,P1,P2,p_trans):
    
    if pconst<= max(p_trans) and pconst>= min(p_trans):
        pts = p_plane_intersection_12(P1,P2,pconst*sqrt(3))
    elif pconst> max(p_trans):
        pts = p_plane_intersection_6(P1,pconst*sqrt(3))
    elif pconst < min(p_trans):
        pts = p_plane_intersection_6(P2,pconst*sqrt(3))
    
    return pts

# From the Plane class we get all necessary data for P1 or P2
class Plane:
    def __init__(self,data):
        if type(data)==tuple:
            self.data = np.concatenate(data,axis=0)
        else:
            self.data = data
        
        self.sig123 = self.data[:,:3].transpose()
        self.pts = np.zeros(self.sig123.shape)
        
        self.pts[0,:] = self.sig123[1,:]
        self.pts[1,:] = self.sig123[2,:]
        self.pts[2,:] = self.sig123[0,:] 
        
        # Get mean stress, deviatoric stress and Lodge angle (in degree)
        self.p = get_p(self.data)
        self.q = get_q(self.data)
        self.t = get_t(self.data)
        self.t = self.t * pi/180
        
        # Computation of the fitting solution
        a = self.q*np.sin(self.t)

        A = np.ndarray((self.t.shape[0],3))
        A[:,0] = self.p
        A[:,1] = a
        A[:,2] = 1

        b = self.q*np.cos(self.t)

        x = np.linalg.pinv(A)@b
        
        # Solutions parameters
        bc = x[2]
        k = x[1]
        self.k = k
        self.Vo = bc/x[0]
        be = 2*bc/(1-sqrt(3)*k)
        self.bc = bc
        self.be = be
        self.phyC = arcsin(3*bc/(6*self.Vo+bc))
        self.phyE = arcsin(3*be/(6*self.Vo-be))
        #if np.isin(60*pi/180, self.t)== True :
        #    self.phyE = arcsin(3*be/(6*self.Vo-be))
        #elif np.isin(60*pi/180, self.t)== False:
        #    self.phyE = self.phyC
        self.sol = np.array([bc, be, k, self.Vo, self.phyC*180/pi, self.phyE*180/pi])
        
        # Computation of Planes Coefficients
        self.cc = (self.bc*(3-sin(self.phyC)))/(6*cos(self.phyC))
        self.ce = (self.be*(3-sin(self.phyE)))/(6*cos(self.phyE))
        
        self.Mc = (1+sin(self.phyC))/(1-sin(self.phyC))
        self.Me = (1+sin(self.phyE))/(1-sin(self.phyE))
        self.Cc = (2*self.cc*cos(self.phyC))/(1-sin(self.phyC))
        self.Ce = (2*self.ce*cos(self.phyE))/(1-sin(self.phyE))
        
        self.mc = (6*sin(self.phyC)/(3-sin(self.phyC)))
        self.me = (6*sin(self.phyE)/(3+sin(self.phyE)))
        
        self.Nc = (1-sin(self.phyC))/(2*sin(self.phyC))
        self.Ne = (1-sin(self.phyE))/(2*sin(self.phyE))
        
        self.A = ((1-sin(self.phyC))/(2*sin(self.phyC)))/self.Vo
        self.B = ((sin(self.phyC)-sin(self.phyE))/(2*sin(self.phyC)*sin(self.phyE)))/self.Vo
        self.C = -((1+sin(self.phyE))/(2*sin(self.phyE)))/self.Vo
              
        
        
    def get_pC(self):
        return get_p(get_C(self.data))
    
    def get_qC(self):
        return get_q(get_C(self.data))
    
    def get_pE(self):
        return get_p(get_E(self.data))
    
    def get_qE(self):
        return get_q(get_E(self.data))
    
    def get_po(self):
        return get_p(get_o(self.data))
    
    def get_qo(self):
        return get_q(get_o(self.data))