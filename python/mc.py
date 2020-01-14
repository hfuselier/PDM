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
            
def qCfit(P,x): # Compute q for plane fitting of compression data (p-q coordinate system)
    phy = P.phy
    c = P.c
    return (6*sin(phy)/(3-sin(phy)))*x+(6*c*cos(phy)/(3-sin(phy)));

def qEfit(P,x): # Compute q for plane fitting of extension data (p-q coordinate system)
    phy = P.phy
    c = P.c
    return (6*sin(phy)/(3+sin(phy)))*x+(6*c*cos(phy)/(3+sin(phy)));

conf_cycle = [1,2,2,3,3,1]
form_cycle = ['x','x','y','y','x','y']
conf_cycle = [1,3,3,2,2,1]
form_cycle = ['y','x','y','y','x','x']
#conf_cycle = [1,2,3,3,2,1]
#form_cycle = ['x','y','x','y','x','y']

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

def p_plane_intersection_6(P,plane_dist):
    p_plane_n = np.ones(3)/np.sqrt(3)
    p_plane_o = p_plane_n*plane_dist

    eq = np.zeros((6,3,3))
    b = np.array([[1,1,np.dot(p_plane_n,p_plane_o)]])
    for i in range(6):
        eq[i] = np.array([get_plane_normal_6_cycle(P,i),get_plane_normal_6_cycle(P,i-1),p_plane_n])
    
    pts = np.linalg.solve(eq,b).transpose()
    return pts

# From the Plane class we get all necessary data for P1 or P2
class Plane_MC:
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
        
        # Fitting parameters and coefficients for Mohr-Coulomb
        self.Co = 43.092
        #self.K = (get_C(data)[:,0]-self.Co)/get_C(data)[:,2]
        #K = []
        #err = []
        #bset = []
        #for i in range(self.K.size):
        #        if np.isnan(self.K[i]) == False and np.isinf(self.K[i]) == False :
        #            K.append(self.K[i])
        #for j in np.linspace(2.0, 5.0, num=50):
        #    bset.append(j)
        #    diff_calc = get_C(data)[:,2]*j+self.Co-get_C(data)[:,0]
        #    diff = np.mean(diff_calc)
        #    err.append(diff)
        #err = np.array(err)
        #Kp = bset[np.argmin(abs(err))]
        self.Kp = 2.119
        self.phi = arcsin((self.Kp-1)/(self.Kp+1))
        self.Vo = self.Co/(self.Kp-1) 
        self.c = self.Co*(1-sin(self.phi))/(2*cos(self.phi))
        self.sol = np.array([self.phi*180/pi,self.c,self.Co,self.Kp,self.Vo])
        self.A = 1/self.Co
        self.B = 0
        self.C = -self.Kp/self.Co
        
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