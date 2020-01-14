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

# From the Plane class we get all necessary data for P1 or P2
class Plane_HB:
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
       
        # Fitting parameters and coefficients for Hoek-Brown
        #self.Co = data[0,0] 
        #self.mc = (self.Co/get_C(data)[:,2])*(np.square((get_C(data)[:,0]-get_C(data)[:,2])/self.Co)-1)
        
        #mtC = []
        #bset = []
        #err = []
        #for i in range(self.mc.size):
        #        if np.isnan(self.mc[i]) == False and np.isinf(self.mc[i]) == False :
        #            mtC.append(self.mc[i])
        #for j in np.linspace(int(np.min(mtC)), int(np.max(mtC)), num=100):
        #    bset.append(j)
        #    diff_calc = get_C(data)[:,2]+self.Co*np.sqrt((j/self.Co)*get_C(data)[:,2]+1)-get_C(data)[:,0]
        #    diff = np.mean(diff_calc)
        #    err.append(diff)
        #err = np.array(err)
        #m = bset[np.argmin(abs(err))]
        
        ## All points
        #self.Co = 42.35 
        #self.m = 3.1025
        
        ## Six points
        self.Co = 41.40
        self.m = 3.6479
        self.Vo = self.Co/self.m
        
        self.sol = np.array([self.Co,self.m,self.Vo])
        
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