from idlsave import  idlsave 
import numpy as np
import healpy
import matplotlib.pypot as plt
import scipy
betw = lambda x,x1,x2 : (x>x1)&(x<x2)

exec idlsave.restore('gaia_scan_512_gal.psav')

#locs=np.r_[np.arange(48)/24.,np.arange(2,365),np.arange(365,5*365,5)]
locs=np.r_[1,10,100,500]

def doit(x):
    return scipy.histogram(x,locs)[0]

superres=[doit(_) for _ in obstimes]
arr=np.array(superrres)
arr=np.array(superres)

for i in range(nx):
    plt.clf();
    healpy.mollview((arr[:,:i]).sum(axis=1),fig=1,nest=True,hold=True,
    	title='%.03f days'%locs[i],min=0,max=max(((arr[:,:i]).sum(axis=1)).max(),1),xsize=2000);
    plt.savefig('xx_%04d.png'%i,dpi=600)
    
