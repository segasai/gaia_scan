

from idlplotInd import plot, oplot, plothist, ploterror ,tvhist2d
from idlplot import tvaxis, contour
from idlsave import  idlsave 
from readcol import readcol
from mpfit import mpfit
from mpfitexpr import mpfitexpr
import numpy,numpy.random,atpy,sqlutil,pyfits,healpy,math,scipy
import scipy.stats,scipy.special,numexpr,euler,h5py
import matplotlib.pyplot as plt
import tables
import numpy as np
import glob
import sphdist
betw = lambda x,x1,x2 : (x>x1)&(x<x2)
def doit2(x):
    d=np.r_[np.diff(np.sort(x)),np.inf]
    x1=np.minimum(d,np.r_[np.inf,d[:-1]])
    return scipy.histogram(x1,locs)[0]

def doit(x):
    return scipy.histogram(np.diff(np.sort(x)),locs)[0]

exec idlsave.restore('gaia_scan_512.psav',printVars=True)


from idlplotInd import plot, oplot, plothist, ploterror ,tvhist2d
from idlplot import tvaxis, contour
from idlsave import  idlsave 
from readcol import readcol
from mpfit import mpfit
from mpfitexpr import mpfitexpr
import numpy,numpy.random,atpy,sqlutil,pyfits,healpy,math,scipy
import scipy.stats,scipy.special,numexpr,euler,h5py
import matplotlib.pyplot as plt
import tables
import numpy as np
import glob
import sphdist
betw = lambda x,x1,x2 : (x>x1)&(x<x2)
exec idlsave.restore('gaia_scan_256.psav',printVars=True)
len(hh[0])
type(hh)
type(hh[0])
type(obstimes)
type(obstimes[0])
sorted(obstimes[0])
type(sorted(obstimes[0]))
np.sort
['2',3]


from idlplotInd import plot, oplot, plothist, ploterror ,tvhist2d
from idlplot import tvaxis, contour
from idlsave import  idlsave 
from readcol import readcol
from mpfit import mpfit
from mpfitexpr import mpfitexpr
import numpy,numpy.random,atpy,sqlutil,pyfits,healpy,math,scipy
import scipy.stats,scipy.special,numexpr,euler,h5py
import matplotlib.pyplot as plt
import tables
import numpy as np
import glob
import sphdist
betw = lambda x,x1,x2 : (x>x1)&(x<x2)
import gaiarot
obstimes[0]
len(obstimes[0])
dtype(obstimes[0])
type(obstimes[0])
obstimes1=[np.array(_) for _ in obstimes]
del obstimes
obstimes=obstimes1
idlsave.save('gaia_scan_512.psav','hh,obstimes',hh,obstimes)
superres1=x[doit2(_) for _ in obstimes]


from idlplotInd import plot, oplot, plothist, ploterror ,tvhist2d
from idlplot import tvaxis, contour
from idlsave import  idlsave 
from readcol import readcol
from mpfit import mpfit
from mpfitexpr import mpfitexpr
import numpy,numpy.random,atpy,sqlutil,pyfits,healpy,math,scipy
import scipy.stats,scipy.special,numexpr,euler,h5py
import matplotlib.pyplot as plt
import tables
import numpy as np
import glob
import sphdist
betw = lambda x,x1,x2 : (x>x1)&(x<x2)


from idlplotInd import plot, oplot, plothist, ploterror ,tvhist2d
from idlplot import tvaxis, contour
from idlsave import  idlsave 
from readcol import readcol
from mpfit import mpfit
from mpfitexpr import mpfitexpr
import numpy,numpy.random,atpy,sqlutil,pyfits,healpy,math,scipy
import scipy.stats,scipy.special,numexpr,euler,h5py
import matplotlib.pyplot as plt
import tables
import numpy as np
import glob
import sphdist
betw = lambda x,x1,x2 : (x>x1)&(x<x2)
exec idlsave.restore('gaia_scan_512_gal.psav',printVars=True)

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
    
