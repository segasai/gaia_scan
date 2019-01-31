#
# Program to generate the gaia scanning law 
#
# Sergey Koposov 2013
# koposov@ast.cam.ac.uk
#
# Non-standard dependencies: sphere_rotate.py and idlsave.py available from
# http://astrolibpy.googlecode.com/git/utils/idlsave.py
# http://astrolibpy.googlecode.com/git/my_utils/sphere_rotate.py
# Standard dependencies numpy, scipy, healpy
#

import scipy, numpy as np,scipy.spatial, multiprocessing as mp, healpy, euler
import sphere_rotate
import pickle 

def pix2radec(nside,pix,nest=True):
    __, _ = healpy.pix2ang(nside, pix, nest=nest)
    ra, dec = np.rad2deg(_), np.rad2deg(np.pi/2 - __)
    return ra,dec

def getmatrices(vec, ang):
    # get the rotation matrix given the axis and the rotation angle
    ux, uy, uz = vec[0], vec[1], vec[2]
    ca = np.cos(ang)
    sa = np.sin(ang)
    Ca = 1 - np.cos(ang)
    a  = {}
    a[0,0] = ca + ux**2 * Ca
    a[1,0] = ux * uy * ca- uz * sa
    a[2,0] = ux * uz * Ca + uy * sa

    a[0,1] = Ca * ux * uy + uz * sa
    a[1,1] = ca + uy**2 * Ca
    a[2,1] = Ca * uy*uz - ux * sa

    a[0,2] = ux * uz*Ca - uy * sa
    a[1,2] = uz * uy*Ca + ux * sa
    a[2,2] = ca + uz**2 * Ca
    return a
        
def rotateVec(mats, vec):
    # perform the rotation
    a = mats
    vec0 = a[0,0]*vec[0]+a[1,0]*vec[1]+a[2,0]*vec[2]
    vec1 = a[0,1]*vec[0]+a[1,1]*vec[1]+a[2,1]*vec[2]
    vec2 = a[0,2]*vec[0]+a[1,2]*vec[1]+a[2,2]*vec[2]
    return np.array([vec0,vec1,vec2])

def xyztoradec(vec):
    # xyz - > (ra,dec)
    x,y,z=vec[0],vec[1],vec[2]
    ra = np.rad2deg(np.arctan2(y,x))
    dec = np.rad2deg(np.arctan2(z,np.sqrt(x**2+y**2)))
    return ra,dec
    
def scanninglaw(times1): # times in days
    # return the XYZs of the scanning law for two gaia fields at the specified times
    ang1 = np.deg2rad(106.5)
    vec1 = [1,0,0]
    vec2 = [np.cos(ang1), np.sin(ang1), 0]
    freq1  = 1 * 60 * 24 # 1 deg /per min  in deg/day
    phase1 = np.deg2rad(freq1 * times1 )
    vecrot1 = [0, 0, 1]  # z axis
    mats = getmatrices(vecrot1, phase1)
    vec1_1 = rotateVec(mats, vec1)
    vec2_1 = rotateVec(mats, vec2)
    # rotation of the satellite itself
    
    vecrot2 = [0,1,0] 
    ang2 = np.deg2rad(45) # solar aspect angle 
    # 45 deg inclination to sun

    freq2 = 360/63. # 360 deg per 63 days
    mats = getmatrices(vecrot2,ang2)
    vec1_2 = rotateVec(mats, vec1_1)
    vec2_2 = rotateVec(mats, vec2_1)


    phase3 = np.deg2rad(freq2*times1)        
    vecrot3=[0,0,1] # sun in the z axis
    mats = getmatrices(vecrot3, phase3)
    vec1_3 = rotateVec(mats, vec1_2)
    vec2_3 = rotateVec(mats, vec2_2)
    # precession around the sun
        
    mats = getmatrices([0, 1, 0],np.deg2rad(90))
    vec1_4 = rotateVec(mats, vec1_3)
    vec2_4 = rotateVec(mats, vec2_3)
    # placing the sun in the plane

    freq3 = 360./365.25 # 360 deg per 365.25 days
    phase5 = np.deg2rad(freq3 * times1)    
    mats = getmatrices([0, 0, 1], phase5)
    vec1_5 = rotateVec(mats, vec1_4)
    vec2_5 = rotateVec(mats, vec2_4)
    # rotation of the sun
    return vec1_5, vec2_5

class si:
    tree1 = None
    tree2 = None
    times = None

def numcovers(xyz):
    # get the number and times of observations of the given point
    aperture = 0.33 # deg
    dist = np.sin(np.deg2rad(aperture))
    ids1 = si.tree1.query_ball_point(xyz, dist, 2)
    ids2 = si.tree2.query_ball_point(xyz, dist, 2)
    num = 0
    obstimes = []
    for i in [ids1, ids2]:
        isort = np.sort(i).astype(np.int)
        mask = np.r_[[-2],np.array(isort)]
        #mask = np.zeros(len(si.times)+1)
        #mask[i+1] = 1
        curt = np.diff(mask)>1
        obstimes.append(si.times[isort[curt]])
        num += curt.sum()
    return num,np.concatenate(obstimes)
    
def coveragemap(nside=64, nyears=5,galactic=False, nthreads=24):
    # returns the number of visits of the healpix map of given resolution 
    # as well as the times of observations for a given healpixel
    # The healpixes are in the NESTED scheme
    tottime = nyears * 365.25 # 5yrs
    step = 0.00005 # days
    times = np.arange(tottime/step) * step
    ntimes = len(times)
    nchunks = 10
    pool = mp.Pool(nthreads)

    #vecs =  scanninglaw(times)
    xtimes = [times[i*(ntimes//nthreads+1):(i+1)*(ntimes//nthreads+1)] for i in range(nthreads)]
    xvecs =  pool.map(scanninglaw,xtimes)
    pool.close()
    pool.join()

    _tmp1 = np.hstack([_[0] for _ in xvecs])
    _tmp2 = np.hstack([_[1] for _ in xvecs])
    vecs = _tmp1,_tmp2    
    cosd = lambda x : np.cos(np.deg2rad(x))
    sind = lambda x : np.sin(np.deg2rad(x))
    getxyz = lambda r, d: [cosd(r)*cosd(d), sind(r)*cosd(d), sind(d)]
    
    tree1 = scipy.spatial.cKDTree(vecs[0].T)
    tree2 = scipy.spatial.cKDTree(vecs[1].T)

    si.tree1 = tree1
    si.tree2 = tree2
    si.times=times        

    ipixes = np.arange(12 * nside**2)
    if galactic:
        ls,bs = pix2radec(nside,ipixes)
        ras,decs = euler.euler(ls,bs,2)
    else:
        ras, decs = pix2radec(nside, ipixes)
    ra0 = 18 * 15
    dec0 = 66 + 33/60. + 38.55/3600. # pole of the ecliptic
    eclipra, eclipdec = sphere_rotate.sphere_rotate(ras,decs,ra0,dec0,0)
    eclipra = (eclipra + 360) % 360     
    pool = mp.Pool(nthreads)
     
    res = pool.map(numcovers, np.array(getxyz(eclipra, eclipdec)).T, nchunks)
    #res = map(numcovers,np.array(getxyz(eclipra,eclipdec)).T)
    pool.close()
    pool.join()
    poss = np.array([_[0] for _ in res])
    obstimes = [_[1] for _ in res]
    return poss, obstimes


def doit():
    # main loop
    hh, obstimes = coveragemap(128)
    #construct the number of coverings map and get the observation times
    with open('gaia_scan.pkl', 'wb') as fp:
        pickle.dump((hh,obstimes),fp)

if __name__ == '__main__':
    doit()
