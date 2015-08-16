from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# Top half C-mesh for Joukowski airfoil
# Based on conformal mappings

#-----------------------------------
def Distance(i, dx_te, ratio):
    return dx_te*(ratio**i - 1) / (ratio-1)

def GradDist(i, dx_te, ratio):
    return dx_te*((i-1)*ratio**i - i*ratio**(i-1) + 1) / (ratio-1)**2


def FindStretching(n, h_min, Hc):
    # Find the ratio of successive cell sizes to get a total length of Hc from
    # n cells, starting with the first cell size = h_min.
    
    # This guess is exact for very large size ratios
    guess = (Hc / h_min)**(1/(n-1))
    finished = False
    while not finished:
        func = Distance(n, h_min, guess) - Hc
        grad = GradDist(n, h_min, guess)
        delta = - func / grad
        guess += delta
        #print guess, ' ', delta
        if (abs(delta) < 1.e-12):
            finished = True

    return guess;

#===============================================================================
def spaceq(re, Q):
    nsub = Q
    nre = len(re) - 1
    nr  = nsub*nre
    r = np.zeros(nr+1)
    for i in xrange(nre):
        for j in xrange(nsub):
            f = j/nsub
            r[i*nsub+j] = re[i]*(1.0-f) + re[i+1]*f
    r[nr] = re[nre]
    
    return r

#-----------------------------------
def Bezier(nn, smax=1, ds0=-0.2, ds1=-0.2):

    s0 = np.linspace(0,smax,nn+1)
    
    #Use a Bezier curve to cluster at LE and TE: ds = -1 gives a linear distribution. Clustering is added as ds->0 from -1
    #ds0 = -0.2
    #ds1 = -0.2
    P0 = 1
    P1 = (3 + ds1)/3
    P2 = -(ds0/3)
    s1 = P0*(1 - s0)**3 + P1*3*s0*(1 - s0)**2 + P2*3*s0**2*(1 - s0)
    return s1

#-----------------------------------
def meshplot(X, Y, edgecolor='k'):
    """Plot a mapped Cartesian grid."""
    plt.clf()
    plt.axis('equal')
    plt.pcolor(X, Y, 0*X, edgecolor=edgecolor, cmap='Greens')
    plt.show()
    plt.draw()

def joukowski_conformal(S, T, joux=0.1):
    """Conformal mapping (S,T) -> (X,Y) for Joukowski mesh.

    joux:   Joukowski x-shift
    """

    # Special case for (0,0)
    zero_ix = np.logical_and(S == 0.0, T == 0.0).ravel().nonzero()

    # Map to circle
    z = S + 1j * T
    w = z**2
    w.ravel()[zero_ix] = 1.0
    A = (w - 1) / w
    B = np.sqrt(A)
    c = (1 + B) / (1 - B)
    c.ravel()[zero_ix] = -1.0

    # Map to Joukowski
    cs = c - complex(joux, 0.0)
    l = 1.0 - joux
    g = cs + l**2 / cs

    # Scale airfoil to x = [0,1]
    endpnts = np.array([-1, 1])
    left_right = endpnts - joux + l**2 / (endpnts - joux)
    g = (g - left_right[0]) / np.diff(left_right)

    # Return physical and parameter meshes
    return np.real(g), np.imag(g)

def joukowski_inverse(X, Y, joux=0.1):
    """Inverse of joukowski_conformal (X,Y) -> (S,T).

    joux:   Joukowski x-shift
    """
    # Maple auto-generated
    zz = X + 1j * Y
    t1 = (joux ** 2)
    t2 = 2 * t1
    t3 = t1 * zz
    t4 = 3 * t3
    t5 = zz ** 2
    t7 = 2 * joux * t5
    t8 = t1 ** 2
    t15 = t5 ** 2
    t17 = t5 * zz
    t22 = np.sqrt((t8 - 3 * zz * t8 - t3 + 3 * t8 * t5 + 3 * t1 * t5 + t1 * t15 - t8 * t17 - 3 * t1 * t17))
    t23 = 0.2e1 * t22
    t31 = 1 / (-3 * t1 + 4 * t3 + 1 + 4 * zz * joux - 2 * joux)
    t33 = np.sqrt(t31 * (-t2 + t4 + t7 + zz + t23))
    t36 = np.sqrt(-t31 * (t2 - t4 - t7 - zz + t23))
    out = np.array([t33, -t33, t36, -t36])

    # Test closest match between positive roots
    map1 = joukowski_conformal(np.real(out[0]), np.imag(out[0]), joux)
    map2 = joukowski_conformal(np.real(out[2]), np.imag(out[2]), joux)

    d1 = (map1[0] - X)**2 + (map1[1] - Y)**2
    d2 = (map2[0] - X)**2 + (map2[1] - Y)**2

    out1 = out[0]
    ix = (d2 < d1).ravel().nonzero()
    out1.ravel()[ix] = out[2].ravel()[ix]
    
    return np.real(out1), np.imag(out1)

def joukowski_parameter(ref, Q, reynolds, growth=1.3, R=100, joux=0.1):
    """Make parameter space mesh (S,T) for Joukowski mapping.

    nchord: Number of streamwise points along the chord (cos-distributed)
    growth: Element size growth ratio in wake + normal directions
    R:      Farfield distance (at least)
    """
    
    nchord=8*2**ref           # number of elements along one side of the airfoil geometry
    nxwake=8*2**ref           # x-wake on centerline
    nnormal=16*2**ref         # points normal to airfoil surface
    
    # Trailing edge spacing
    if (reynolds > 5e5):
        # Turbulent. 
        ds0 = -0.1
    else:
        # Laminar.  
        ds0 = -0.0
        ds1 = -0.05

    # Chord distribution
    #phi = np.linspace(np.pi, 0.0, nchord+1)
    #sAf = (np.cos(phi) + 1) / 2
    sAf = 1-Bezier(nchord,ds0=ds0,ds1=ds1)

    sAf_half = 1-sAf[:nchord/2-1:-1]
    
    ds = sAf_half[-1] - sAf_half[-2]
    
    sx = np.zeros(nchord+nxwake+1)
    sx[0:nchord+1] = sAf
    sx[nchord:nchord+len(sAf_half)] = 1+sAf_half
    nWake = nxwake+1-len(sAf_half)
    ratio = FindStretching(nWake, ds, np.sqrt((R + 1.5)/1.05)-1.5)
    for i in xrange(1,nWake+1):
        sx[nchord+len(sAf_half)+i-1] = 1.5 + Distance(i, ds, ratio)

    sx = spaceq(sx, Q)

    # Wake distribution
    #sx = sAf
    #sx = np.append(sx, 1.0 + sAf_half[1:])
    #while sx[-1]**2 < (R + 1)/1.5:
    #    sx = np.append(sx, sx[-1] + growth * (sx[-1] - sx[-2]))
    
    
    sy = np.zeros(nnormal+1)
    sy[0:len(sAf_half)] = sAf_half
    nNormal = nnormal+1-len(sAf_half)
    ratio = FindStretching(nNormal, ds, np.sqrt(R/1.05)-0.5)
    for i in xrange(nNormal+1):
        sy[len(sAf_half)+i-1] = 0.5 + Distance(i, ds, ratio)

    sy = spaceq(sy, Q)
    
    # Normal distribution
    #sy = sAf_half.copy()
    #growth_normal = 1.0 + (growth - 1.0) / 2.0  # empirical
    #while sy[-1]**2 < R/1.5:
    #    sy = np.append(sy, sy[-1] + growth_normal * (sy[-1] - sy[-2]))


    lx0 = sx / sx.max()
    ly0 = sy / sy.max()
    lx = -2 * lx0**3 + 3 * lx0**2
    ly = -2 * ly0**3 + 3 * ly0**2

    # Bottom and left
    bottom = [sx, 0*sx]
    left = [0*sy, sy]

    # Find parameters for straight vertical outflow boundary
    xright = joukowski_conformal(sx[-1], 0.0, joux)[0]
    yright = ly0 * -joukowski_conformal(0.0*sy[-1], sy[-1], joux)[0]
    right = joukowski_inverse(xright + 0*sy, yright, joux)
    right_eps = joukowski_inverse(xright, yright[-1] - 1e-5, joux)

    # Top boundary

    # Straight
    #top = [np.linspace(left[0][-1], right[0][-1], sx.shape[0]),
    #       np.linspace(left[1][-1], right[1][-1], sx.shape[0])]
    # Hermite
    lxtop = lx.copy()
    # Some smoothing, but smoothing removes grid nesting
    #for i in range(1000): # empirical
    #    lxtop[1:-1] = (lxtop[0:-2] + lxtop[2:]) / 2.0
    lxtop = np.linspace(0, lxtop[-1], lxtop.shape[0]) # straight
    lxtop1 = -2 * lxtop**3 + 3 * lxtop**2
    # Make orthogonal to right boundary
    #slope = (right[0][-1] - right[0][-2]) / (right[1][-1] - right[1][-2])
    slope = (right[0][-1] - right_eps[0]) / (right[1][-1] - right_eps[1])
    newslope = slope * right[0][-1]
    lxtop2 = lxtop**3 - lxtop**2
    top = [lxtop * right[0][-1], \
           left[1][-1]*(1-lxtop1) + right[1][-1]*lxtop1 - newslope*lxtop2]
    
    if False:  #debugging
        plt.clf()
        plt.plot(bottom[0], bottom[1], '.-')
        plt.plot(top[0], top[1], '.-')
        plt.plot(left[0], left[1], '.-')
        plt.plot(right[0], right[1], '.-')
        plt.axis('equal')
        plt.draw()

    # TFI mapping
    X1 = np.outer(1-lx, left[0]) + np.outer(lx, right[0])
    Y1 = np.outer(1-lx, left[1]) + np.outer(lx, right[1])

    X2 = np.outer(bottom[0], 1-ly) + np.outer(top[0], ly)
    Y2 = np.outer(bottom[1], 1-ly) + np.outer(top[1], ly)

    X12 = np.outer(1-lx, 1-ly) * X1[0,0] + np.outer(lx, 1-ly) * X1[-1,0] + \
          np.outer(1-lx, ly) * X1[0,-1] + np.outer(lx, ly) * X1[-1,-1]
    Y12 = np.outer(1-lx, 1-ly) * Y1[0,0] + np.outer(lx, 1-ly) * Y1[-1,0] + \
          np.outer(1-lx, ly) * Y1[0,-1] + np.outer(lx, ly) * Y1[-1,-1]

    X = X1 + X2 - X12
    Y = Y1 + Y2 - Y12

    return X, Y

def make_joukowski_classic(ref, Q, reynolds=1.e6):
    S, T = joukowski_parameter(ref, Q, reynolds)
    X, Y = joukowski_conformal(S, T)
    
    X = np.concatenate(( np.flipud(np.delete(X, 0, axis=0)), X), axis=0)
    Y = np.concatenate((-np.flipud(np.delete(Y, 0, axis=0)), Y), axis=0)

    return X, Y

if __name__ == "__main__":
    X, Y = make_joukowski_classic(0, 1, 1e6)
    meshplot(X, Y)
