#cython: language_level=3,boundscheck=False, wraparound=False,nonecheck=False
#import cython
import numpy as np
cimport numpy as np

ctypedef np.float64_t Dtype
#ctypedef double Dtype

#-------------------------------------------------------------
cdef extern from "fcic.h":
    void forward_cic(Dtype *cic_in,Dtype *x1_in,Dtype *x2_in,Dtype bsx,Dtype bsy,int nx1,int nx2,int npp,Dtype *cic_out)

#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.nonecheck(False)

def call_forward_cic(int nx1,int nx2,Dtype boxsize,
                     np.ndarray[Dtype, ndim=1, mode="c"] yif1,
                     np.ndarray[Dtype, ndim=1, mode="c"] yif2):

    cdef int ylen = len(yif1)
    cdef np.ndarray img_in = np.array(np.ones(ylen),dtype=np.double)
    cdef np.ndarray img_out = np.zeros((nx1,nx2),dtype=np.double)
    forward_cic(<Dtype *>img_in.data,<Dtype *>yif1.data,<Dtype *>yif2.data,boxsize,boxsize,nx1,nx2,ylen,<Dtype *>img_out.data)
    return img_out.T
#-------------------------------------------------------------

cdef extern from "cv_test.h":
    int sign(Dtype x)
    Dtype deg2rad(Dtype pha)
    void xy_rotate(Dtype *x1_in,Dtype *x2_in,int nx1,int nx2,Dtype xc1,Dtype xc2,Dtype pha,Dtype *x_out,Dtype *y_out)
    void gauss_2d(Dtype *x1,Dtype *x2, int nx1,int nx2,Dtype *par,Dtype *res)
    void tophat_2d(Dtype *x1,Dtype *x2, int nx1,int nx2,Dtype *par,Dtype *res)
    void lq_nie(Dtype *x1,Dtype *x2,int nx1,int nx2,Dtype *lpar,Dtype *alpha1,Dtype *alpha2)
    void find_critical_curve(Dtype *mu,int nx1,int nx2,Dtype* res)

cdef int sign_cy(Dtype x):
    return sign(x)

cdef Dtype deg2rad_cy(Dtype pha):
    return deg2rad(pha)

cdef call_xy_rotate(np.ndarray[Dtype, ndim=2, mode="c"] x1_in,
                    np.ndarray[Dtype, ndim=2, mode="c"] x2_in,
                    Dtype xc1,Dtype xc2,Dtype pha):
    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(x1_in)
    cdef np.ndarray x1_out = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray x2_out = np.zeros((nx1,nx2),dtype=np.double)
    xy_rotate(<Dtype *>x1_in.data,<Dtype *>x2_in.data,nx1,nx2,xc1,xc2,pha,<Dtype *>x1_out.data,<Dtype *>x2_out.data) 
    return x1_out,x2_out
     
#def call_xy_rotate(x, y, xcen, ycen, phi):
#
#    phirad = np.deg2rad(phi)
#    xnew = (x - xcen) * np.cos(phirad) + (y - ycen) * np.sin(phirad)
#    ynew = (y - ycen) * np.cos(phirad) - (x - xcen) * np.sin(phirad)
#    return (xnew,ynew)
##-------------------------------------------------------------
cdef call_gauss_2d(np.ndarray[Dtype, ndim=2, mode="c"] x1,
                  np.ndarray[Dtype, ndim=2, mode="c"] x2,
                  np.ndarray[Dtype, ndim=1, mode="c"] par):
    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(x1)
    cdef np.ndarray res = np.zeros((nx1,nx2),dtype=np.double)
    gauss_2d(<Dtype *>x1.data,<Dtype *>x2.data,nx1,nx2,<Dtype *>par.data,<Dtype *>res.data)
    return res

#def call_gauss_2d(x, y, par):
#    #gpars = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis])
#    #gpars = np.asarray([aps,l_sigs,ylcs,xlcs,qls,phis])
#
#    (xnew,ynew) = call_xy_rotate(x, y, par[0], par[1], par[5])
#    r_ell_sq = ((xnew**2)*par[2] + (ynew**2)/par[2]) / np.abs(par[4])**2
#    return par[3] * np.exp(-0.5*r_ell_sq)
##-------------------------------------------------------------
cdef call_tophat_2d(np.ndarray[Dtype, ndim=2, mode="c"] x1,
                   np.ndarray[Dtype, ndim=2, mode="c"] x2,
                   np.ndarray[Dtype, ndim=1, mode="c"] par):
    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(x1)
    cdef np.ndarray res = np.zeros((nx1,nx2),dtype=np.double)
    tophat_2d(<Dtype *>x1.data,<Dtype *>x2.data,nx1,nx2,<Dtype *>par.data,<Dtype *>res.data)
    return res

#def call_tophat_2d(x, y, par):
#    (xnew,ynew) = call_xy_rotate(x, y, par[0], par[1], par[5])
#    r_ell = np.sqrt(((xnew**2)*par[2] + (ynew**2)/par[2]))
#    res = r_ell*0.0
#    res[r_ell>=par[4]] = -1.0
#    res[r_ell<par[4]] = 10000.0
#    return res
##-------------------------------------------------------------
cdef call_lq_nie(np.ndarray[Dtype, ndim=2, mode="c"] x1,
                np.ndarray[Dtype, ndim=2, mode="c"] x2,
                np.ndarray[Dtype, ndim=1, mode="c"] par):
    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(x1)
    cdef np.ndarray alpha1 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray alpha2 = np.zeros((nx1,nx2),dtype=np.double)
    lq_nie(<Dtype *>x1.data,<Dtype *>x2.data,nx1,nx2,<Dtype *>par.data,<Dtype *>alpha1.data,<Dtype *>alpha2.data)
    return alpha1,alpha2

#def call_lq_nie(x1,x2,lpar):
#
#    xc1 = lpar[0]
#    xc2 = lpar[1]
#    q   = lpar[2]
#    rc  = lpar[3]
#    re  = lpar[4]
#    pha = lpar[5]
#
#    phirad = np.deg2rad(pha)
#    cosa = np.cos(phirad)
#    sina = np.sin(phirad)
#
#    xt1 = (x1-xc1)*cosa+(x2-xc2)*sina
#    xt2 = (x2-xc2)*cosa-(x1-xc1)*sina
#
#    phi = np.sqrt(xt2*xt2+xt1*q*xt1*q+rc*rc)
#    sq = np.sqrt(1.0-q*q)
#    pd1 = phi+rc/q
#    pd2 = phi+rc*q
#    fx1 = sq*xt1/pd1
#    fx2 = sq*xt2/pd2
#    qs = np.sqrt(q)
#
#    a1 = qs/sq*np.arctan(fx1)
#    a2 = qs/sq*np.arctanh(fx2)
#
#    res1 = (a1*cosa-a2*sina)*re
#    res2 = (a2*cosa+a1*sina)*re
#    return res1,res2
##--------------------------------------------------------------------
def call_find_critical_curve(np.ndarray[Dtype, ndim=2, mode="c"] mu):
    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(mu)
    cdef np.ndarray res = np.zeros((nx1,nx2),dtype=np.double)
    find_critical_curve(<Dtype *>mu.data,nx1,nx2,<Dtype *>res.data)
    return res

#def call_find_critical_curve(mu):
#    rows,cols = np.indices(np.shape(mu))
#    cdtn = np.sign(mu)*(np.sign(mu[rows-1,cols])+np.sign(mu[rows,cols-1])+np.sign(mu[(rows+1)%len(rows),cols])+np.sign(mu[rows,(cols+1)%len(cols)]))
#
#    res = mu*0.0
#    res[cdtn<4] = 1
#    res[cdtn>=4] = 0
#
#    return res
#--------------------------------------------------------------------
#def tot_lq(xi1,xi2,lpar,lpars):
#
#    al1,al2 = call_lq_nie(xi1,xi2,lpar)
#    for i in lpars:
#        al1s,al2s = call_lq_nie(xi1,xi2,i)
#        al1 = al1+al1s
#        al2 = al2+al2s
#
#    yi1 = xi1-al1
#    yi2 = xi2-al2
#
#    return yi1,yi2

cdef tot_lq(np.ndarray[Dtype, ndim=2, mode="c"] x1,
            np.ndarray[Dtype, ndim=2, mode="c"] x2,
            np.ndarray[Dtype, ndim=1, mode="c"] lpar,
            list lpars):

    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(x1)
    cdef np.ndarray al1 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray al2 = np.zeros((nx1,nx2),dtype=np.double)
    lq_nie(<Dtype *>x1.data,<Dtype *>x2.data,nx1,nx2,<Dtype *>lpar.data,<Dtype *>al1.data,<Dtype *>al2.data)

    cdef np.ndarray al1s = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray al2s = np.zeros((nx1,nx2),dtype=np.double)

    llen = len(lpars)
    cdef np.ndarray lpars_i = np.zeros((nx1,nx2),dtype=np.double)
    for i in xrange(llen):
        lpars_i = lpars[i]
        lq_nie(<Dtype *>x1.data,<Dtype *>x2.data,nx1,nx2,<Dtype *>lpars_i.data,<Dtype *>al1s.data,<Dtype *>al2s.data)
        al1 = al1+al1s
        al2 = al2+al2s

    cdef np.ndarray y1 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray y2 = np.zeros((nx1,nx2),dtype=np.double)
    y1 = x1-al1
    y2 = x2-al2

    return y1,y2
#@profile
#def refine_critical(lpar,lpars,critical,xi1,xi2,dsx,nfiner=4):
#    dsf = dsx/nfiner/2
#    x1t = []
#    x2t = []
#    for i in xrange(nfiner):
#        for j in xrange(nfiner):
#            x1tmp = xi1[critical>0]+(dsf*(1-nfiner)*0.5)+dsf*i
#            x2tmp = xi2[critical>0]+(dsf*(1-nfiner)*0.5)+dsf*j
#            x1t.append(x1tmp)
#            x2t.append(x2tmp)
#
#    yift1,yift2 = tot_lq(x1t,x2t,lpar,lpars)
#
#    return yift1,yift2
cdef refine_critical(np.ndarray[Dtype, ndim=2, mode="c"] xi1,
                     np.ndarray[Dtype, ndim=2, mode="c"] xi2,
                     np.ndarray[Dtype, ndim=1, mode="c"] lpar,
                     list lpars,
                     np.ndarray[Dtype, ndim=2, mode="c"] critical):

    cdef nfiner = 4
    cdef double dsx = xi1[1,1]-xi1[0,0]
    cdef dsf = dsx/nfiner/4
    cdef int clen
    clen = len(critical[critical>0])

    cdef np.ndarray x1 = np.zeros((clen),dtype=np.double)
    cdef np.ndarray x2 = np.zeros((clen),dtype=np.double)
    cdef np.ndarray xt1 = np.zeros((clen,nfiner*nfiner),dtype=np.double)
    cdef np.ndarray xt2 = np.zeros((clen,nfiner*nfiner),dtype=np.double)

    x1 = xi1[critical>0]
    x2 = xi2[critical>0]

    cdef int i,j
    for i in xrange(nfiner):
        for j in xrange(nfiner):
            xt1[:,i*nfiner+j] = x1+(dsf*(1-nfiner)*0.5)+dsf*i
            xt2[:,i*nfiner+j] = x2+(dsf*(1-nfiner)*0.5)+dsf*j

    cdef np.ndarray yift1 = np.zeros((clen,nfiner*nfiner),dtype=np.double)
    cdef np.ndarray yift2 = np.zeros((clen,nfiner*nfiner),dtype=np.double)
    yift1,yift2 = tot_lq(xt1,xt2,lpar,lpars)

    return yift1.flatten(),yift2.flatten()

def lensed_images(np.ndarray[Dtype, ndim=2, mode="c"] xi1,
                  np.ndarray[Dtype, ndim=2, mode="c"] xi2,
                  np.ndarray[Dtype, ndim=1, mode="c"] spar,
                  np.ndarray[Dtype, ndim=1, mode="c"] lpar,
                  list lpars):

    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(xi1)
    cdef np.ndarray al1 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray al2 = np.zeros((nx1,nx2),dtype=np.double)
    lq_nie(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>lpar.data,<Dtype *>al1.data,<Dtype *>al2.data)

    cdef np.ndarray al1s = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray al2s = np.zeros((nx1,nx2),dtype=np.double)

    llen = len(lpars)
    cdef np.ndarray lpars_i = np.zeros((nx1,nx2),dtype=np.double)
    for i in xrange(llen):
        lpars_i = lpars[i]
        lq_nie(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>lpars_i.data,<Dtype *>al1s.data,<Dtype *>al2s.data)
        al1 = al1+al1s
        al2 = al2+al2s

    cdef np.ndarray yi1 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray yi2 = np.zeros((nx1,nx2),dtype=np.double)
    yi1 = xi1-al1
    yi2 = xi2-al2

    cdef double dsx = xi1[1,1]-xi1[0,0]

    cdef np.ndarray a11 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray a12 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray a21 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray a22 = np.zeros((nx1,nx2),dtype=np.double)

    a12,a11 = np.gradient(al1,dsx)
    a22,a21 = np.gradient(al2,dsx)

    cdef np.ndarray mu = np.zeros((nx1,nx2),dtype=np.double)
    mu = 1.0/(1.0-(a11+a22)+a11*a22-a12*a21)

    cdef np.ndarray s_image = np.zeros((nx1,nx2),dtype=np.double)
    gauss_2d(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>spar.data,<Dtype *>s_image.data)

    cdef np.ndarray g_lensimage = np.zeros((nx1,nx2),dtype=np.double)
    gauss_2d(<Dtype *>yi1.data,<Dtype *>yi2.data,nx1,nx2,<Dtype *>spar.data,<Dtype *>g_lensimage.data)

    return s_image,g_lensimage,mu,yi1,yi2

#def lensed_images(np.ndarray[Dtype, ndim=2, mode="c"] xi1,
#                   np.ndarray[Dtype, ndim=2, mode="c"] xi2,
#                   np.ndarray[Dtype, ndim=1, mode="c"] spar,
#                   np.ndarray[Dtype, ndim=1, mode="c"] lpar,
#                   list lpars):
#
#    cdef int nx1
#    cdef int nx2
#    nx1,nx2 = np.shape(xi1)
#    cdef np.ndarray al1 = np.zeros((nx1,nx2),dtype=np.double)
#    cdef np.ndarray al2 = np.zeros((nx1,nx2),dtype=np.double)
#    lq_nie(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>lpar.data,<Dtype *>al1.data,<Dtype *>al2.data)
#
#    cdef np.ndarray al1s = np.zeros((nx1,nx2),dtype=np.double)
#    cdef np.ndarray al2s = np.zeros((nx1,nx2),dtype=np.double)
#
#    llen = len(lpars)
#    cdef np.ndarray lpars_i = np.zeros((nx1,nx2),dtype=np.double)
#    for i in xrange(llen):
#        lpars_i = lpars[i]
#        lq_nie(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>lpars_i.data,<Dtype *>al1s.data,<Dtype *>al2s.data)
#        al1 = al1+al1s
#        al2 = al2+al2s
#
#    cdef np.ndarray yi1 = np.zeros((nx1,nx2),dtype=np.double)
#    cdef np.ndarray yi2 = np.zeros((nx1,nx2),dtype=np.double)
#    yi1 = xi1-al1
#    yi2 = xi2-al2
#
#    cdef float dsx = xi1[1,1]-xi1[0,0]
#
#    cdef np.ndarray a11 = np.zeros((nx1,nx2),dtype=np.double)
#    cdef np.ndarray a12 = np.zeros((nx1,nx2),dtype=np.double)
#    cdef np.ndarray a21 = np.zeros((nx1,nx2),dtype=np.double)
#    cdef np.ndarray a22 = np.zeros((nx1,nx2),dtype=np.double)
#
#    a12,a11 = np.gradient(al1,dsx)
#    a22,a21 = np.gradient(al2,dsx)
#
#    cdef np.ndarray mu = np.zeros((nx1,nx2),dtype=np.double)
#    mu = 1.0/(1.0-(a11+a22)+a11*a22-a12*a21)
#
#    cdef np.ndarray s_image = np.zeros((nx1,nx2),dtype=np.double)
#    gauss_2d(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>spar.data,<Dtype *>s_image.data)
#
#    cdef np.ndarray g_lensimage = np.zeros((nx1,nx2),dtype=np.double)
#    gauss_2d(<Dtype *>yi1.data,<Dtype *>yi2.data,nx1,nx2,<Dtype *>spar.data,<Dtype *>g_lensimage.data)
#
#    return s_image,g_lensimage,mu,yi1,yi2

##@profile
#def lensed_images(xi1,xi2,spar,lpar,lpars):
#
#    dsx = xi1[1,1]-xi1[0,0]
#    al1,al2 = call_lq_nie(xi1,xi2,lpar)
#    for i in lpars:
#        al1s,al2s = call_lq_nie(xi1,xi2,i)
#        al1 = al1+al1s
#        al2 = al2+al2s
#
#    a12,a11 = np.gradient(al1,dsx)
#    a22,a21 = np.gradient(al2,dsx)
#
#    mu = 1.0/(1.0-(a11+a22)+a11*a22-a12*a21)
#
#    s_image = call_gauss_2d(xi1,xi2,spar)
#
#    yi1 = xi1-al1
#    yi2 = xi2-al2
#
#    g_lensimage = call_gauss_2d(yi1,yi2,spar)
#
#    return s_image,g_lensimage,mu,yi1,yi2

#def lens_images(xi1,xi2,gpar,gpars):
#
#    g_lens = call_gauss_2d(xi1,xi2,gpar)
#    for i in gpars:
#        g_lens_subs = call_gauss_2d(xi1,xi2,i)
#        g_lens = g_lens + g_lens_subs
#    return g_lens
#
#def mmbr_images(xi1,xi2,gpar,gpars):
#
#    g_lens = call_tophat_2d(xi1,xi2,gpar)
#    g_edge = call_find_critical_curve(g_lens)
#    for i in gpars:
#        g_lens_subs = call_tophat_2d(xi1,xi2,i)
#        g_edge_subs = call_find_critical_curve(g_lens_subs)
#        g_edge = g_edge+g_edge_subs
#    g_edge[g_edge>0.0] = 1.0
#    return g_edge

def lens_images(np.ndarray[Dtype, ndim=2, mode="c"] xi1,
                np.ndarray[Dtype, ndim=2, mode="c"] xi2,
                np.ndarray[Dtype, ndim=1, mode="c"] gpar,
                list gpars):

    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(xi1)
    cdef np.ndarray g_lens = np.zeros((nx1,nx2),dtype=np.double)
    gauss_2d(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>gpar.data,<Dtype *>g_lens.data)

    llen = len(gpars)
    cdef np.ndarray gpars_i = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray g_lens_subs = np.zeros((nx1,nx2),dtype=np.double)
    for i in xrange(llen):
        gpars_i = gpars[i]
        gauss_2d(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>gpars_i.data,<Dtype *>g_lens_subs.data)
        g_lens = g_lens + g_lens_subs

    return g_lens

def mmbr_images(np.ndarray[Dtype, ndim=2, mode="c"] xi1,
                np.ndarray[Dtype, ndim=2, mode="c"] xi2,
                np.ndarray[Dtype, ndim=1, mode="c"] gpar,
                list gpars):

    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(xi1)
    cdef np.ndarray g_lens = np.zeros((nx1,nx2),dtype=np.double)
    tophat_2d(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>gpar.data,<Dtype *>g_lens.data)
    cdef np.ndarray g_edge = np.zeros((nx1,nx2),dtype=np.double)
    find_critical_curve(<Dtype *>g_lens.data,nx1,nx2,<Dtype *>g_edge.data)

    llen = len(gpars)
    cdef np.ndarray gpars_i = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray g_lens_subs = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray g_edge_subs = np.zeros((nx1,nx2),dtype=np.double)
    for i in xrange(llen):
        gpars_i = gpars[i]
        tophat_2d(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>gpars_i.data,<Dtype *>g_lens_subs.data)
        find_critical_curve(<Dtype *>g_lens_subs.data,nx1,nx2,<Dtype *>g_edge_subs.data)
        g_edge = g_edge+g_edge_subs

    g_edge[g_edge>0.0] = 1.0
    return g_edge

#def find_caustic(lpar,lpars,critical,xi1,xi2,dsx):
#
#    cdef int nx1
#    cdef int nx2
#    nx1,nx2 = np.shape(xi1)
#    cdef double bsz;
#    bsz = dsx*nx1
#    yif1,yif2 = refine_critical(xi1,xi2,lpar,lpars,critical)
#    caustic = call_forward_cic(nx1,nx2,bsz,yif1.flatten(),yif2.flatten())
#    caustic[caustic>0]=1
#    return caustic

def all_about_lensing(np.ndarray[Dtype, ndim=2, mode="c"] xi1,
                      np.ndarray[Dtype, ndim=2, mode="c"] xi2,
                      np.ndarray[Dtype, ndim=1, mode="c"] spar,
                      np.ndarray[Dtype, ndim=1, mode="c"] lpar,
                      list lpars):

    cdef double dsx = xi1[1,1]-xi1[0,0]
    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(xi1)
    cdef np.ndarray al1 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray al2 = np.zeros((nx1,nx2),dtype=np.double)
    lq_nie(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>lpar.data,<Dtype *>al1.data,<Dtype *>al2.data)

    cdef np.ndarray al1s = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray al2s = np.zeros((nx1,nx2),dtype=np.double)

    cdef int i,j
    llen = len(lpars)
    cdef np.ndarray lpars_i = np.zeros((nx1,nx2),dtype=np.double)
    for i in xrange(llen):
        lpars_i = lpars[i]
        lq_nie(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>lpars_i.data,<Dtype *>al1s.data,<Dtype *>al2s.data)
        al1 = al1+al1s
        al2 = al2+al2s

    cdef np.ndarray yi1 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray yi2 = np.zeros((nx1,nx2),dtype=np.double)
    yi1 = xi1-al1
    yi2 = xi2-al2

    cdef np.ndarray a11 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray a12 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray a21 = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray a22 = np.zeros((nx1,nx2),dtype=np.double)

    a12,a11 = np.gradient(al1,dsx)
    a22,a21 = np.gradient(al2,dsx)

    #cdef np.ndarray mu = np.zeros((nx1,nx2),dtype=np.double)
    #mu = 1.0/(1.0-(a11+a22)+a11*a22-a12*a21)
    cdef np.ndarray imu = np.zeros((nx1,nx2),dtype=np.double)
    imu = (1.0-(a11+a22)+a11*a22-a12*a21)

#------------------------------------------------------------------------
    cdef np.ndarray s_image = np.zeros((nx1,nx2),dtype=np.double)
    gauss_2d(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>spar.data,<Dtype *>s_image.data)

    cdef np.ndarray g_lensimage = np.zeros((nx1,nx2),dtype=np.double)
    gauss_2d(<Dtype *>yi1.data,<Dtype *>yi2.data,nx1,nx2,<Dtype *>spar.data,<Dtype *>g_lensimage.data)
#------------------------------------------------------------------------
    cdef np.ndarray critical = np.zeros((nx1,nx2),dtype=np.double)
    find_critical_curve(<Dtype *>imu.data,nx1,nx2,<Dtype *>critical.data)
#------------------------------------------------------------------------
    yift1,yift2 = refine_critical(xi1,xi2,lpar,lpars,critical)
    cdef np.ndarray yif1 = np.zeros((len(yift1)),dtype=np.double)
    cdef np.ndarray yif2 = np.zeros((len(yift2)),dtype=np.double)
    yif1 = yift1
    yif2 = yift2

    cdef double bsz;
    bsz = dsx*nx1
    cdef int ylen = len(yif1.flatten())
    cdef np.ndarray img_in = np.array(np.ones(ylen),dtype=np.double)
    cdef np.ndarray caustic = np.zeros((nx1,nx2),dtype=np.double)
    forward_cic(<Dtype *>img_in.data,<Dtype *>yif1.data,<Dtype *>yif2.data,bsz,bsz,nx1,nx2,ylen,<Dtype *>caustic.data)
    caustic[caustic>0]=1

    return s_image,g_lensimage,critical,caustic.T
