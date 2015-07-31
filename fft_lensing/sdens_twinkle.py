#!/usr/bin/env python
import numpy as np
import libfft_lensing as lf
import pylab as pl
import scipy.signal as ss

#def re0_sigma(sigma):
#    cv = 3e5
#    Dds = 1.0
#    Ds = 2.0
#    res = 4.0*np.pi*(sigma/cv)**2.0*Dds/Ds
#    return res

def hfunc(x1,x2,rcore,qe):
    res = np.sqrt(qe*qe*(rcore*rcore+x1*x1)+x2*x2)
    return res


def nie_phi(x1,x2,re0,rcore,qe):

    res0 = re0/np.sqrt(1-qe*qe)
    al1 = res0*np.arctan(x1*np.sqrt(1-qe*qe)/(hfunc(x1,x2)+rcore))
    al2 = res0*np.arctanh(x2*np.sqrt(1-qe*qe)/(hfunc(x1,x2)+rcore*qe*qe))

    res1 = x1*al1+x2*al2
    res2 = re0*rcore*np.log(np.sqrt((hfunc(x1,x2)+rcore)**2.0+(1-qe*qe)*x1*x1))
    res = res1-res2
    return res

def nie_alphas(x1,x2,re0,rcore,qe):
    res0 = re0/np.sqrt(1-qe*qe)
    al1 = np.arctan(x1*np.sqrt(1-qe*qe)/(hfunc(x1,x2)+rcore))
    al2 = np.arctanh(x2*np.sqrt(1-qe*qe)/(hfunc(x1,x2)+rcore*qe*qe))
    return res0*al1,res0*al2

def nie_kappa(x1,x2,re0,rcore,qe):
    res = re0/(2.0*np.sqrt(qe*qe*(rcore+x1*x1)+x2*x2))
    return res

def lpar_nie_kappa(xi1,xi2,lpar):
    xc1 = lpar[0]
    xc2 = lpar[1]
    b = lpar[2]
    s = lpar[3]
    q = lpar[4]
    rot = lpar[5]

    x1,x2 = xy_rotate(xi1,xi2,xc1,xc2,rot)

    wx = np.sqrt(q*q*(x1*x1+s*s)+x2*x2)

    kappa = b/(2.0*wx)

    return kappa

def nie_mu(x1,x2,re0,rcore,qe):
    res = 1.0/(1.0-re0/hfunc(x1,x2)-re0*re0*rcore/(hfunc(x1,x2)*((hfunc(x1,x2)+rcore)**2+(1-qe*qe)*x1*x1)))
    return res

def new_nie_all(xi1,xi2,lpar):
    xc1 = lpar[0]
    xc2 = lpar[1]
    b = lpar[2]
    s = lpar[3]
    q = lpar[4]
    rot = lpar[5]

    dsx = xi1[1,1]-xi1[0,0]

    x1,x2 = xy_rotate(xi1,xi2,xc1,xc2,rot)

    wx = np.sqrt(q*q*(x1*x1+s*s)+x2*x2)

    a1 = b/np.sqrt(1-q*q)*np.arctan(x1*np.sqrt(1-q*q)/(wx+s))
    a2 = b/np.sqrt(1-q*q)*np.arctanh(x2*np.sqrt(1-q*q)/(wx+q*q*s))

    hx = np.sqrt((wx+s)**2.0+(1-q*q)*x1*x1)
    phi = x1*a1+x2*a2-b*s*np.log(hx)+b*q*s*np.log((1+q)*s)

    ai2,ai1 = np.gradient(phi,dsx)

    return phi,ai1,ai2#,kappa,mu,y1,y2

def nie_all(xi1,xi2,xc1,xc2,b,s,q,rot,ys1,ys2):

    x1,x2 = xy_rotate(xi1,xi2,xc1,xc2,rot)

    wx = np.sqrt(q*q*(x1*x1+s*s)+x2*x2)

    al1 = b/np.sqrt(1-q*q)*np.arctan(x1*np.sqrt(1-q*q)/(wx+s))
    al2 = b/np.sqrt(1-q*q)*np.arctanh(x2*np.sqrt(1-q*q)/(wx+q*q*s))

    kappa = b/(2.0*wx)

    hx = np.sqrt((wx+s)**2.0+(1-q*q)*x1*x1)
    phi = x1*al1+x2*al2-b*s*np.log(hx)+b*q*s*np.log((1+q)*s)

    Kc = 1.0
    #Kc = (1.0+zl)/c*(Dl*Ds/Dls)
    td = Kc*(0.5*((al1)**2.0+(al2)**2.0)-phi)
    #td = Kc*(0.5*((x1-ys1)**2.0+(x2-ys2)**2.0)-phi)

    y1 = xi1-al1
    y2 = xi2-al2

    #y1,y2 = xy_rotate(y1,y2,xc1,xc2,-rot)

#------------------------------------------------------------------
    demon1 = ((wx+s)**2+(1.0-q*q)*x1*x1)*wx
    demon2 = (((wx+q*q*s)**2-(1.0-q*q)*x2*x2)*wx)
    y11 = 1-b*(wx*(wx+s)-q*q*x1*x1)/demon1
    y22 = 1-b*(wx*(wx+q*q*s)-x2*x2)/demon2
    y12 = -b*x1*x2/demon1
    y21 = -b*x1*x2*q*q/demon2

    mu = 1.0/(y11*y22-y12*y21)

    return phi,td,al1,al2#,kappa,mu,y1,y2

def multiple_nie_all(xi1,xi2,lpars_list):
    phi = xi1*0.0
    al1 = xi1*0.0
    al2 = xi1*0.0
    for i in lpars_list:
        phi_tmp,al1_tmp,al2_tmp = lpar_nie_all(xi1,xi2,i)
        phi = phi + phi_tmp
        al1 = al1 + al1_tmp
        al2 = al2 + al2_tmp

    return phi,al1,al2

def multiple_new_nie_all(xi1,xi2,lpars_list):
    phi = xi1*0.0
    al1 = xi1*0.0
    al2 = xi1*0.0
    for i in lpars_list:
        phi_tmp,al1_tmp,al2_tmp = new_nie_all(xi1,xi2,i)
        phi = phi + phi_tmp
        al1 = al1 + al1_tmp
        al2 = al2 + al2_tmp

    return phi,al1,al2

def lpar_nie_all(xi1,xi2,lpar):

    xc1 = lpar[0]
    xc2 = lpar[1]
    b = lpar[2]
    s = lpar[3]
    q = lpar[4]
    rot = lpar[5]

    x1,x2 = xy_rotate(xi1,xi2,xc1,xc2,rot)

    wx = np.sqrt(q*q*(x1*x1+s*s)+x2*x2)

    al1 = b/np.sqrt(1-q*q)*np.arctan(x1*np.sqrt(1-q*q)/(wx+s))
    al2 = b/np.sqrt(1-q*q)*np.arctanh(x2*np.sqrt(1-q*q)/(wx+q*q*s))

    hx = np.sqrt((wx+s)**2.0+(1-q*q)*x1*x1)
    phi = x1*al1+x2*al2-b*s*np.log(hx)+b*q*s*np.log((1+q)*s)

    return phi,al1,al2

#--------------------------------------------------------------------
def green_iso(Nc,dsx):
    boxsize = Nc*dsx
    x1 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,Nc)+0.5*dsx
    x2 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,Nc)+0.5*dsx
    x1,x2 = np.meshgrid(x1,x2)

    epsilon= 0.0001
    r = np.sqrt(x1**2.0+x2**2.0+epsilon**2.0)
    res = 1.0/np.pi*np.log(r)
    return res

def alphas_iso(Nc,dsx):
    boxsize = Nc*dsx
    x1 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,Nc)+0.5*dsx
    x2 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,Nc)+0.5*dsx
    x1,x2 = np.meshgrid(x1,x2)

    epsilon= 0.0001
    r = np.sqrt(x1**2.0+x2**2.0+epsilon**2.0)
    res1 = x1/(r**2.0*np.pi)
    res2 = x2/(r**2.0*np.pi)
    return res1,res2

def fft_lensing_signals(kappac,green_in,dsx):

    phi = ss.fftconvolve(kappac,green_in,mode='same')*(dsx*dsx)

    alpha2,alpha1 = np.gradient(phi,dsx)
    phi12,phi11 = np.gradient(alpha1,dsx)
    phi22,phi21 = np.gradient(alpha2,dsx)

    kappas = (phi11+phi22)*0.5
    #shear1 = (phi11-phi22)*0.5
    #shear2 = (phi12+phi21)*0.5

    mu = 1.0/(1.0-phi11-phi22+phi11*phi22-phi12*phi21)

    Kc = 1.0
    #Kc = (1.0+zl)/c*(Dl*Ds/Dls)
    td = Kc*(0.5*((alpha1)**2.0+(alpha2)**2.0)-phi)

    nx,ny = np.shape(phi)

    return phi[nx/4:3*nx/4,ny/4:3*ny/4],alpha1[nx/4:3*nx/4,ny/4:3*ny/4],alpha2[nx/4:3*nx/4,ny/4:3*ny/4],td[nx/4:3*nx/4,ny/4:3*ny/4],mu[nx/4:3*nx/4,ny/4:3*ny/4],kappas[nx/4:3*nx/4,ny/4:3*ny/4]

def lensed_images(xi1,xi2,yi1,yi2,gpar):

    g_image = gauss_2d(xi1,xi2,gpar)
    g_lensimage = gauss_2d(yi1,yi2,gpar)

    return g_image,g_lensimage

def xy_rotate(x, y, xcen, ycen, phi):

    phirad = np.deg2rad(phi)
    xnew = (x - xcen) * np.cos(phirad) + (y - ycen) * np.sin(phirad)
    ynew = (y - ycen) * np.cos(phirad) - (x - xcen) * np.sin(phirad)
    return (xnew,ynew)

def gauss_2d(x, y, par):

    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    r_ell_sq = ((xnew**2)*par[4] + (ynew**2)/par[4]) / np.abs(par[1])**2
    return par[0] * np.exp(-0.5*r_ell_sq)

def find_critical_curve(mu):
    rows,cols = np.indices(np.shape(mu))
    cdtn = np.sign(mu)*(np.sign(mu[rows-1,cols])+np.sign(mu[rows,cols-1])+np.sign(mu[(rows+1)%len(rows),cols])+np.sign(mu[rows,(cols+1)%len(cols)]))

    res = mu*0
    res[cdtn<4] = 1
    res[cdtn>=4] = 0

    return res

def lens_galaxies(xi1,xi2,glpar):

    g_lens = gauss_2d(xi1,xi2,glpar)

    return g_lens

@profile
def main():

    nnn = 512
    boxsize = 4.0
    zl = 0.1
    zs = 1.0
    p_mass = 1.0
    dsx = boxsize/nnn
    xi1 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi2 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi1,xi2 = np.meshgrid(xi1,xi2)
    #----------------------------------------------------
    # lens parameters for main halo
    xlc1 = 0.0
    xlc2 = 0.0
    ql0 = 0.999999999999
    rc0 = 0.000000000001
    re0 = 1.0
    phi0 = 0.0
    lpar = np.asarray([xlc1, xlc2, re0, rc0, ql0, phi0])

    lpars_list = []
    lpars_list.append(lpar)
    #----------------------------------------------------
    sdens = lpar_nie_kappa(xi1,xi2,lpar)
    #pii,aii1,aii2 = multiple_new_nie_all(xi1,xi2,lpars_list)

    phi,phi1,phi2,td = lf.call_all_about_lensing(sdens,nnn,zl,zs,p_mass,dsx)
    #print np.shape(phi)

    #phi2,phi1 = np.gradient(phi,dsx)

    #phi12,phi11 = np.gradient(phi2,dsx)
    #phi22,phi21 = np.gradient(phi1,dsx)
    #kappac = 0.5*(phi11+phi22)
    ##----------------------------------------------------
    ## lens parameters for main halo
    #xls1 = 0.7
    #xls2 = 0.8
    #qls = 0.999999999999
    #rcs = 0.000000000001
    #res = 0.5
    #phis = 0.0
    #lpars = np.asarray([xls1, xls2, res, rcs, qls, phis])
    #lpars_list.append(lpars)

    #sdens = lpar_nie_kappa(xi1,xi2,lpar)

    #pii,pii1,pii2 = multiple_new_nie_all(xi1,xi2,lpars_list)

    #phi,alpha1,alpha2,td = lf.call_all_about_lensing(sdens,nnn,zl,zs,p_mass,dsx)

    #phi12,phi11 = np.gradient(alpha2,dsx)
    #phi22,phi21 = np.gradient(alpha1,dsx)
    #kappai = 0.5*(phi11+phi22)
#--------------------------------------------------------------------
    #sdens_pad = np.zeros((nnn*2,nnn*2))
    #sdens_pad[nnn/2:nnn/2*3,nnn/2:nnn/2*3] = sdens
    #green_in = green_iso(nnn*2,dsx)
    #phi,alpha1,alpha2,td,mu,kappas = fft_lensing_signals(sdens_pad,green_in,dsx)
#--------------------------------------------------------------------

    #Kc = 1.0
    ##Kc = (1.0+zl)/c*(Dl*Ds/Dls)
    #td = Kc*(0.5*((phi1)**2.0+(phi2)**2.0)-pii)
    #tdi = Kc*(0.5*((pii1)**2.0+(pii2)**2.0)-pii)

    ##levels = [-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6]
    #pl.figure()
    #pl.contourf(xi1,xi2,np.log10(kappai))
    #pl.colorbar()
    #pl.figure()
    #pl.contourf(xi1,xi2,np.log10(kappas))
    #pl.colorbar()

    #pl.figure()
    #pl.imshow((np.sqrt(phi1*phi1+phi2*phi2)-np.sqrt(aii1*aii1+aii2*aii2))/np.sqrt(phi1*phi1+phi2*phi2), aspect='auto', cmap=pl.get_cmap(pl.cm.jet),vmin=-0.01, vmax=0.01 )
    #pl.colorbar()

    ##levels = [-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0]
    #pl.figure()
    ##pl.contourf(xi1,xi2,np.log10(kappac),levels)
    ##pl.colorbar()
    ##print np.max((kappac-sdens)/sdens),np.mean((kappac-sdens)/sdens)
    #pl.imshow((kappac-sdens)/sdens, aspect='auto', cmap=pl.get_cmap(pl.cm.jet),vmin=-0.01, vmax=0.01)
    #pl.colorbar()
    ##pl.contour(xi1,xi2,np.log10(kappac),levels,colors=['k',])
    ##pl.contour(xi1,xi2,np.log10(sdens),levels,colors=['r',])
    ##pl.show()

    ##levels = [3.0,2.5,2.0,1.5,1.0,0.5,0.0,-0.5]
    #pl.figure()
    #pl.contour(xi1,xi2,phi,levels,colors=['k',])
    #pl.contour(xi1,xi2,pii,levels,colors=['r',])
    ##pl.imshow((phi-(np.median(phi-pii))-pii)/phi, aspect='auto', cmap=pl.get_cmap(pl.cm.jet),vmin=-0.01, vmax=0.01)
    ##pl.imshow((phi-pii)/phi, aspect='auto', cmap=pl.get_cmap(pl.cm.jet))#,vmin=-0.01, vmax=0.01)
    #pl.colorbar()
    #pl.show()

    #levels = [-1.6,-1.2,-0.8,0.4,0.0,0.4,0.8,1.2,1.6]
    #pl.figure()
    ##pl.contour(xi1,xi2,np.sqrt(alpha2**2.0+alpha1**2.0),levels,colors=['k',])
    ##pl.contour(xi1,xi2,np.sqrt(phi2**2.0+phi1**2.0),levels,colors=['r',])
    #pl.contour(xi1,xi2,np.sqrt(alpha2**2.0+alpha1**2.0),levels)
    #pl.contour(xi1,xi2,np.sqrt(phi2**2.0+phi1**2.0),levels)
    #pl.colorbar()

    #levels = [-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0]
    #pl.figure()
    #pl.contourf(xi1,xi2,td,levels)#,colors=['k',])
    #pl.colorbar()
    #pl.figure()
    #pl.contourf(xi1,xi2,tdi,levels)#,colors=['r',])
    #pl.colorbar()


    ##----------------------------------------------------
    #data = np.fromfile("./out_iso.bin",dtype=np.double)
    #data = data.reshape((np.sqrt(len(data)),np.sqrt(len(data))))
    #pl.figure()
    #pl.contourf(data)
    #pl.colorbar()

    #pl.show()


if __name__ == '__main__':
    main()


