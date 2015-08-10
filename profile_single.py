#!/usr/bin/env python

from libv4_cv import *
#from libv3_cv import np, call_lens_images,call_mmbr_images,call_all_about_lensing
#from libv3_cv import *

@profile
def main():
    nnn = 1024
    boxsize = 4.0
    dsx = boxsize/nnn
    xi1 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi2 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi1,xi2 = np.meshgrid(xi1,xi2)

    baset = np.zeros((nnn,nnn,3),'uint8')
    base0 = np.zeros((nnn,nnn,3),'uint8')
    base1 = np.zeros((nnn,nnn,3),'uint8')
    base2 = np.zeros((nnn,nnn,3),'uint8')
    base3 = np.zeros((nnn,nnn,3),'uint8')
    base4 = np.zeros((nnn,nnn,3),'uint8')

    #----------------------------------------------------
    # parameters of source
    x = 0
    y = 0
    #step = 1
    gr_sig = 0.02
    gr_eq = 1.0
    gr_pa = 0.0

    #----------------------------------------------------
    # lens parameters for mainhalo
    xlc0 = 0.0
    ylc0 = 0.0
    ql0 = 0.7
    rc0 = 0.1
    re0 = 1.0
    phi0 = 0.0
    #----------------------------------------------------
    # lens parameters for subhalo
    xlcs = 0.7
    ylcs = 0.77
    qls = 0.99999
    rcs = 0.00001
    res = 0.05
    phis = 0.0
    #----------------------------------------------------
    # parameters of NIE model (lens model, deflection angles)
    #----------------------------------------------------
    # 1, y position of center
    # 2, x position of center
    # 3, minor-to-major axis ratio
    # 4, size of flat core
    # 5, Einstein radius (lensing strength)
    # 6, major-axis position angle (degrees) c.c.w. from y axis
    lpar_sub = np.asarray([ylcs,xlcs,qls,rcs,res,phis])
    lpars = [lpar_sub]

    #----------------------------------------------------
    # luminosity parameters for mainbhalo
    ap0 = 1.0
    l_sig0 = 0.5
    #----------------------------------------------------
    # luminosity parameters for subhalo
    aps = 0.4
    l_sigs = 0.05
    #----------------------------------------------------
    # Parameters of Gaussian model (luminosity distribution of lenses)
    #----------------------------------------------------
    # 1, peak brightness value
    # 2, Gaussian "sigma" (i.e., size)
    # 3, y position of center
    # 4, x position of center
    # 5, minor-to-major axis ratio
    # 6, major-axis position angle (degrees) c.c.w. from y axis

    gpars_sub = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis])
    gpars = [gpars_sub]
    #---------------------------------------------------

    lpar = np.asarray([ylc0,xlc0,ql0,rc0,re0,phi0])
    gpar = np.asarray([ylc0,xlc0,ql0,ap0,l_sig0,phi0])

    #----------------------------------------------
    #parameters of source galaxies.
    #----------------------------------------------
    g_amp = 1.0         # peak brightness value
    g_sig = gr_sig          # Gaussian "sigma" (i.e., size)
    g_ycen = y*2.0/nnn  # y position of center
    g_xcen = x*2.0/nnn  # x position of center
    g_axrat = gr_eq       # minor-to-major axis ratio
    g_pa = gr_pa          # major-axis position angle (degrees) c.c.w. from y axis
    spar = np.asarray([g_ycen,g_xcen,g_axrat,g_amp,g_sig,g_pa])
    spars = []
    #----------------------------------------------

    g_lenses = call_lens_images(xi1,xi2,gpar,gpars)
    g_shapes = call_mmbr_images(xi1,xi2,gpar,gpars)
    s_image,g_lensimage,criticals,caustics = call_all_about_lensing(xi1,xi2,spar,spars,lpar,lpars)

    baset[:,:,0] = g_shapes*255
    baset[:,:,1] = g_shapes*255
    baset[:,:,2] = g_shapes*255

    base0[:,:,0] = g_lenses*255
    base0[:,:,1] = g_lenses*127
    base0[:,:,2] = g_lenses*0

    base1[:,:,0] = s_image*255
    base1[:,:,1] = s_image*255
    base1[:,:,2] = s_image*255

    base2[:,:,0] = g_lensimage*102
    base2[:,:,1] = g_lensimage*178
    base2[:,:,2] = g_lensimage*255

    base3[:,:,0] = criticals*255
    base3[:,:,1] = criticals*0
    base3[:,:,2] = criticals*0

    base4[:,:,0] = caustics*0
    base4[:,:,1] = caustics*255
    base4[:,:,2] = caustics*0

    wf = base1+base2+base3+base4
    #print np.max(wf)

if __name__ == '__main__':
    main()
