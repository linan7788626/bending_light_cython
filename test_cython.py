#!/usr/bin/env python

import pygame
from pygame.locals import *
from sys import exit
from libv2_cv import *

def main():
    nnn = 512
    nnw = 512
    boxsize = 4.0
    dsx = boxsize/nnn
    xi1 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi2 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi1,xi2 = np.meshgrid(xi1,xi2)

    pygame.init()
    FPS = 15
    fpsClock = pygame.time.Clock()

    screen = pygame.display.set_mode((nnw, nnw), pygame.RESIZABLE, 32)
    #screen = pygame.display.set_mode((nnw, nnw), pygame.RESIZABLE| pygame.OPENGLBLIT | pygame.HWSURFACE | pygame.OPENGL | pygame.DOUBLEBUF)

    pygame.display.set_caption("Gravitational Lensing Toy")

    mouse_cursor = pygame.Surface((nnn,nnn))


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

    delta = 1e-8
    #lineThickness = 1

    #LeftButton=0

    ## Define some colors
    #BLACK = (0, 0, 0)
    #WHITE = (255, 255, 255)
    #GREEN = (0, 255, 0)
    #RED = (255, 0, 0)

    pygame.RESIZABLE


    while True:
        for event in pygame.event.get():
            if event.type == QUIT:
                exit()

        pos = pygame.mouse.get_pos()
        rotation=pygame.mouse.get_rel()
        buttonpress=pygame.mouse.get_pressed()
        keys = pygame.key.get_pressed()  #checking pressed keys

        if buttonpress[2] and keys[pygame.K_EQUALS]:
            xlcs = (pos[0]*nnn/nnw-nnn/2.0)*dsx
            ylcs = (pos[1]*nnn/nnw-nnn/2.0)*dsx
            lpar_sub = np.asarray([ylcs,xlcs,qls,rcs,res,phis])
            lpars.append(lpar_sub)
            gpar_sub = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis])
            gpars.append(gpar_sub)
        if buttonpress[2] and keys[pygame.K_MINUS]:
            if len(lpars) > 0:
                del lpars[-1]
            if len(gpars) > 0:
                del gpars[-1]
        for i in xrange(len(lpars)):
            kp = "K_"+str(i+2)
            kpv = getattr(pygame, kp)

            if rotation[0] and buttonpress[2] and keys[pygame.K_s] and keys[kpv] :

                lpars[i][3]=lpars[i][3]+rotation[0]*0.002
                if lpars[i][3] <= 0:
                    lpars[i][3] = delta
                if lpars[i][3] >= 1:
                    lpars[i][3] = 1.0-delta

                lpars[i][4]=lpars[i][4]-rotation[1]*0.005
                if lpars[i][4] <= 0:
                    lpars[i][4] = delta
                #-----------------------------------------------------
                gpars[i][3]=gpars[i][3]+rotation[0]*0.002
                if gpars[i][3] <= 0:
                    gpars[i][3] = delta
                if gpars[i][3] >= 1:
                    gpars[i][3] = 1.0-delta

                gpars[i][4]=gpars[i][4]-rotation[1]*0.005
                if gpars[i][4] <= 0:
                    gpars[i][4] = delta
                #-----------------------------------------------------

            if rotation[0] and buttonpress[2] and keys[pygame.K_w] and keys[kpv] :
                lpars[i][1]=lpars[i][1]+rotation[0]*0.01
                lpars[i][0]=lpars[i][0]+rotation[1]*0.01

                #-----------------------------------------------------
                gpars[i][1]=gpars[i][1]+rotation[0]*0.01
                gpars[i][0]=gpars[i][0]+rotation[1]*0.01
                #-----------------------------------------------------

            if rotation[0] and buttonpress[2] and keys[pygame.K_e] and keys[kpv] :
                lpars[i][5]=lpars[i][5]+rotation[0]

                lpars[i][2]=lpars[i][2]+rotation[1]*0.002
                if lpars[i][2] <= 0.3:
                    lpars[i][2] = 0.3
                if lpars[i][2] >= 1:
                    lpars[i][2] = 1.0-delta

                #-----------------------------------------------------
                gpars[i][5]=gpars[i][5]+rotation[0]

                gpars[i][2]=gpars[i][2]+rotation[1]*0.002
                if gpars[i][2] <= 0.3:
                    gpars[i][2] = 0.3
                if gpars[i][2] >= 1:
                    gpars[i][2] = 1.0-delta
                #-----------------------------------------------------

            #lpars[i] =  np.asarray([ylcs,xlcs,qls,rcs,res,phis])

        #----------------------------------------------------
        if rotation[0] and buttonpress[0] and keys[pygame.K_s]:
            gr_sig=gr_sig-rotation[1]*0.001
            if gr_sig <= 0:
                gr_sig = delta
            if gr_sig >= 1:
                gr_sig = 1.0

        if rotation[0] and buttonpress[0] and keys[pygame.K_w]:
            x += rotation[0]
            y += rotation[1]

        if rotation[0] and buttonpress[0] and keys[pygame.K_e]:
            gr_pa=gr_pa+rotation[0]

            gr_eq=gr_eq+rotation[1]*0.002
            if gr_eq <= 0.1:
                gr_eq = 0.1
            if gr_eq >= 1:
                gr_eq = 1.0-delta


        #----------------------------------------------------
        if rotation[0] and buttonpress[2] and keys[pygame.K_s] and keys[pygame.K_1]:

            rc0=rc0+rotation[0]*0.002

            if rc0 <= 0:
                rc0 = delta
            if rc0 >= 1:
                rc0 = 1.0-delta

            #l_sig0 = l_sig0-rotation[1]*0.001
            #if l_sig0 <= 0:
            #    l_sig0 = delta

            re0=re0-rotation[1]*0.005
            if re0 <= 0:
                re0 = delta

            l_sig0 = l_sig0-rotation[1]*0.005
            if l_sig0 <= 0:
                l_sig0 = delta

        if rotation[0] and buttonpress[2] and keys[pygame.K_w] and keys[pygame.K_1]:
            xlc0=xlc0+rotation[0]*0.01
            ylc0=ylc0+rotation[1]*0.01

        if rotation[0] and buttonpress[2] and keys[pygame.K_e] and keys[pygame.K_1]:
            phi0=phi0+rotation[0]

            ql0=ql0+rotation[1]*0.002
            if ql0 <= 0.3:
                ql0 = 0.3
            if ql0 >= 1:
                ql0 = 1.0-delta


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
        #----------------------------------------------

        g_lenses = lens_images(xi1,xi2,gpar,gpars)
        g_shapes = mmbr_images(xi1,xi2,gpar,gpars)

        baset[:,:,0] = g_shapes*255
        baset[:,:,1] = g_shapes*255
        baset[:,:,2] = g_shapes*255

        s_image,g_lensimage,critical,caustic = all_about_lensing(xi1,xi2,spar,lpar,lpars)

        base0[:,:,0] = g_lenses*255
        base0[:,:,1] = g_lenses*127
        base0[:,:,2] = g_lenses*0

        base1[:,:,0] = s_image*255
        base1[:,:,1] = s_image*255
        base1[:,:,2] = s_image*255

        base2[:,:,0] = g_lensimage*102
        base2[:,:,1] = g_lensimage*178
        base2[:,:,2] = g_lensimage*255

        base3[:,:,0] = critical*255
        base3[:,:,1] = critical*0
        base3[:,:,2] = critical*0

        base4[:,:,0] = caustic*0
        base4[:,:,1] = caustic*255
        base4[:,:,2] = caustic*0

        if keys[pygame.K_t]:
            base = baset+base1+base2

        elif keys[pygame.K_g]:
            wf = base0+base1+base2

            idx1 = wf>=base0
            idx2 = wf<base0
            base = base0*0
            base[idx1] = wf[idx1]
            base[idx2] = base0[idx2]
        else:
            wf = base1+base2+base3+base4
            base = wf

        pygame.surfarray.blit_array(mouse_cursor,base)


        screen.blit(pygame.transform.scale(mouse_cursor,(nnw,nnw)), (0, 0))

        #font=pygame.font.SysFont(None,30)
        #text = font.render("( "+str(x)+", "+str(-y)+" )", True, (255, 255, 255))
        #screen.blit(text,(10, 10))
        pygame.display.update()
        #pygame.display.flip()
        fpsClock.tick(FPS)

if __name__ == '__main__':
    main()
