#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int sign(double x) {
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

double deg2rad(double pha) {
	double res = 0;
	res = pha*M_PI/180.0;
	return res;
}

void forward_cic(double *cic_in,double *x_in,double *y_in,double bsx,double bsy,int nx,int ny,int np,double *cic_out) {
    double dx = bsx/nx;
    double dy = bsy/ny;
    double xc = bsx/2.0;
    double yc = bsy/2.0;
    double wx,wy;
    double xp,yp,zp;

    int i;
    int ip,jp;

    for (i=0;i<np;i++) {
        xp = (x_in[i]+xc)/dx-0.5;
        yp = (y_in[i]+yc)/dy-0.5;
		zp = cic_in[i];

        ip = (int)xp;
        jp = (int)yp;

		if (ip<0||ip>(nx-2)||jp<0||jp>(ny-2)) continue;
        wx = 1.0-(xp-(double)ip);
        wy = 1.0-(yp-(double)jp);

        cic_out[ip*ny+jp] += wx*wy*zp;
        cic_out[ip*ny+(jp+1)] += wx*(1.0-wy)*zp;
        cic_out[(ip+1)*ny+jp] += (1.0-wx)*wy*zp;
        cic_out[(ip+1)*ny+(jp+1)] += (1.0-wx)*(1.0-wy)*zp;
    }
}

//--------------------------------------------------------------------
void lanczos_diff_2_tag(double *m1, double *m2, double *m11, double *m12, double *m21, double *m22, double Dcell, int Ncc, int dif_tag) {
    int i_m3,i_p3,j_m3,j_p3,i_m2,i_p2,j_m2,j_p2,i_m1,j_m1,i_p1,j_p1,i,j;
    int index;

    for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		i_m1 = i-1;i_m2 = i-2;i_m3 = i-3;
		i_p1 = i+1;i_p2 = i+2;i_p3 = i+3;
		j_m1 = j-1;j_m2 = j-2;j_m3 = j-3;
		j_p1 = j+1;j_p2 = j+2;j_p3 = j+3;

		if (j_m3<0||j_p3>(Ncc-1)||i_m3<0||i_p3>(Ncc-1)) continue;

        //if (i==0) {i_m1 = Ncc-1;i_m2 = Ncc-2;i_m3 = Ncc-3;}
        //else if (i==1) {i_m1 = 0;i_m2 = Ncc-1;i_m3 = Ncc-2;}
        //else if (i==2) {i_m1 = 1;i_m2 = 0;i_m3 = Ncc-1;}
        //else {i_m1 = i-1;i_m2 = i-2;i_m3 = i-3;}
        //if (j==0) {j_m1 = Ncc-1;j_m2 = Ncc-2;j_m3 = Ncc-3;}
        //else if (j==1) {j_m1 = 0;j_m2 = Ncc-1;j_m3 = Ncc-2;}
        //else if (j==2) {j_m1 = 1;j_m2 = 0;j_m3 = Ncc-1;}
        //else {j_m1 = j-1;j_m2 = j-2;j_m3 = j-3;}
        //if (i==Ncc-1) {i_p1 = 0;i_p2 = 1;i_p3 = 2;}
        //else if (i==Ncc-2) {i_p1 = Ncc-1;i_p2 = 0;i_p3 = 1;}
        //else if (i==Ncc-3) {i_p1 = Ncc-2;i_p2 = Ncc-1;i_p3 = 0;}
        //else {i_p1 = i+1;i_p2 = i+2;i_p3 = i+3;}
        //if (j==Ncc-1) {j_p1 = 0;j_p2 = 1;j_p3 = 2;}
        //else if (j==Ncc-2) {j_p1 = Ncc-1;j_p2 = 0;j_p3 = 1;}
        //else if (j==Ncc-2) {j_p1 = Ncc-2;j_p2 = Ncc-1;j_p3 = 0;}
        //else {j_p1 = j+1;j_p2 = j+2;j_p3 = j+3;}

        index = i*Ncc+j;
        if (dif_tag==-1) {
            m11[index] = (m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])/(2.0*Dcell);
            m12[index] = (m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])/(2.0*Dcell);
            m21[index] = (m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])/(2.0*Dcell);
            m22[index] = (m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])/(2.0*Dcell);
        }

        if (dif_tag==0) {
            m11[index] = (m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])*2.0/3.0/Dcell
            - (m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])/12.0/Dcell;
            m22[index] = (m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])*2.0/3.0/Dcell
            - (m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])/12.0/Dcell;
            m21[index] = (m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])*2.0/3.0/Dcell
            - (m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])/12.0/Dcell;
            m12[index] = (m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])*2.0/3.0/Dcell
            - (m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])/12.0/Dcell;
        }

        if (dif_tag==1) {
            m11[index] =(1.0*(m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])
                         + 2.0*(m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])
                         + 3.0*(m1[i_p3*Ncc+j]-m1[i_m3*Ncc+j]))/(28.0*Dcell);
            m22[index] =(1.0*(m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])
                         + 2.0*(m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])
                         + 3.0*(m2[i*Ncc+j_p3]-m2[i*Ncc+j_m3]))/(28.0*Dcell);
            m12[index] =(1.0*(m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])
                         + 2.0*(m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])
                         + 3.0*(m1[i*Ncc+j_p3]-m1[i*Ncc+j_m3]))/(28.0*Dcell);
            m21[index] =(1.0*(m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])
                         + 2.0*(m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])
                         + 3.0*(m2[i_p3*Ncc+j]-m2[i_m3*Ncc+j]))/(28.0*Dcell);
        }

        if (dif_tag==2) {
            m11[index] = (5.0*(m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])
                          + 4.0*(m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])
                          + 1.0*(m1[i_p3*Ncc+j]-m1[i_m3*Ncc+j]))/(32.0*Dcell);
            m22[index] = (5.0*(m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])
                          + 4.0*(m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])
                          + 1.0*(m2[i*Ncc+j_p3]-m2[i*Ncc+j_m3]))/(32.0*Dcell);
            m12[index] = (5.0*(m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])
                          + 4.0*(m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])
                          + 1.0*(m1[i*Ncc+j_p3]-m1[i*Ncc+j_m3]))/(32.0*Dcell);
            m21[index] = (5.0*(m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])
                          + 4.0*(m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])
                          + 1.0*(m2[i_p3*Ncc+j]-m2[i_m3*Ncc+j]))/(32.0*Dcell);
        }

        if (dif_tag==3) {
            m11[index] = (58.0*(m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])
                          + 67.0*(m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])
                          + 22.0*(m1[i_p3*Ncc+j]-m1[i_m3*Ncc+j]))/(252.0*Dcell);
            m22[index] = (58.0*(m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])
                          + 67.0*(m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])
                          - 22.0*(m2[i*Ncc+j_p3]-m2[i*Ncc+j_m3]))/(252.0*Dcell);
            m12[index] = (58.0*(m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])
                          + 67.0*(m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])
                          - 22.0*(m1[i*Ncc+j_p3]-m1[i*Ncc+j_m3]))/(252.0*Dcell);
            m21[index] = (58.0*(m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])
                          + 67.0*(m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])
                          - 22.0*(m2[i_p3*Ncc+j]-m2[i_m3*Ncc+j]))/(252.0*Dcell);
        }
    }
}

void xy_rotate(double *x1_in,double *x2_in,int nx1,int nx2,double xc1,double xc2,double pha,double *x1_out,double *x2_out) {
    double phirad = deg2rad(pha);
	int i,j,index;
	for (i=0;i<nx1;++i) for (j=0;j<nx2;++j){
		index = i*nx2+j;
		x1_out[index] = (x1_in[index] - xc1)*cos(phirad)+(x2_in[index]-xc2)*sin(phirad);
    	x2_out[index] = (x2_in[index] - xc2)*cos(phirad)-(x1_in[index]-xc1)*sin(phirad);
	}
}

void gauss_2d(double *x1,double *x2,int nx1,int nx2,double *par,double *res) {
    //gpars = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis]);
	double *x1new = (double *)malloc(nx1*nx2*sizeof(double));
	double *x2new = (double *)malloc(nx1*nx2*sizeof(double));
    xy_rotate(x1,x2,nx1,nx2,par[0], par[1],par[5],x1new,x2new);
	int i,j,index;
	for (i=0;i<nx1;++i) for (j=0;j<nx2;++j) {
		index = i*nx2+j;
    	res[index] = par[3]*exp(-0.5*((x1new[index]*x1new[index])*par[2]+(x2new[index]*x2new[index])/par[2])/(par[4]*par[4]));
	}
	free(x1new);
	free(x2new);
}

void tophat_2d(double *x1,double *x2, int nx1, int nx2,double *par,double *res) {
    //gpars = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis]);
	double *x1new = (double *)malloc(nx1*nx2*sizeof(double));
	double *x2new = (double *)malloc(nx1*nx2*sizeof(double));
    xy_rotate(x1,x2,nx1,nx2,par[0], par[1], par[5],x1new,x2new);
	int i,j,index;
	double r_ell;
	for (i=0;i<nx1;++i) for (j=0;j<nx2;++j){
		index = i*nx2+j;
    	r_ell = sqrt((x1new[index]*x1new[index])*par[2]+(x2new[index]*x2new[index])/par[2]);
		if (r_ell>=par[4]) {
			res[index] = -1.0;
		}
		else {
			res[index] = 10000.0;
		}
	}
	free(x1new);
	free(x2new);
}

void lq_nie(double *x1,double *x2,int nx1,int nx2,double *lpar,double *alpha1,double *alpha2) {

    double xc1 = lpar[0];
    double xc2 = lpar[1];
    double q   = lpar[2];
    double rc  = lpar[3];
    double re  = lpar[4];
    double pha = lpar[5];

    double phirad = deg2rad(pha);
    double cosa = cos(phirad);
    double sina = sin(phirad);
	double phi,a1,a2;
	double *xt1 = (double *)malloc(sizeof(double)*nx1*nx2);
	double *xt2 = (double *)malloc(sizeof(double)*nx1*nx2);

	int i,j,index;
	for (i = 0;i<nx1;++i) for (j = 0;j<nx2;++j) {
		index = i*nx2+j;
		xt1[index] = (x1[index]-xc1)*cosa+(x2[index]-xc2)*sina;
    	xt2[index] = (x2[index]-xc2)*cosa-(x1[index]-xc1)*sina;
		phi = sqrt(xt2[index]*xt2[index]+xt1[index]*q*xt1[index]*q+rc*rc);

    	a1 = sqrt(q)/sqrt(1.0-q*q)*atan(sqrt(1.0-q*q)*xt1[index]/(phi+rc/q));
    	a2 = sqrt(q)/sqrt(1.0-q*q)*atanh(sqrt(1.0-q*q)*xt2[index]/(phi+rc*q));

    	alpha1[index] = (a1*cosa-a2*sina)*re;
    	alpha2[index] = (a2*cosa+a1*sina)*re;
	}
	free(xt1);
	free(xt2);
}

void find_critical_curve(double *mu,int nx,int ny,double* res) {

	int i,j,index,sign_t=0;
	int im1,ip1,jm1,jp1;
	for (i = 0; i < nx; ++i) for (j = 0; j < ny; ++j) {
		index = i*ny+j;
		im1 = i-1;
		ip1 = i+1;
		jm1 = j-1;
		jp1 = j+1;

		if (im1<0||jm1<0||ip1>(nx-1)||jp1>(ny-1)) continue;

		sign_t = sign(mu[index])*(sign(mu[im1*ny+j])
								 +sign(mu[i*ny+jm1])
								 +sign(mu[ip1*ny+j])
								 +sign(mu[i*ny+jp1]));
		if (sign_t < 4) {
			res[index] = 1.0;
		}
		else {
			res[index] = 0.0;
		}
	}
}
void tot_lq(double *x1, double *x2,int nx1,int nx2,double *lpar, int npars, double *lpars, int nsubs, double *y1, double *y2) {

	int i,j,k,l,index;
    double * al1 = (double *)malloc(sizeof(double)*nx1*nx2);
    double * al2 = (double *)malloc(sizeof(double)*nx1*nx2);
    lq_nie(x1,x2,nx1,nx2,lpar,al1,al2);

    double * als1 = (double *)malloc(sizeof(double)*nx1*nx2);
    double * als2 = (double *)malloc(sizeof(double)*nx1*nx2);

    double * lpars_i = (double *)malloc(sizeof(double)*npars);
	for (i = 0; i < nsubs; ++i) {
		for (j = 0; j < npars;++j) {
			lpars_i[j] = lpars[i*npars+j];
		}
        lq_nie(x1,x2,nx1,nx2,lpars_i,als1,als2);
		for (k = 0; k < nx1; k++) for (l = 0; l < nx2; l++) {
			index = k*nx2+l;
			al1[index] = al1[index]+als1[index];
			al2[index] = al2[index]+als2[index];
		}
	}

	for (k = 0; k < nx1; k++) for (l = 0; l < nx2; l++) {
		index = k*nx2+l;
		y1[index] = x1[index]-al1[index];
		y2[index] = x2[index]-al2[index];
	}
    free(al1);
    free(al2);
    free(als1);
    free(als2);
    free(lpars_i);
}

//    int clen = 0;
//
//	for (i = 0; i < nx1; ++i) for (j = 0; j < nx2; ++j){
//		index = i*nx2+j;
//		if (critical[index]>0) {
//			clen = clen+1;
//		}
//	}
void refine_critical(double * xi1,double * xi2,int nx1,int nx2,double * lpar,int npars,double * lpars, int nsubs,double * critical,int clen, int nfiner, double * yi1,double *yi2) {

	int i,j,k=0,m,n,index;
    double dsx = xi1[nx2+1]-xi1[0];
    double dsf = dsx/nfiner;

    double * xt1 = (double *)malloc(clen*nfiner*nfiner*sizeof(double));
    double * xt2 = (double *)malloc(clen*nfiner*nfiner*sizeof(double));

	for (i = 0; i < nx1; ++i) for (j = 0; j < nx2; ++j){
		index = i*nx2+j;
		if (critical[index]>0) {
			for (m = 0; m < nfiner; ++m) for (n = 0; n < nfiner; ++n){
        	    xt1[k*nfiner*nfiner+m*nfiner+n] = xi1[index]+(dsf*(1-nfiner)*0.5)+dsf*m;
        	    xt2[k*nfiner*nfiner+m*nfiner+n] = xi2[index]+(dsf*(1-nfiner)*0.5)+dsf*n;
			}
			k = k+1;
		}
	}
	tot_lq(xt1,xt2,clen,nfiner*nfiner,lpar,npars,lpars,nsubs,yi1,yi2);
	free(xt1);
	free(xt2);
}
void lens_images(double *xi1,double *xi2,int nx1,int nx2,double *gpar,int npars,double *gpars,int nsubs,double *g_lens) {
	int i,j,k,l,index;
    gauss_2d(xi1,xi2,nx1,nx2,gpar,g_lens);
    double * gpars_i = (double *)malloc(npars*sizeof(double));
    double * g_lens_subs = (double *)malloc(nx1*nx2*sizeof(double));
	for (i = 0; i < nsubs; ++i) {
		for (j = 0; j < npars; ++j) {
			gpars_i[j] = gpars[i*npars+j];
		}
		gauss_2d(xi1,xi2,nx1,nx2,gpars_i,g_lens_subs);
		for (k = 0; k < nx1; ++k) for (l = 0; l < nx2; ++l){
			index = k*nx2+l;
			g_lens[index] = g_lens[index] + g_lens_subs[index];
		}
	}
    free(gpars_i);
    free(g_lens_subs);
}

void mmbr_images(double *xi1,double *xi2,int nx1,int nx2,double *gpar,int npars,double *gpars,int nsubs,double *g_edge) {

	int i,j,k,l,index;
    double * g_lens = (double *)malloc(nx1*nx2*sizeof(double));
    tophat_2d(xi1,xi2,nx1,nx2,gpar,g_lens);
    find_critical_curve(g_lens,nx1,nx2,g_edge);

    double * gpars_i = (double *)malloc(npars*sizeof(double));
    double * g_lens_subs = (double *)malloc(nx1*nx2*sizeof(double));
    double * g_edge_subs = (double *)malloc(nx1*nx2*sizeof(double));

	for (i = 0; i < nsubs; ++i) {
		for (j = 0; j < npars; ++j) {
			gpars_i[j] = gpars[i*npars+j];
		}
		tophat_2d(xi1,xi2,nx1,nx2,gpars_i,g_lens_subs);
        find_critical_curve(g_lens_subs,nx1,nx2,g_edge_subs);
		for (k = 0; k < nx1; ++k) for (l = 0; l < nx2; ++l){
			index = k*nx2+l;
			g_edge[index] = g_edge[index] + g_edge_subs[index];
		}
	}
    free(gpars_i);
    free(g_lens_subs);

	for (k = 0; k < nx1; ++k) for (l = 0; l < nx2; ++l){
		index = k*nx2+l;
		if (g_edge[index]>0.0) {
			g_edge[index] = 1.0;
		}
	}
}
void all_about_lensing(double *xi1,double *xi2,int nx1,int nx2,double * spar,double * lpar,int npars,double * lpars,int nsubs,double *s_image,double *g_lensimage,double *critical,double *caustic){
	int i,j,k,l,index;
    double * al1 = (double *)malloc(sizeof(double)*nx1*nx2);
    double * al2 = (double *)malloc(sizeof(double)*nx1*nx2);
    lq_nie(xi1,xi2,nx1,nx2,lpar,al1,al2);

    double * als1 = (double *)malloc(sizeof(double)*nx1*nx2);
    double * als2 = (double *)malloc(sizeof(double)*nx1*nx2);

    double * lpars_i = (double *)malloc(sizeof(double)*npars);
	for (i = 0; i < nsubs; ++i) {
		for (j = 0; j < npars;++j) {
			lpars_i[j] = lpars[i*npars+j];
		}
        lq_nie(xi1,xi2,nx1,nx2,lpars_i,als1,als2);
		for (k = 0; k < nx1; k++) for (l = 0; l < nx2; l++) {
			index = k*nx2+l;
			al1[index] = al1[index]+als1[index];
			al2[index] = al2[index]+als2[index];
		}
	}
    free(als1);
    free(als2);
    free(lpars_i);
//------------------------------------------------------------------------
    double dsx = xi1[nx2+1]-xi1[0];
    double * a11 = (double *)malloc(nx1*nx2*sizeof(double));
    double * a12 = (double *)malloc(nx1*nx2*sizeof(double));
    double * a21 = (double *)malloc(nx1*nx2*sizeof(double));
    double * a22 = (double *)malloc(nx1*nx2*sizeof(double));

	lanczos_diff_2_tag(al1,al2,a21,a22,a11,a12,dsx,nx1,-1);

    double * imu = (double *)malloc(nx1*nx2*sizeof(double));
	for (k = 0; k < nx1; k++) for (l = 0; l < nx2; l++) {
		index = k*nx2+l;
		imu[index] = (1.0-(a11[index]+a22[index])+a11[index]*a22[index]-a12[index]*a21[index]);
	}

	free(a11);
	free(a12);
	free(a21);
	free(a22);

    find_critical_curve(imu,nx1,nx2,critical);
	free(imu);

//------------------------------------------------------------------------
    double * yi1 = (double *)malloc(sizeof(double)*nx1*nx2);
    double * yi2 = (double *)malloc(sizeof(double)*nx1*nx2);
	for (k = 0; k < nx1; k++) for (l = 0; l < nx2; l++) {
		index = k*nx2+l;
		yi1[index] = xi1[index]-al1[index];
		yi2[index] = xi2[index]-al2[index];
	}
    free(al1);
    free(al2);
    gauss_2d(xi1,xi2,nx1,nx2,spar,s_image);
    gauss_2d(yi1,yi2,nx1,nx2,spar,g_lensimage);
	free(yi1);
	free(yi2);
//------------------------------------------------------------------------
    int clen = 0;
	int nfiner = 1;

	for (i = 0; i < nx1; ++i) for (j = 0; j < nx2; ++j){
		index = i*nx2+j;
		if (critical[index]>0) {
			clen = clen+1;
		}
	}
	int ylen = clen*nfiner*nfiner;

    double * yif1 = (double *)malloc(ylen*sizeof(double));
    double * yif2 = (double *)malloc(ylen*sizeof(double));

	refine_critical(xi1,xi2,nx1,nx2,lpar,npars,lpars,nsubs,critical,clen,nfiner,yif1,yif2);

    double bsz;
    bsz = dsx*nx1;
    double * img_in = (double *)malloc(ylen*sizeof(double));

	for (i = 0; i < ylen; ++i) {
		img_in[i] = 1.0;
	}
    forward_cic(img_in,yif1,yif2,bsz,bsz,nx1,nx2,ylen,caustic);

	free(yif1);
	free(yif2);
	free(img_in);

	for (i = 0; i < nx1; ++i) for (j = 0; j < nx2; ++j){
		index = i*nx2+j;
		if (caustic[index]>0) {
			caustic[index] = 1;
		}
	}
}

//int main(int argc, const char *argv[])
//{
//	all_about_lensing(double *xi1,double *xi2,int nx1,int nx2,double * spar,double * lpar,int npars,double * lpars,int nsubs,double *s_image,double *g_lensimage,double *critical,double *caustic){
//
//	return 0;
//}
//
