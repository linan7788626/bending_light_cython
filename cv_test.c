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
