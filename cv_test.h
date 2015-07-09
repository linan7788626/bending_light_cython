int sign(double x);
double deg2rad(double pha);
void xy_rotate(double *x_in,double *y_in,int nx,int ny,double xcen,double ycen,double pha,double *x_out,double *y_out);
void gauss_2d(double *x1,double *x2, int nx1,int nx2,double *par,double *res);
void tophat_2d(double *x1,double *x2, int nx1,int nx2,double *par,double *res);
void lq_nie(double *x1,double *x2,int nx1,int nx2,double *lpar,double *alpha1,double *alpha2);
//void refine_critical();
void find_critical_curve(double *mu,int nx1,int nx2,double* res);
