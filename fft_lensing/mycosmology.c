#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "mycosmology.h"
//-----------------------------------------------------------------------------
//print '---------------------------------------------------------'
//print 'Cosmology : h			=', h
//print 'Cosmology : OmegaM0		=', Om0
//print 'Cosmology : OmegaL0		=', Ol0
//print 'Lensing   : zl			=', zl
//print 'Lensing   : zs			=', zs
//print 'Constant  : arcsec per rad	=', apr
//print 'Constant  : Hubble costant	= %.2e  km/s/(Mpc/h)' %H0
//print 'Constant  : Hubble distance	= %.2e  Mpc/h' %Dh()
//print 'Constant  : Hubble time		= %.2e  Gyr/h' %Th()
//print 'Constant  : G			= %.2e  Mpc/Msun/(km/s)^2' % G
//print 'Constant  : c			= %.2e  km/s' % vc
//print 'Constant  : Msun 		= %.2e  kg' % Msun
//print 'Constant  : pc			= %.2e  m' % pc
//print '---------------------------------------------------------'
//print 'For zl = %.3f:' % (zl)
//print 'Omega_m(z)			= %.2e' % Omz(zl)
//print 'Omega_Lambda(z)			= %.2e' % Olz(zl)
//print 'H(z)				= %.2e km/s/(Mpc/h)' % Hz(zl)
//print 'Comoving L.O.S. distance 	= %.2e Mpc/h' % Dc(zl)
//print 'Angular diameter distance	= %.2e Mpc/h' % Da(zl)
//print 'Luminosity distance		= %.2e Mpc/h' % Dl(zl)
//print 'arcsec per Mpc/h 		= %.2e ' % (1.0/(Dl(zl)/apr))
//print 'Critical Density 		= %.2e (Msun/h)/(Mpc/h)^3' % rho_crit(zl)
//print 'Average Density			= %.2e (Msun/h)/(Mpc/h)^3' % rho_bar(zl)
//print 'Lookback time			= %.2e Gyr/h' % Tl(zl)
//print 'Age of the Universe		= %.2e Gyr/h' % age_z(zl)
//print 'Scale Factor			= %.2e ' % a(zl)
//print 'Distance modulus 		= %.2e ' % DistMod(zl)
//print '---------------------------------------------------------'
//print 'For zs = %.3f:' % (zs)
//print 'Omega_m(z)			= %.2e' % Omz(zs)
//print 'Omega_Lambda(z)			= %.2e' % Olz(zs)
//print 'H(z)				= %.2e km/s/(Mpc/h)' % Hz(zs)
//print 'Comoving L.O.S. distance 	= %.2e Mpc/h' % Dc(zs)
//print 'Angular diameter distance	= %.2e Mpc/h' % Da(zs)
//print 'Luminosity distance		= %.2e Mpc/h' % Dl(zs)
//print 'arcsec per Mpc/h 		= %.2e ' % (1.0/(Dl(zs)/apr))
//print 'Critical Density 		= %.2e (Msun/h)/(Mpc/h)^3' % rho_crit(zs)
//print 'Average Density			= %.2e (Msun/h)/(Mpc/h)^3' % rho_bar(zs)
//print 'Lookback time			= %.2e Gyr/h' % Tl(zs)
//print 'Age of the Universe		= %.2e Gyr/h' % age_z(zs)
//print 'Scale Factor			= %.2e ' % a(zs)
//print 'Distance modulus 		= %.2e ' % DistMod(zs)
//print '---------------------------------------------------------'
//print 'Critical surface density 	= %.2e (Msun/h)/(Mpc/h)^2' % sigma_crit(zl,zs)
//print 'Dals				= %.2e Mpc/h' % (Da2(zl,zs))
//print '---------------------------------------------------------'
//-----------------------------------------------------------------------------
double tzf(double x){
	double res;
	res = ezf(x)/(1.0+x);
	return res;
}
double tz(double x,void *params){
	double alpha;
	alpha = *(double *) params;
	alpha = 1.0;
	double res;
	res = tzf(x)*alpha;
	return res;
}
double Hz(double x){
	double res;
	res = H0/ezf(x);
	return res;
}
double scale_factor(double x){
	double res;
	res = 1.0/(1.0+x);
	return res;
}
double Dh(){
	double res;
	res = vc/H0;
	return res;
}
double Th(){
	double res;
	res = 1.0/H0*1e3*kpc/Gyr;
	return res;
}
double Tl(double x){
	double res;
	double abr;
	int_l1p(0.0,x,tz,&res,&abr);
	return res*Th();
}
double age_z(double x){
	double res;
	double abr;
	int_su1p(x,tz,&res,&abr);
	return res*Th();
}
double Dc(double x){
	double res;
	double abr;
	int_l1p(0.0,x,ez,&res,&abr);
	return res*Dh();
}
double Dm(double x){
	double res;
	double sOk0;
	sOk0 = sqrt(fabs(Ok0));
	if(Ok0 > 0){
		res = Dh()/sOk0*sinh(sOk0*Dc(x)/Dh());
	}
	else if(Ok0 == 0.0){
		res = Dc(x);
	}
	else{
		res = Dh()/sOk0*sin(sOk0*Dc(x)/Dh());
	}
	return res;
}
double Dl(double x){
	double res;
	res = Dc(x)*(1.0+x);
	return res;
}
double Dp(double x1,double x2){
	double res;
	double abr;
	int_l1p(x1,x2,tz,&res,&abr);
	return res*Dh();
}
double DistMod(double x){
	double res;
	res = 5.0*log10(Dl(x)*(Mpc/h)/pc/10.0);
	return res;
}
double ez(double x,void *params){
	double alpha;
	alpha = *(double *) params;
	alpha = 1.0;
	double res;
	res = 1.0/sqrt(Ol0+Ok0*(1.0+x)*(1.0+x)+Om0*(1.0+x)*(1.0+x)*(1.0+x))*alpha;
	return res;
}
double ezf(double x){
	double res;
	res = 1.0/sqrt(Ol0+Ok0*(1.0+x)*(1.0+x)+Om0*(1.0+x)*(1.0+x)*(1.0+x));
	return res;
}
double ez_integral(double x){
	double res;
	double abr;
	int_l1p(0.0,x,ez,&res,&abr);
	return res;
}
double Omz(double x){
	double res = 0.0;
	res = Om0*pow((1.0+x),3.0)*ezf(x)*ezf(x);
	return res;
}
double Olz(double x){
	double res = 0.0;
	res = Ol0*pow((1.0+x),0.0)*ezf(x)*ezf(x);
	return res;
}
double Okz(double x){
	double res = 0.0;
	res = Ok0*pow((1.0+x),2.0)*ezf(x)*ezf(x);
	return res;
}
double rho_crit(double x){
	double res = 0.0;
	//rho_crit0 = 1.0e-29*1.0e-33*2.937999e+73
	res =  rho_crit0/(ezf(x)*ezf(x));
	return res;
}
double rho_bar(double x){
	double res = 0.0;
	res = rho_crit(x)*Omz(x);
	return res;
}
int retdz12(double x1,double x2,double *dz12){
	*dz12 = vc/(1.0+x2)/H0*(ez_integral(x2)-ez_integral(x1));
	return 0;
}
double Da(double x){
	double res;
	retdz12(0.0,x,&res);
	return res;
}
double Da2(double x1,double x2){
	double res;
	retdz12(x1,x2,&res);
	return res;
}
double dv(double z){
	double res = 0.0;
	double ov = 0.0;
        ov = 1.0/Omz(z)-1.0;
        res = 18.8*M_PI*M_PI*(1.0+0.4093*pow(ov,0.9052));
        return res;
}
//----------------------------------------------------------------------------
double frand (void){
	double value;
	value =((double)rand()/(RAND_MAX));
	return value;
}
//----------------------------------------------------------------------------
double sign(double x){
	double value;
	value = x/fabs(x);
	return value;
}
//-------------------------------1 parameter----------------------------------
double sp_func(double x, void *params){
	double alpha = *(double *) params;
	alpha = 1.0;
	double res;
	//res = log(alpha*x)/sqrt(x);
	res = x*x*alpha;
	return res;
}
//----------------------------------------------------------------------------
int dif1p(double x,double (*fp)(double,void*),double * result,double * abserr){
	double alpha =1.0;
	gsl_function F;
	F.function = fp;
	F.params = &alpha;

 	gsl_deriv_central (&F, x, 1e-8, result, abserr);
 	//gsl_deriv_forward (&F, x, 1e-8, result, abserr);
 	//gsl_deriv_backward (&F, x, 1e-8, result, abserr);

	return 0;
}
//----------------------------------------------------------------------------
int int_u1p(double (*fp)(double, void*),double * result,double * error){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	double alpha = 1.0;

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagi (&F, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_sd1p(double a,double (*fp)(double, void*),double * result,double * error){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	double alpha = 1.0;

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagil (&F, a, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_su1p(double a,double (*fp)(double, void*),double * result,double * error){

	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	double alpha = 1.0;

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagiu (&F, a, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_l1p(double a,double b,double (*fp)(double, void*),double * result,double * error){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	double alpha = 1.0;

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qags (&F, a, b, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int f_root1p(double a,double b,double (*fp)(double, void*),double * result){
	int status;
	int iter = 0, max_iter = 1000;
	double r = 0.0;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;

	double alpha = 1.0;

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s,&F,a,b);

	do{
	    iter++;
	    status = gsl_root_fsolver_iterate(s);
	    r = gsl_root_fsolver_root(s);
	    a = gsl_root_fsolver_x_lower(s);
	    b = gsl_root_fsolver_x_upper(s);
	    status = gsl_root_test_interval(a,b,0,0.001);
	}while (status == GSL_CONTINUE && iter < max_iter);
	*result = r;

	gsl_root_fsolver_free (s);

	return 0;
}
//-------------------------------2 parameters---------------------------------
double dp_func(double x, void *params)
{
	double alpha = *(double *) params;
	double res;
	//res = log(alpha*x)/sqrt(x);
	res = alpha*x*x;
	return res;
}
//----------------------------------------------------------------------------
int dif2p(double x,double alpha,double (*fp)(double,void*),double * result,double * abserr){

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

 	gsl_deriv_central (&F, x, 1e-8, result, abserr);
 	//gsl_deriv_forward (&F, x, 1e-8, result, abserr);
 	//gsl_deriv_backward (&F, x, 1e-8, result, abserr);

  	return 0;
}
//----------------------------------------------------------------------------
int int_u2p(double alpha,double (*fp)(double, void*),double * result,double * error){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagi (&F, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_sd2p(double alpha,double a,double (*fp)(double, void*),double * result,double * error){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagil (&F, a, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_su2p(double alpha,double a,double (*fp)(double, void*),double * result,double * error){

	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagiu (&F, a, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_l2p(double alpha,double a,double b,double (*fp)(double, void*),double * result,double * error){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qags (&F, a, b, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int f_root2p(double alpha,double a,double b,double (*fp)(double, void*),double * result){
	int status;
	int iter = 0, max_iter = 1000;
	double r = 0;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s,&F,a,b);

	do{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		a = gsl_root_fsolver_x_lower(s);
		b = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(a,b,0,0.001);
	}while (status == GSL_CONTINUE && iter < max_iter);
	*result = r;

	gsl_root_fsolver_free (s);

       return 0;
}
//-------------------------------3 parameters---------------------------------
double tp_func(double x, void *params){
	double res = 0.0;
	struct f_params_2 *p
	  = (struct f_params_2 *) params;

	double a = p->a;
	double b = p->b;

	res = a*x*x+b;

	return res;
}
//----------------------------------------------------------------------------
int dif3p(double x,double p1,double p2,double (*fp)(double,void*),
double * result,double * abserr){

	struct f_params_2 alpha = {p1,p2};
	gsl_function F;
	F.function = fp;
	F.params = &alpha;

 	gsl_deriv_central (&F, x, 1e-8, result, abserr);
 	//gsl_deriv_forward (&F, x, 1e-8, result, abserr);
 	//gsl_deriv_backward (&F, x, 1e-8, result, abserr);

	return 0;
}
//----------------------------------------------------------------------------
int int_u3p(double p1,double p2,double (*fp)(double, void*),
double * result,double * error){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	struct f_params_2 alpha = {p1,p2};
	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagi (&F, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_sd3p(double p1,double p2,double a,double (*fp)(double, void*),
double * result,double * error){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	struct f_params_2 alpha = {p1,p2};
	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagil (&F, a, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_su3p(double p1,double p2,double a,double (*fp)(double, void*),
double * result,double * error){

	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	struct f_params_2 alpha = {p1,p2};
	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagiu (&F, a, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_l3p(double p1,double p2,double a,double b,double (*fp)(double, void*),
double * result,double * error){

	struct f_params_2 alpha = {p1,p2};

	gsl_integration_workspace * w
	 = gsl_integration_workspace_alloc (1000);

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qags (&F, a, b, 0, 1e-7, 1000,w, result, error);
	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int f_root3p(double p1,double p2,double a,double b,double (*fp)(double, void*),double * result){
	struct f_params_2 alpha = {p1,p2};
	int status;
	int iter = 0, max_iter = 1000;
	double r = 0.0;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s,&F,a,b);

	do{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		a = gsl_root_fsolver_x_lower(s);
		b = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(a,b,0,0.001);
	}while (status == GSL_CONTINUE && iter < max_iter);
	*result=r;

	gsl_root_fsolver_free (s);

	return 0;
}
//-------------------------------4 parameters---------------------------------
double qp_func(double x, void *params){
	double res = 0.0;
	struct f_params_3 *p
	  = (struct f_params_3 *) params;

	double a = p->a;
	double b = p->b;
	double c = p->c;

	res = a*pow(x,b)+c;
	return res;
}
//----------------------------------------------------------------------------
int dif4p(double x,double p1,double p2,double p3,double (*fp)(double,void*),
double * result,double * abserr){

	struct f_params_3 alpha = {p1,p2,p3};
	gsl_function F;
	F.function = fp;
	F.params = &alpha;

 	gsl_deriv_central (&F, x, 1e-8, result, abserr);
 	//gsl_deriv_forward (&F, x, 1e-8, result, abserr);
 	//gsl_deriv_backward (&F, x, 1e-8, result, abserr);

	return 0;
}
//----------------------------------------------------------------------------
int int_u4p(double p1,double p2,double p3,double (*fp)(double, void*),
double * result,double * error){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	struct f_params_3 alpha = {p1,p2,p3};
	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagi (&F, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_sd4p(double p1,double p2,double p3,double a,double (*fp)(double, void*),
double * result,double * error){
	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	struct f_params_3 alpha = {p1,p2,p3};
	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagil (&F, a, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_su4p(double p1,double p2,double p3,double a,double (*fp)(double, void*),
double * result,double * error){

	gsl_integration_workspace * w
	  = gsl_integration_workspace_alloc (1000);

	struct f_params_3 alpha = {p1,p2,p3};
	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	gsl_integration_qagiu (&F, a, 0, 1e-7, 1000,w, result, error);

	gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int int_l4p(double p1,double p2,double p3,double a,double b,double (*fp)(double, void*),
double * result,double * error){

        struct f_params_3 alpha = {p1,p2,p3};

        gsl_integration_workspace * w
         = gsl_integration_workspace_alloc (1000);

        gsl_function F;
        F.function = fp;
        F.params = &alpha;

        gsl_integration_qags (&F, a, b, 0, 1e-7, 1000,w, result, error);
        gsl_integration_workspace_free (w);

	return 0;
}
//----------------------------------------------------------------------------
int f_root4p(double p1,double p2,double p3,double a,double b,double (*fp)(double, void*),
double * result){
	struct f_params_3 alpha = {p1,p2,p3};
	int status;
	int iter = 0, max_iter = 1000;
	double r = 0.0;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;

	gsl_function F;
	F.function = fp;
	F.params = &alpha;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s,&F,a,b);

	do{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		a = gsl_root_fsolver_x_lower(s);
		b = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(a,b,0,0.001);
	}while (status == GSL_CONTINUE && iter < max_iter);
	*result=r;

	gsl_root_fsolver_free (s);

	return 0;
}
//----------------------------------------------------------------------------
double sigma_crit(double z1,double z2){
	double res = 0.0;
	res = vc*vc*Da(z2)/(4.0*M_PI*M_G*Da(z1)*Da2(z1,z2));
	return res;
}
////----------------------------------------------------------------------------
//void main(){
//	double res,abr;
//	double a,b,c;
//
//	//----------------------1p-----------------------
//	//dif1p(1.0,sp_func, &res,&abr);
//	//int_l1p(0.0,3.0,sp_func,&res,&abr);
//	//f_root1p(0.0,1.0,sp_func-0.25,&res);
//	//----------------------2p-----------------------
//	//a = 1.0;
//	//dif2p(1.0,1.0,dp_func,&res,&abr);
//	//int_l2p(1.0,0.0,3.0,dp_func,&res,&abr);
//	//f_root2p(1.0,0.0,1.0,dp_func,&res);
//	//printf("%lf %lf\n",res,abr);
//	//----------------------3p-----------------------
//	//a = 1.0;
//	//b = 0.0
//	//dif3p(1.0,0.0,1.0,tp_func,&res,&abr);
//	//int_l3p(1.0,0.0,0.0,3.0,tp_func,&res,&abr);
//	//f_root3p(1.0,-0.25,0.0,1.0,tp_func,&res);
//	//----------------------4p-----------------------
//	a = 1.0;
//	b = 2.0;
//	c = 0.0;
//	dif4p(a,b,c,1.0,qp_func,&res,&abr);
//	int_l4p(a,b,c,0.0,3.0,qp_func,&res,&abr);
//	f_root4p(a,b,c-0.25,0.0,1.0,qp_func,&res);
//	//-----------------------------------------------
//	printf("%lf %lf\n",res,abr);
//}
