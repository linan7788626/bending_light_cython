#define h 0.71
#define Om0 0.27
#define Ol0 0.73
#define Otot Om0+Ol0
#define Ok0 0.0
//#define w 0.0-1.0
#define rho_crit0 2.78e11 //M_sun Mpc^-3 *h*h
#define rho_bar0 rho_crit0*Om0
#define sigma8 0.801 //wmap 7th
//----------------------------------------------------------------------------
#define apr 206269.43				//1/1^{''}
#define vc 2.9970e5				//km/s
#define M_G 4.3e-9				//(Mpc/h)^1 (Msun/h)^-1 (km/s)^2
#define H0 100.0				//km/s/(Mpc/h)
#define pc 3.085677e16				//meter
#define kpc 3.085677e19				//meter
#define Mpc 3.085677e22				//meter
#define Msun 1.98892e30				//kg
#define yr 31536000.0/365.0*365.25		//second
#define Gyr yr*1e9				//second
//double M_G = 6.67259e-11;			//m^3/kg/s^2
//double 1Gly = 9.461e26cm = 9.461e24m = 9.461e21km;
//double 1Mpc = 3.08568e24cm = 3.261566ly;
//----------------------------------------------------------------------------
double frand (void);
double sign(double x);
double ez(double x,void *params);
double ezf(double x);
double ez_integral(double x);
double tzf(double x);
double tz(double x,void *params);
double Hz(double x);
double scale_factor(double x);
double Dh();
double Th();
double Tl(double x);
double age_z(double x);
double Dc(double x);
double Dm(double x);
double Dl(double x);
double Dp(double x1,double x2);
double DistMod(double x);
double Omz(double x);
double Olz(double x);
double Okz(double x);
double rho_crit(double x);
double rho_bar(double x);
double Da(double x);
double Da2(double x1,double x2);
double dv(double z);
int retdz12(double x1,double x2,double *dz12);
double sigma_crit(double x1,double x2);
//----------------------------------------------------------------------------
struct f_params_2{double a,b;};
struct f_params_3{double a,b,c;};
//-------------------------------Def GSL fun----------------------------------
double sp_func(double x, void *params);
double dp_func(double x, void *params);
double tp_func(double x, void *params);
double qp_func(double x, void *params);
//-------------------------------1 parameter----------------------------------
int dif1p(double x,double (*fp)(double,void*),double * result,double * abserr);
int int_u1p(double (*fp)(double, void*),double * result,double * error);
int int_sd1p(double a,double (*fp)(double, void*),double * result,double * error);
int int_su1p(double a,double (*fp)(double, void*),double * result,double * error);
int int_l1p(double a,double b,double (*fp)(double, void*),double * result,double * error);
int f_root1p(double a,double b,double (*fp)(double, void*),double * result);
//-------------------------------2 parameters---------------------------------
int dif2p(double x,double alpha,double (*fp)(double,void*),double * result,double * abserr);
int int_u2p(double alpha,double (*fp)(double, void*),double * result,double * error);
int int_sd2p(double alpha,double a,double (*fp)(double, void*),double * result,double * error);
int int_su2p(double alpha,double a,double (*fp)(double, void*),double * result,double * error);
int int_l2p(double alpha,double a,double b,double (*fp)(double, void*),double * result,double * error);
int f_root2p(double alpha,double a,double b,double (*fp)(double, void*),double * result);
//-------------------------------3 parameters---------------------------------
int dif3p(double x,double p1,double p2,double (*fp)(double,void*),
double * result,double * abserr);
int int_u3p(double p1,double p2,double (*fp)(double, void*),
double * result,double * error);
int int_sd3p(double p1,double p2,double a,double (*fp)(double, void*),
double * result,double * error);
int int_su3p(double p1,double p2,double a,double (*fp)(double, void*),
double * result,double * error);
int int_l3p(double p1,double p2,double a,double b,double (*fp)(double, void*),
double * result,double * error);
int f_root3p(double p1,double p2,double a,double b,double (*fp)(double, void*),double * result);
//-------------------------------4 parameters---------------------------------
int dif4p(double x,double p1,double p2,double p3,double (*fp)(double,void*),
double * result,double * abserr);
int int_u4p(double p1,double p2,double p3,double (*fp)(double, void*),
double * result,double * error);
int int_sd4p(double p1,double p2,double p3,double a,double (*fp)(double, void*),
double * result,double * error);
int int_su4p(double p1,double p2,double p3,double a,double (*fp)(double, void*),
double * result,double * error);
int int_l4p(double p1,double p2,double p3,double a,double b,double (*fp)(double, void*),
double * result,double * error);
int f_root4p(double p1,double p2,double p3,double a,double b,double (*fp)(double, void*),
double * result);
