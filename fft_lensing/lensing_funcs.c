#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "mycosmology.h"
//--------------------------------------------------------------------
void lanczos_diff_1_tag(double *mi, double *m1, double *m2, double Dcell, int Ncc, int dif_tag) {
	long i_m3,i_p3,j_m3,j_p3,i_m2,i_p2,j_m2,j_p2,i_m1,j_m1,i_p1,j_p1,i,j;
	long index;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if (i==0) {i_m1 = Ncc-1;i_m2 = Ncc-2;i_m3 = Ncc-3;}
		else if (i==1) {i_m1 = 0;i_m2 = Ncc-1;i_m3 = Ncc-2;}
		else if (i==2) {i_m1 = 1;i_m2 = 0;i_m3 = Ncc-1;}
		else {i_m1 = i-1;i_m2 = i-2;i_m3 = i-3;}
		if (j==0) {j_m1 = Ncc-1;j_m2 = Ncc-2;j_m3 = Ncc-3;}
		else if (j==1) {j_m1 = 0;j_m2 = Ncc-1;j_m3 = Ncc-2;}
		else if (j==2) {j_m1 = 1;j_m2 = 0;j_m3 = Ncc-1;}
		else {j_m1 = j-1;j_m2 = j-2;j_m3 = j-3;}
		if (i==Ncc-1) {i_p1 = 0;i_p2 = 1;i_p3 = 2;}
		else if (i==Ncc-2) {i_p1 = Ncc-1;i_p2 = 0;i_p3 = 1;}
		else if (i==Ncc-3) {i_p1 = Ncc-2;i_p2 = Ncc-1;i_p3 = 0;}
		else {i_p1 = i+1;i_p2 = i+2;i_p3 = i+3;}
		if (j==Ncc-1) {j_p1 = 0;j_p2 = 1;j_p3 = 2;}
		else if (j==Ncc-2) {j_p1 = Ncc-1;j_p2 = 0;j_p3 = 1;}
		else if (j==Ncc-2) {j_p1 = Ncc-2;j_p2 = Ncc-1;j_p3 = 0;}
		else {j_p1 = j+1;j_p2 = j+2;j_p3 = j+3;}

		index = i*Ncc+j;

        if (dif_tag==-1) {
            m1[index] = (mi[i_p1*Ncc+j]-mi[i_m1*Ncc+j])/(2.0*Dcell);
            m2[index] = (mi[i*Ncc+j_p1]-mi[i*Ncc+j_m1])/(2.0*Dcell);
        }

		if (dif_tag==0) {
			m1[index] = (mi[i_p1*Ncc+j]-mi[i_m1*Ncc+j])*2.0/3.0/Dcell
			    	  - (mi[i_p2*Ncc+j]-mi[i_m2*Ncc+j])/12.0/Dcell;
			m2[index] = (mi[i*Ncc+j_p1]-mi[i*Ncc+j_m1])*2.0/3.0/Dcell
					  - (mi[i*Ncc+j_p2]-mi[i*Ncc+j_m2])/12.0/Dcell;
		}

		if (dif_tag==1) {
			m1[index] =(1.0*(mi[i_p1*Ncc+j]-mi[i_m1*Ncc+j])
			  		  + 2.0*(mi[i_p2*Ncc+j]-mi[i_m2*Ncc+j])
			   		  + 3.0*(mi[i_p3*Ncc+j]-mi[i_m3*Ncc+j]))/(28.0*Dcell);
			m2[index] =(1.0*(mi[i*Ncc+j_p1]-mi[i*Ncc+j_m1])
			  		  + 2.0*(mi[i*Ncc+j_p2]-mi[i*Ncc+j_m2])
			   		  + 3.0*(mi[i*Ncc+j_p3]-mi[i*Ncc+j_m3]))/(28.0*Dcell);
		}

		if (dif_tag==2) {
			m1[index] =(5.0*(mi[i_p1*Ncc+j]-mi[i_m1*Ncc+j])
			  	   	  + 4.0*(mi[i_p2*Ncc+j]-mi[i_m2*Ncc+j])
			     	  + 1.0*(mi[i_p3*Ncc+j]-mi[i_m3*Ncc+j]))/(32.0*Dcell);
			m2[index] =(5.0*(mi[i*Ncc+j_p1]-mi[i*Ncc+j_m1])
			  		  + 4.0*(mi[i*Ncc+j_p2]-mi[i*Ncc+j_m2])
			  		  + 1.0*(mi[i*Ncc+j_p3]-mi[i*Ncc+j_m3]))/(32.0*Dcell);
		}

		if (dif_tag==3) {
			m1[index] = (58.0*(mi[i_p1*Ncc+j]-mi[i_m1*Ncc+j])
			            + 67.0*(mi[i_p2*Ncc+j]-mi[i_m2*Ncc+j])
			  	   	 	+ 22.0*(mi[i_p3*Ncc+j]-mi[i_m3*Ncc+j]))/(252.0*Dcell);
			m2[index] = (58.0*(mi[i*Ncc+j_p1]-mi[i*Ncc+j_m1])
					    + 67.0*(mi[i*Ncc+j_p2]-mi[i*Ncc+j_m2])
						- 22.0*(mi[i*Ncc+j_p3]-mi[i*Ncc+j_m3]))/(252.0*Dcell);
		}
	}
}
//--------------------------------------------------------------------
void lanczos_diff_2_tag(double *m1, double *m2, double *m11, double *m12, double *m21, double *m22, double Dcell, int Ncc, int dif_tag) {
    int i_m3,i_p3,j_m3,j_p3,i_m2,i_p2,j_m2,j_p2,i_m1,j_m1,i_p1,j_p1,i,j;
    int index;

    for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {

		//i_m1 = i-1;i_m2 = i-2;i_m3 = i-3;
		//i_p1 = i+1;i_p2 = i+2;i_p3 = i+3;
		//j_m1 = j-1;j_m2 = j-2;j_m3 = j-3;
		//j_p1 = j+1;j_p2 = j+2;j_p3 = j+3;

		//if (j_m3<0||j_p3>(Ncc-1)||i_m3<0||i_p3>(Ncc-1)) continue;
		//if (j_m2<0||j_p2>(Ncc-1)||i_m2<0||i_p2>(Ncc-1)) continue;
		//if (j_m1<0||j_p1>(Ncc-1)||i_m1<0||i_p1>(Ncc-1)) continue;

        if (i==0) {i_m1 = Ncc-1;i_m2 = Ncc-2;i_m3 = Ncc-3;}
        else if (i==1) {i_m1 = 0;i_m2 = Ncc-1;i_m3 = Ncc-2;}
        else if (i==2) {i_m1 = 1;i_m2 = 0;i_m3 = Ncc-1;}
        else {i_m1 = i-1;i_m2 = i-2;i_m3 = i-3;}
        if (j==0) {j_m1 = Ncc-1;j_m2 = Ncc-2;j_m3 = Ncc-3;}
        else if (j==1) {j_m1 = 0;j_m2 = Ncc-1;j_m3 = Ncc-2;}
        else if (j==2) {j_m1 = 1;j_m2 = 0;j_m3 = Ncc-1;}
        else {j_m1 = j-1;j_m2 = j-2;j_m3 = j-3;}
        if (i==Ncc-1) {i_p1 = 0;i_p2 = 1;i_p3 = 2;}
        else if (i==Ncc-2) {i_p1 = Ncc-1;i_p2 = 0;i_p3 = 1;}
        else if (i==Ncc-3) {i_p1 = Ncc-2;i_p2 = Ncc-1;i_p3 = 0;}
        else {i_p1 = i+1;i_p2 = i+2;i_p3 = i+3;}
        if (j==Ncc-1) {j_p1 = 0;j_p2 = 1;j_p3 = 2;}
        else if (j==Ncc-2) {j_p1 = Ncc-1;j_p2 = 0;j_p3 = 1;}
        else if (j==Ncc-2) {j_p1 = Ncc-2;j_p2 = Ncc-1;j_p3 = 0;}
        else {j_p1 = j+1;j_p2 = j+2;j_p3 = j+3;}

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
//--------------------------------------------------------------------
void write_2_signals(char *out1,char *out2,double *in1,double *in2, int Ncc) {
    long i,j;
    long index;
    FILE *f1,*f2;

    f1 =fopen(out1,"wb");
    f2 =fopen(out2,"wb");

    for (i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
        index = i*Ncc+j;

        fwrite(&in1[index],sizeof(double),1,f1);
        fwrite(&in2[index],sizeof(double),1,f2);
    }
    fclose(f1);
    fclose(f2);
}
//--------------------------------------------------------------------
void write_1_signal(char *out1,double *in1, int Ncc) {
    long i,j;
    long index;
    FILE *f1;

    f1 =fopen(out1,"wb");

    for (i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
        index = i*Ncc+j;
        fwrite(&in1[index],sizeof(double),1,f1);
    }
    fclose(f1);
}
//----------------------------------------------------------------------------------
void Loadin_grids_mesh(double boxsize, double xc1, double xc2, int Ncc, double *posx1, double *posx2) {
	double dsx = boxsize/(double)Ncc;

	int i,j;
	int index;
	for (i=0; i<Ncc; i++) for (j=0; j<Ncc; j++) {
		index = i*Ncc+j;
		posx1[index] = dsx*(double)(i)-boxsize*0.5+0.5*dsx+xc1;
		posx2[index] = dsx*(double)(j)-boxsize*0.5+0.5*dsx+xc2;
	}
}
//----------------------------------------------------------------------------------
void sdens_to_kappa(double p_mass_in, double* sdens_in, int Ncc, double Dcell, double zl, double zs, double *kappa) {

	int Nc2 = Ncc*2;
	double scrit;
	scrit = sigma_crit(zl,zs);//*(apr*apr/(Da(zl)*Da(zl)));
	//double kappa_max = 0.0;

	int i,j,index1,index2;
	for (i=Nc2/4;i<Nc2*3/4;i++) for (j=Nc2/4;j<Nc2*3/4;j++) {
		index1 = i*Nc2+j;
		index2 = (i-Nc2/4)*Nc2/2+(j-Nc2/4);
		kappa[index1] = p_mass_in*sdens_in[index2];///scrit;
	}
}
//----------------------------------------------------------------------------------
void sdens_to_kappac(double p_mass_in, double* sdens_in, int Ncc, double Dcell, double zl, double zs, double *kappac) {

	double scrit;
	scrit = sigma_crit(zl,zs);//*(apr*apr/(Da(zl)*Da(zl)));
	//double kappa_max = 0.0;

	int i,j,index;
	for (i=0;i<Ncc;i++) for (j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		kappac[index] = p_mass_in*sdens_in[index];///scrit;
	}
}
//--------------------------------------------------------------------
void kernel_green_iso(int Ncc, double *in, double Dcell) {
	int i,j;
	double x,y,r;
	double epsilon = Dcell*0.0000000001;
	double halfbox = Dcell*(double)Ncc/2.0;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if(i <=(Ncc/2)  && j <=(Ncc/2)) {

			x = (double)(i)*Dcell;
			y = (double)(j)*Dcell;
			r = sqrt(x*x+y*y+epsilon*epsilon);
			in[i*Ncc+j] = 1.0/M_PI*log(r);

		    //if(r > halfbox) {
			//	in[i*Ncc+j] = 0.0;
		    //}
		}
		else {
			if(i <= Ncc/2 && j > (Ncc/2)) {
				in[i*Ncc+j] = in[i*Ncc+(Ncc-j)];
			}
			if(i > (Ncc/2) && j <= (Ncc/2)) {
				in[i*Ncc+j] = in[(Ncc-i)*Ncc+j];
			}

			if(i > (Ncc/2) && j > (Ncc/2)) {
				in[i*Ncc+j] = in[(Ncc-i)*Ncc+(Ncc-j)];
			}
		}
	}
	//write_1_signal("out_iso.bin",in,Ncc);
}
//--------------------------------------------------------------------
void kernel_shears_iso(int Ncc,double *in1,double *in2,double Dcell) {
	int i,j;
	double x,y,r;
	double epsilon = Dcell*0.000000000000;
	double halfbox = Dcell*(double)Ncc/2.0;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if(i <(Ncc/2)  && j <(Ncc/2)) {
			x = (double)(i)*Dcell;
			y = (double)(j)*Dcell;
			r = sqrt(x*x+y*y+epsilon*epsilon);

			in1[i*Ncc+j] =  (y*y-x*x)/(M_PI*r*r*r*r);
			in2[i*Ncc+j] = (-2.0*x*y)/(M_PI*r*r*r*r);

		    if(r > halfbox) {
				in1[i*Ncc+j] = 0.0;
				in2[i*Ncc+j] = 0.0;
		    }
		}

		else {
			if(i < Ncc/2 && j >= (Ncc/2)) {
				in1[i*Ncc+j]  =  in1[i*Ncc+(Ncc-j)];
				in2[i*Ncc+j]  = -in2[i*Ncc+(Ncc-j)];
			}
			if(i >= (Ncc/2) && j < (Ncc/2)) {
				in1[i*Ncc+j]  =  in1[(Ncc-i)*Ncc+j];
				in2[i*Ncc+j]  = -in2[(Ncc-i)*Ncc+j];
			}

			if(i >= (Ncc/2) && j >= (Ncc/2)) {
				in1[i*Ncc+j]  =  in1[(Ncc-i)*Ncc+(Ncc-j)];
				in2[i*Ncc+j]  =  in2[(Ncc-i)*Ncc+(Ncc-j)];
			}
		}
	}
}
//--------------------------------------------------------------------
void kernel_alphas_iso(int Ncc,double *in1,double *in2,double Dcell) {
	int i,j;
	double x,y,r;
	double epsilon = Dcell*0.00000000001;
	double halfbox = Dcell*(double)Ncc/2.0;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if(i <=(Ncc/2)  && j <=(Ncc/2)) {
			x = (double)(i)*Dcell;
			y = (double)(j)*Dcell;
			r = sqrt(x*x+y*y+epsilon*epsilon);

			in1[i*Ncc+j] = x/(M_PI*r*r);
			in2[i*Ncc+j] = y/(M_PI*r*r);

		    if(r > halfbox) {
				in1[i*Ncc+j] = 0.0;
				in2[i*Ncc+j] = 0.0;
		    }
		}
		else {
			if(i <= Ncc/2 && j > (Ncc/2)) {
				in1[i*Ncc+j]  =  in1[i*Ncc+(Ncc-j)];
				in2[i*Ncc+j]  = -in2[i*Ncc+(Ncc-j)];
			}
			if(i > (Ncc/2) && j <= (Ncc/2)) {
				in1[i*Ncc+j]  = -in1[(Ncc-i)*Ncc+j];
				in2[i*Ncc+j]  =  in2[(Ncc-i)*Ncc+j];
			}

			if(i > (Ncc/2) && j > (Ncc/2)) {
				in1[i*Ncc+j]  = -in1[(Ncc-i)*Ncc+(Ncc-j)];
				in2[i*Ncc+j]  = -in2[(Ncc-i)*Ncc+(Ncc-j)];
			}
		}
	}
}
//--------------------------------------------------------------------
void kernel_smooth_iso(double sigma,int Ncc,double *in,double Dcell) {
	int i,j;
	double x,y,r;
	double epsilon = 0.00000001*Dcell;
	double cnorm = 0.0;
	double halfbox = Dcell*(double)Ncc/2.0;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if(i <(Ncc/2)  && j <(Ncc/2)) {
			x = (double)(i)*Dcell;
			y = (double)(j)*Dcell;
			r = sqrt(x*x+y*y+epsilon*epsilon);

			in[i*Ncc+j] = 1.0/(2.0*M_PI*sigma*sigma)*exp(-(r*r)/(2.0*sigma*sigma));

		    //if(r > halfbox) {
			//	in[i*Ncc+j] = 0.0;
		    //}
		}
		else {
			if(i < Ncc/2 && j >= (Ncc/2)) {
				in[i*Ncc+j] = in[i*Ncc+(Ncc-j)];
			}
			if(i >= (Ncc/2) && j < (Ncc/2)) {
				in[i*Ncc+j] = in[(Ncc-i)*Ncc+j];
			}

			if(i >= (Ncc/2) && j >= (Ncc/2)) {
				in[i*Ncc+j] = in[(Ncc-i)*Ncc+(Ncc-j)];
			}
		}
		cnorm += in[i*Ncc+j]*Dcell*Dcell;
	}

	double ctotal = 0.0;
	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		in[i*Ncc+j] = in[i*Ncc+j]/cnorm;
		ctotal += in[i*Ncc+j]*Dcell*Dcell;
	}
}
//--------------------------------------------------------------------
void fftw_r2c_2d(double *in_real, fftw_complex *in_fft, long Ncell, double Dcell) {
	long i,j;
	long nh = Ncell/2+1;
	long index;

	fftw_complex *in_fft_tmp = (fftw_complex *)fftw_malloc(Ncell*nh*sizeof(fftw_complex));
	fftw_plan plfd;
	plfd = fftw_plan_dft_r2c_2d(Ncell,Ncell,in_real,in_fft_tmp,FFTW_ESTIMATE);
	fftw_execute(plfd);

	for(i=0;i<Ncell;i++) for(j=0;j<nh;j++) {
		index = i*nh+j;
		in_fft[index][0] = in_fft_tmp[index][0];
		in_fft[index][1] = in_fft_tmp[index][1];
	}

	fftw_destroy_plan(plfd);
	fftw_free(in_fft_tmp);
}
//--------------------------------------------------------------------
void fftw_c2r_2d(fftw_complex *in_fft, double *in_real, long Ncell, double Dcell) {

	int i,j;
	//int nh = Ncell/2+1;
	int index;

	double *in_real_tmp = (double *)malloc(Ncell*Ncell*sizeof(double));
	fftw_plan plbd;
	plbd = fftw_plan_dft_c2r_2d(Ncell, Ncell, in_fft, in_real_tmp, FFTW_ESTIMATE);
	fftw_execute(plbd);

	for(i=0;i<Ncell;i++) for(j=0;j<Ncell;j++) {
		index = i*Ncell+j;
		in_real[index] = in_real_tmp[index];
	}

	fftw_destroy_plan(plbd);
	free(in_real_tmp);
}
//--------------------------------------------------------------------
void convolve_fft(double *in1, double *in2, double *out, long Ncell, double Dcell) {
	int i,j;
	int nh = Ncell/2+1;
	int index;
	double tmpr,tmpi;

	double in1_fftr = 0.0;
	double in1_ffti = 0.0;
	double in2_fftr = 0.0;
	double in2_ffti = 0.0;

	fftw_complex *in1_fft = (fftw_complex *)fftw_malloc(Ncell*nh*sizeof(fftw_complex));
	fftw_complex *in2_fft = (fftw_complex *)fftw_malloc(Ncell*nh*sizeof(fftw_complex));

	fftw_r2c_2d(in1, in1_fft, Ncell, Dcell);
	fftw_r2c_2d(in2, in2_fft, Ncell, Dcell);

	fftw_complex *out_fft = (fftw_complex *)fftw_malloc(Ncell*nh*sizeof(fftw_complex));
	for (i=0;i<Ncell;i++) for(j=0;j<nh;j++) {
		index = i*nh+j;
		in1_fftr = in1_fft[index][0];
		in1_ffti = in1_fft[index][1];

		in2_fftr = in2_fft[index][0];
		in2_ffti = in2_fft[index][1];
		tmpr = in1_fftr*in2_fftr-in1_ffti*in2_ffti;
		tmpi = in1_fftr*in2_ffti+in1_ffti*in2_fftr;
		out_fft[index][0] = tmpr;
		out_fft[index][1] = tmpi;
	}


	double *out_tmp = (double *)malloc(Ncell*Ncell*sizeof(double));
	fftw_c2r_2d(out_fft, out_tmp, Ncell, Dcell);

	for(i=0;i<Ncell;i++) for(j=0;j<Ncell;j++) {
		index = i*Ncell+j;
		out[index] = out_tmp[index]/(Ncell*Ncell)*Dcell*Dcell;
	}

	fftw_free(in1_fft);
	fftw_free(in2_fft);
	fftw_free(out_fft);

	free(out_tmp);
}
//--------------------------------------------------------------------
void kappa_to_phi(double *kappa, double *phi,int Ncc2,double Dcell) {

	double *green_iso = calloc(Ncc2*Ncc2,sizeof(double));

	kernel_green_iso(Ncc2,green_iso,Dcell);

	double *phi_tmp = malloc(Ncc2*Ncc2*sizeof(double));

	convolve_fft(kappa,green_iso,phi_tmp,Ncc2,Dcell);

	int i,j,index1,index2;
	for (i=Ncc2/4;i<Ncc2*3/4;i++) for (j=Ncc2/4;j<Ncc2*3/4;j++) {
		index1 = i*Ncc2+j;
		index2 = (i-Ncc2/4)*Ncc2/2+(j-Ncc2/4);
		phi[index2] = phi_tmp[index1];
	}

	free(phi_tmp);
}
//--------------------------------------------------------------------
void kappa_to_alphas(double *kappa,double *alpha1,double *alpha2,int Ncc2,double Dcell) {

	double *alpha1_iso = calloc(Ncc2*Ncc2,sizeof(double));
	double *alpha2_iso = calloc(Ncc2*Ncc2,sizeof(double));

	kernel_alphas_iso(Ncc2,alpha1_iso,alpha2_iso,Dcell);

	double *alpha1_tmp = malloc(Ncc2*Ncc2*sizeof(double));
	double *alpha2_tmp = malloc(Ncc2*Ncc2*sizeof(double));

	convolve_fft(kappa,alpha1_iso,alpha1_tmp,Ncc2,Dcell);
	convolve_fft(kappa,alpha2_iso,alpha2_tmp,Ncc2,Dcell);

	int i,j,index1,index2;
	for (i=Ncc2/4;i<Ncc2*3/4;i++) for (j=Ncc2/4;j<Ncc2*3/4;j++) {
		index1 = i*Ncc2+j;
		index2 = (i-Ncc2/4)*Ncc2/2+(j-Ncc2/4);
		alpha1[index2] = alpha1_tmp[index1];
		alpha2[index2] = alpha2_tmp[index1];
	}

	free(alpha1_tmp);
	free(alpha2_tmp);
}
//--------------------------------------------------------------------
void kappa_to_shears(double *kappa,double *shear1,double *shear2, int Ncc2,double Dcell) {

	double *shear1_iso = malloc(Ncc2*Ncc2*sizeof(double));
	double *shear2_iso = malloc(Ncc2*Ncc2*sizeof(double));

	kernel_shears_iso(Ncc2,shear1_iso,shear2_iso,Dcell);

	double *shear1_tmp = malloc(Ncc2*Ncc2*sizeof(double));
	double *shear2_tmp = malloc(Ncc2*Ncc2*sizeof(double));

	convolve_fft(kappa,shear1_iso,shear1_tmp,Ncc2,Dcell);
	convolve_fft(kappa,shear2_iso,shear2_tmp,Ncc2,Dcell);

	int i,j,index1,index2;
	for (i=Ncc2/4;i<Ncc2*3/4;i++) for (j=Ncc2/4;j<Ncc2*3/4;j++) {
		index1 = i*Ncc2+j;
		index2 = (i-Ncc2/4)*Ncc2/2+(j-Ncc2/4);
		shear1[index2] = shear1_tmp[index1];
		shear2[index2] = shear2_tmp[index1];
	}

	free(shear1_tmp);
	free(shear2_tmp);
	free(shear1_iso);
	free(shear2_iso);
}
//----------------------------------------------------------------------------------
void sdens_to_kappai(double p_mass_in, double* sdens_in, int Ncc, double Dcell, double zl, double zs, double *kappai) {

	int Ncc2 = Ncc*2;
	double * phi = (double *)malloc(Ncc*Ncc*sizeof(double));
	double * kappa = (double *)malloc(Ncc2*Ncc2*sizeof(double));
	sdens_to_kappa(p_mass_in,sdens_in,Ncc,Dcell,zl,zs,kappa);
	kappa_to_phi(kappa,phi,Ncc2,Dcell);
	free(kappa);

	double * phi1 = (double *)malloc(Ncc*Ncc*sizeof(double));
	double * phi2 = (double *)malloc(Ncc*Ncc*sizeof(double));
	lanczos_diff_1_tag(phi,phi1,phi2,Dcell,Ncc,-1);
	free(phi);

    double * phi11 = (double *)malloc(Ncc*Ncc*sizeof(double));
    double * phi12 = (double *)malloc(Ncc*Ncc*sizeof(double));
    double * phi21 = (double *)malloc(Ncc*Ncc*sizeof(double));
    double * phi22 = (double *)malloc(Ncc*Ncc*sizeof(double));

	lanczos_diff_2_tag(phi1,phi2,phi11,phi12,phi21,phi22,Dcell,Ncc,-1);
	free(phi1);
	free(phi2);

	int i,j,index;
	for (i=0;i<Ncc;i++) for (j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		kappai[index] = 0.5*(phi11[index]+phi22[index]);
	}

    free(phi11);
    free(phi12);
    free(phi21);
    free(phi22);
}
//--------------------------------------------------------------------
void calculate_mu(double *kappac,double *shear1,double *shear2, int Ncc, double *mu) {

	int i,j,index;
	for (i=0;i<Ncc;i++) for (j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		mu[index] = 1.0/((1.0-kappac[index])*(1.0-kappac[index])
				-shear1[index]*shear1[index]-shear2[index]*shear2[index]);
	}
}
//--------------------------------------------------------------------
void calculate_mu_mp(double *alpha1,double *alpha2,int Ncc,double Dcell,double *mu) {

    double * phi11 = (double *)malloc(Ncc*Ncc*sizeof(double));
    double * phi12 = (double *)malloc(Ncc*Ncc*sizeof(double));
    double * phi21 = (double *)malloc(Ncc*Ncc*sizeof(double));
    double * phi22 = (double *)malloc(Ncc*Ncc*sizeof(double));

	lanczos_diff_2_tag(alpha1,alpha2,phi11,phi12,phi21,phi22,Dcell,Ncc,1);

	int i,j,index;
	for (i=0;i<Ncc;i++) for (j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		mu[index] = 1.0/(1.0-(phi11[index]+phi22[index])+phi11[index]*phi22[index]-phi12[index]*phi21[index]);
	}

	free(phi11);
	free(phi12);
	free(phi21);
	free(phi22);
}
//--------------------------------------------------------------------
void calculate_td(double *phi,double *alpha1,double *alpha2, int Ncc, double *td) {
	double Kc = 1.0;
	//Kc = (1.0+zl)/c*(Dl*Ds/Dls)
	//td = Kc*(0.5*((al1)**2.0+(al2)**2.0)-phi)
	//td = Kc*(0.5*((x1-ys1)**2.0+(x2-ys2)**2.0)-phi)

	int i,j,index;
	for (i=0;i<Ncc;i++) for (j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		td[index] = Kc*((alpha1[index]*alpha1[index]+alpha2[index]*alpha2[index])-phi[index]);
	}
}
//----------------------------------------------------------------------------------
void sdens_to_wl(double p_mass_in, double* sdens_in, int Ncc, double Dcell, double zl, double zs, double *kappai, double *shear1, double *shear2) {

	int Ncc2 = Ncc*2;
	double * phi = (double *)malloc(Ncc*Ncc*sizeof(double));
	double * kappa = (double *)malloc(Ncc2*Ncc2*sizeof(double));
	sdens_to_kappa(p_mass_in,sdens_in,Ncc,Dcell,zl,zs,kappa);
	kappa_to_phi(kappa,phi,Ncc2,Dcell);
	free(kappa);

	double * phi1 = (double *)malloc(Ncc*Ncc*sizeof(double));
	double * phi2 = (double *)malloc(Ncc*Ncc*sizeof(double));
	lanczos_diff_1_tag(phi,phi1,phi2,Dcell,Ncc,0);
	free(phi);

    double * phi11 = (double *)malloc(Ncc*Ncc*sizeof(double));
    double * phi12 = (double *)malloc(Ncc*Ncc*sizeof(double));
    double * phi21 = (double *)malloc(Ncc*Ncc*sizeof(double));
    double * phi22 = (double *)malloc(Ncc*Ncc*sizeof(double));

	lanczos_diff_2_tag(phi1,phi2,phi11,phi12,phi21,phi22,Dcell,Ncc,0);
	free(phi1);
	free(phi2);

	int i,j,index;
	for (i=0;i<Ncc;i++) for (j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		kappai[index] = 0.5*(phi11[index]+phi22[index]);
		shear1[index] = 0.5*(phi11[index]-phi22[index]);
		shear2[index] = 0.5*(phi21[index]+phi12[index]);
	}

    free(phi11);
    free(phi12);
    free(phi21);
    free(phi22);
}
//----------------------------------------------------------------------------------
void sdens_to_sl(double p_mass_in, double* sdens_in, int Ncc, double Dcell, double zl, double zs, double *td,double *alpha1,double *alpha2,double * mu) {

	int Ncc2 = Ncc*2;
	double * phi = (double *)malloc(Ncc*Ncc*sizeof(double));
	double * kappa = (double *)malloc(Ncc2*Ncc2*sizeof(double));
	sdens_to_kappa(p_mass_in,sdens_in,Ncc,Dcell,zl,zs,kappa);
	kappa_to_phi(kappa,phi,Ncc2,Dcell);
	kappa_to_alphas(kappa,alpha1,alpha2,Ncc2,Dcell);
	free(kappa);

	//double * kappai = (double *)malloc(Ncc*Ncc*sizeof(double));
	//sdens_to_kappai(p_mass_in,sdens_in,Ncc,Dcell,zl,zs,kappai);
	//double * shear1 = (double *)malloc(Ncc*Ncc*sizeof(double));
	//double * shear2 = (double *)malloc(Ncc*Ncc*sizeof(double));
	//kappa_to_shears(kappa,shear1,shear2,Ncc2,Dcell);

	double * kappai = (double *)malloc(Ncc*Ncc*sizeof(double));
	double * shear1 = (double *)malloc(Ncc*Ncc*sizeof(double));
	double * shear2 = (double *)malloc(Ncc*Ncc*sizeof(double));
	sdens_to_wl(p_mass_in,sdens_in,Ncc,Dcell,zl,zs,kappai,shear1,shear2);
	calculate_mu(kappai,shear1,shear2,Ncc,mu);
	free(kappai);
	free(shear1);
	free(shear2);

	calculate_td(phi,alpha1,alpha2,Ncc,td);
	free(phi);
}
//----------------------------------------------------------------------------------
