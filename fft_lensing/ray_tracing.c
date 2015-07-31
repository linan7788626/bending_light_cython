#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <fitsio.h>
//--------------------------------------------------------------------
void Interplation_on_source_plane(double *source_map, double *posy1, double *posy2, double ysc1, double ysc2,double dsi, int nsx, int nsy, int nlx, int nly, double *lensed_map) {

	int i1,j1,i,j;
	int index;
	double xb1,xb2;
	double ww1,ww2,ww3,ww4,wx,wy;

	for(i=0;i<nlx;i++) for(j=0;j<nly;j++) {

		index = i*nlx+j;

		xb2 = (posy1[index]-ysc1)/dsi+(double)nsx/2.0-0.5;
		xb1 = (posy2[index]-ysc2)/dsi+(double)nsy/2.0-0.5;

		i1 = (int)xb1;
		j1 = (int)xb2;

		wx = 1.-(xb1-(double)(i1));
		wy = 1.-(xb2-(double)(j1));

		ww1 = wx*wy;
		ww2 = wx*(1.0-wy);
		ww3 = (1.0-wx)*wy;
		ww4 = (1.0-wx)*(1.0-wy);

		if (i1<0||i1>nsx-2||j1<0||j1>nsy-2) {
			ww1 = 0.0;ww2 = 0.0;ww3 = 0.0;ww4 = 0.0;
		}

		//if (i1<0) i1 = 0;
		//if (i1>nsx-2) i1 = nsx-2;
		//if (j1<0) j1 = 0;
		//if (j1>nsy-2) j1 = nsy-2;

		if (i1<0) i1 = nsx-2;
		if (i1>nsx-2) i1 = 0;
		if (j1<0) j1 = nsy-2;
		if (j1>nsy-2) j1 = 0;

		lensed_map[index] = ww1*source_map[i1*nsx+j1]
						  + ww2*source_map[i1*nsx+j1+1]
						  + ww3*source_map[(i1+1)*nsx+j1]
						  + ww4*source_map[(i1+1)*nsx+j1+1];
	}
}
//----------------------------------------------------------------------------------
void sfits_to_lfits(char *source_fits,double *posx1, double *posx2, double *alpha1, double *alpha2,double ysc1,double ysc2, double dsi, int nlx, int nly, char *lensed_fits) {

    fitsfile *fptr_in;
	int status,  nfound, anynull;
    long naxes[2],fpixel,npixels, i, j, index;

    float nullval;
    status = 0;

    fits_open_file(&fptr_in, source_fits, READONLY, &status);
    fits_read_keys_lng(fptr_in, "NAXIS", 1, 2, naxes, &nfound, &status);

	long nsx = naxes[0];
	long nsy = naxes[1];


    npixels  = nsx*nsy;
    fpixel   = 1;
    nullval  = 0;

	double *source_map = calloc(npixels,sizeof(double));
    fits_read_img(fptr_in, TDOUBLE, fpixel, npixels, &nullval,
                  source_map, &anynull, &status);

    fits_close_file(fptr_in, &status);

	double *posy1 = calloc(nlx*nly,sizeof(double));
	double *posy2 = calloc(nlx*nly,sizeof(double));

	for(i=0;i<nlx;i++) for(j=0;j<nly;j++){
		index = i*nlx+j;
		posy1[index] = posx1[index]-alpha1[index];
		posy2[index] = posx2[index]-alpha2[index];
	}

	double *lensed_map = calloc(nlx*nly,sizeof(double));
	Interplation_on_source_plane(source_map,posy1,posy2,ysc1,ysc2,dsi,nsx,nsy,nlx,nly,lensed_map);
    fitsfile *fptr_out;
    int fpixel_out = fpixel;
    long bitpix_out =  DOUBLE_IMG;
    long naxis_out = 2;
    long npixels_out = nlx*nly;
    long naxes_out[naxis_out];
    naxes_out[0] = nlx;
    naxes_out[1] = nly;

    status = 0;

    remove(lensed_fits);
    fits_create_file(&fptr_out, lensed_fits, &status);
    fits_create_img(fptr_out,  bitpix_out, naxis_out, naxes_out, &status);
    fits_write_img(fptr_out, TDOUBLE, fpixel_out, npixels_out, lensed_map, &status);
    fits_close_file(fptr_out, &status);

	free(posy1);
	free(posy2);
}
//----------------------------------------------------------------------------------
void smap_to_lmap(double *source_map, int nsx, int nsy, double *posx1, double *posx2, double *alpha1, double *alpha2,double ysc1,double ysc2, double dsi, int nlx, int nly, double *lensed_map) {

	int i,j,index;

	double *posy1 = calloc(nlx*nly,sizeof(double));
	double *posy2 = calloc(nlx*nly,sizeof(double));

	for(i=0;i<nlx;i++) for(j=0;j<nly;j++){
		index = i*nlx+j;
		posy1[index] = posx1[index]-alpha1[index];
		posy2[index] = posx2[index]-alpha2[index];
	}

	Interplation_on_source_plane(source_map,posy1,posy2,ysc1,ysc2,dsi,nsx,nsy,nlx,nly,lensed_map);

	free(posy1);
	free(posy2);
}
//----------------------------------------------------------------------------------
void read_from_fits(char *source_fits, long *nsx, long *nsy, double *source_map) {

    fitsfile *fptr_in;
	int status,  nfound, anynull;
    long naxes[2],fpixel,npixels;

    float nullval;
    status = 0;

    fits_open_file(&fptr_in, source_fits, READONLY, &status);
    fits_read_keys_lng(fptr_in, "NAXIS", 1, 2, naxes, &nfound, &status);

	*nsx = naxes[0];
	*nsy = naxes[1];

    npixels  = naxes[0]*naxes[1];
    fpixel   = 1;
    nullval  = 0;

	source_map = calloc(npixels,sizeof(double));
    fits_read_img(fptr_in, TDOUBLE, fpixel, npixels, &nullval,
                  source_map, &anynull, &status);

    fits_close_file(fptr_in, &status);
}
//----------------------------------------------------------------------------------
void save_to_fits(double *lensed_map, int nlx, int nly, char *lensed_fits) {

    fitsfile *fptr_out;
    int fpixel_out = 1;
    long bitpix_out =  DOUBLE_IMG;
    long naxis_out = 2;
    long npixels_out = nlx*nly;
    long naxes_out[naxis_out];
    naxes_out[0] = nlx;
    naxes_out[1] = nly;

    int status = 0;

    remove(lensed_fits);
    fits_create_file(&fptr_out, lensed_fits, &status);
    fits_create_img(fptr_out,  bitpix_out, naxis_out, naxes_out, &status);
    fits_write_img(fptr_out, TDOUBLE, fpixel_out, npixels_out, lensed_map, &status);
    fits_close_file(fptr_out, &status);
}
//----------------------------------------------------------------------------------
void simple_lensed_img(char *source_fits,double *posx1, double *posx2, double *alpha1, double *alpha2,double ysc1,double ysc2, double dsi, int nlx, int nly, char *lensed_fits) {

	long nsx, nsy;
	double *source_map = NULL;

	read_from_fits(source_fits, &nsx, &nsy, source_map);

	double *lensed_map = calloc(nlx*nly,sizeof(double));
	smap_to_lmap(source_map, nsx, nsy, posx1, posx2, alpha1, alpha2, ysc1, ysc2, dsi, nlx, nly, lensed_map);
	save_to_fits(lensed_map, nlx, nly, lensed_fits);

	free(source_map);
	free(lensed_map);
}
//----------------------------------------------------------------------------------
