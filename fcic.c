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
