void Interplation_on_source_plane(double *source_map, double *posy1, double *posy2, double ysc1, double ysc2,double dsi, int nsx, int nsy, int nlx, int nly, double *lensed_map);
void sfits_to_lfits(char *source_fits,double *posx1, double *posx2, double *alpha1, double *alpha2,double ysc1,double ysc2, double dsi, int nlx, int nly, char *lensed_fits);
void smap_to_lmap(double *source_map, int nsx, int nsy, double *posx1, double *posx2, double *alpha1, double *alpha2,double ysc1,double ysc2, double dsi, int nlx, int nly, double *lensed_map);
void read_from_fits(char *source_fits, long *nsx, long *nsy, double *source_map);
void save_to_fits(double *lensed_map, int nlx, int nly, char *lensed_fits);
void simple_lensed_img(char *source_fits,double *posx1, double *posx2, double *alpha1, double *alpha2,double ysc1,double ysc2, double dsi, int nlx, int nly, char *lensed_fits);
