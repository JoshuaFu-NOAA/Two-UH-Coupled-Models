#ifndef PREP_H
#define PREP_H

#include "pp.h"
#include "ini_codes.h"


int header_info (FILE * file, HEADER * info);
int read_header_records (FILE * file, HEADER * current, int nrecl1);

void extract_block_1 (UINT8 * ddr1, TIMES * all, LABELS * run);
void extract_block_3 (UINT8 * ddr3, DIMENSIONS * model,
		      GRID * vert, char *def);
void extract_block_4 (unsigned long long *ddr4, DIMENSIONS * model,
		      char *def);


double *read_sp (FILE * file, int ilen, HEADER * info, size_t * rb);
double *read_gl (FILE * file, int ilen, HEADER * info, DIMENSIONS * model,
		 size_t * rb);
double *read_g3 (int no, long data_offset, FILE * file, int ilen,
		 HEADER * info, DIMENSIONS * model, size_t * rb);
double *read_sgl (int no, long data_offset, FILE * file, HEADER * info,
		  DIMENSIONS * model, size_t * rb);
double *read_f1 (int no, long data_offset, FILE * file, DIMENSIONS * model,
		 HEADER * info, size_t * rb);
double *read_f3 (int no, long data_offset, FILE * file, DIMENSIONS * model,
		 HEADER * info, size_t * rb);
double *read_f4 (int no, long data_offset, FILE * file, DIMENSIONS * model,
		 HEADER * info, size_t * rb);
double *read_g1r (long data_offset, FILE * file, int ilenpp, int ilenr,
		  int start, HEADER * info, DIMENSIONS * model, size_t * rb);
double *read_g2r (int no, long data_offset, FILE * file, int ilen,
		  HEADER * info, DIMENSIONS * model, size_t * rb);
double *read_g3r (long data_offset, FILE * file, int ilenpp, int ilenr,
		  int start, HEADER * info, DIMENSIONS * model, size_t * rb, int convtoieee);

void compare (DIMENSIONS * d1, DIMENSIONS * d2, TIMES * t1, TIMES * t2,
	      GRID * g1, GRID * g2);


void make_netcdf_31 (char *argv[], netCDF_file * restart, HEADER * info,
		     DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		     char *ofile, double *fsvo, double *favo, double *fsu,
		     double *fau, double *fsv, double *fav, double *fsd,
		     double *fad, double *fstp, double *fatp, double *fstpm,
		     double *fatpm, double *fsu0, double *fau0, double *fsdu0,
		     double *fadu0);

void make_netcdf_32 (char *argv[], netCDF_file * restart, HEADER * info,
		     DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		     char *ofile,
		     double *q, double *x, double *xt,
		     double *lammp, double *phimp, double *sigmp);

/* void make_netcdf_33 (char *argv[], netCDF_file *restart, HEADER *info, 
   DIMENSIONS *d, GRID *g, LABELS *run, TIMES *t, 
   char *ofile,
   double *fsd, double *fad); 
   
   void make_netcdf_34 (char *argv[], netCDF_file *restart, HEADER *info, 
   DIMENSIONS *d, GRID *g, LABELS *run, TIMES *t, 
   char *ofile, 
   double *fstp, double *fatp,
   double *fstpm, double *fatpm, 
   double *fsu0, double *fau0,
   double *fsdu0, double *fadu0);
   */

void make_netcdf_35 (char *argv[], netCDF_file * restart, HEADER * info,
		     DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		     char *ofile, double *vom1, double *dm1, double *tm1,
		     double *alpsm1, double *qm1, double *xm1, double *xtm1);

void make_netcdf_36 (char *argv[], netCDF_file * restart, HEADER * info,
		     DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		     char *ofile, double *um1, double *vm1);

void make_netcdf_37 (char *argv[], netCDF_file * restart, HEADER * info,
		     DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		     char *ofile, double *geospm, double *tsm, double *wsm,
		     double *wlm, double *snm, double *slmm, double *az0m,
		     double *albm, double *varpm, double *varorm,
		     double *forestm, double *vgratm, double *vltm,
		     double *wsmxm, double *faom, double *tdm, double *tdclm,
		     double *apsm, double *aprlm, double *aprcm,
		     double *aprsm, double *ustrgwm, double *vstrgwm,
		     double *vdisgwm, double *aclcovm, double *temp2m,
		     double *dew2m, double *wind10m, double *u10m,
		     double *v10m, double *sradsm, double *tradsm,
		     double *srad0m, double *trad0m, double *vdism,
		     double *ustrm, double *vstrm, double *ahfsm,
		     double *evapm, double *ahflm, double *wlm1m,
		     double *tsm1m, double *tdm1m, double *wsm1m,
		     double *wdm1m, double *snm1m, double *emterm,
		     double *trsolm, double *runoffm, double *srad0um,
		     double *sradsum, double *tradsum, double *albedom,
		     double *tsurfm, double *tsnm, double *td3m, double *td4m,
		     double *td5m, double *tsnm1m, double *td3m1m,
		     double *td4m1m, double *td5m1m, double *tdclm1m,
		     double *seaicem, double *sicedm, double *ustar3m,
		     double *teffm, double *glacm, double *aclcm,
		     double *aclcacm, double *snmelm, double *runtocm,
		     double *tslinm, double *dsnacm, double *t2maxm,
		     double *t2minm, double *tsmaxm, double *tsminm,
		     double *wimaxm, double *topmaxm, double *aclcvm,
		     double *qvim, double *alwcvim, double *runlndm,
		     double *rgcgnm, double *sodifm, double *srafsm,
		     double *trafsm, double *sraf0m, double *traf0m,
		     double *emtefm, double *trsofm, double *tclfsm,
		     double *sclfsm, double *tclf0m, double *sclf0m,
		     double *tkem, double *tkem1m, double *auxil1m,
		     double *auxil2m);



void write_anc (netCDF_file * initial, HEADER * info,
		DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		char *ofile, double *svo, double *sd, double *stp, double *q,
		int format);

void write_sncold (netCDF_file * initial, HEADER * info,
		   DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		   char *ofile, double *geosp, double *ts, double *ws,
		   double *wl, double *sn, double *slm, double *az0,
		   double *alb, double *varp, double *varor, double *forest,
		   double *vgrat, double *vlt, double *wsmx, double *fao);
void write_snc (GRID * g, char *tfile, char *ofile, SINITIAL * vars,
		int format, DIMENSIONS * d, LABELS * run, TIMES * t);
void read_snc (GRID * g, TIMES * t, char *ifile, SINITIAL * vars);

void write_ync (GRID * grid, char *ofile, double *year, int iyear, int icode, int fromcray);

double *read_ync (GRID * grid, char *ifilename, char *date_unit, int *icode);
double *read_onc (GRID * grid, char *ifilename, char *date_unit, int *icode);
void write_onc (char *argv[], GRID * grid, char *ofile, double *ozon,
		int ocode, char *dsource);

double *read_year3_cray (FILE * file, GRID * grid, size_t * rb);
double *read_year_cray (FILE * file, GRID * grid, size_t * rb);
double *read_ozon_cray (FILE * file, GRID * grid, size_t * rb);
double *read_fld_cray (FILE * file, GRID * grid, size_t * rb);
double *read_fnc (GRID * grid, char *ifilename, int *icode);
void write_fnc (GRID * grid, char *ofile, double *fld, int icode);
#endif
