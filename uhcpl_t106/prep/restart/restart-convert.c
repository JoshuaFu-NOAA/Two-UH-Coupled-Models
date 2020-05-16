#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "prep.h"

int read_special_header (FILE * file, HEADER * info, DIMENSIONS * model,
			 int nrecl1);

int
main (int argc, char *argv[])
{
  HEADER header;
  HEADER header31;
  DIMENSIONS model;
  TIMES all;
  LABELS run;
  GRID grid;
  int nlon, nlon31=0, ib, convtoieee;

  int ilenb, start, ilen;

/* global file pointer for read data */

  FILE *file[7];

/* global accessible netCDF file information */

  netCDF_file restart[7];

  double *fsvo=0, *favo=0;
  double *fsu=0, *fau=0;
  double *fsv=0, *fav=0;

  double *q, *x, *xt=0, *lammp, *phimp, *sigmp;

  double *fsd=0, *fad=0;

  double *fstp, *fatp;
  double *fstpm, *fatpm;
  double *fsu0, *fau0;
  double *fsdu0, *fadu0;

  double *vom1, *dm1, *tm1, *alpsm1, *qm1, *xm1, *xtm1=0;

  double *um1, *vm1;

  double *geospm, *tsm, *wsm, *wlm, *snm, *slmm, *az0m, *albm, *varpm;
  double *varorm, *forestm, *vgratm, *vltm, *wsmxm, *faom, *tdm, *tdclm;
  double *apsm, *aprlm, *aprcm, *aprsm, *ustrgwm, *vstrgwm, *vdisgwm;
  double *aclcovm, *temp2m, *dew2m, *wind10m, *u10m, *v10m, *sradsm;
  double *tradsm, *srad0m, *trad0m, *vdism, *ustrm, *vstrm, *ahfsm;
  double *evapm, *ahflm, *wlm1m, *tsm1m, *tdm1m, *wsm1m, *wdm1m, *snm1m;
  double *emterm, *trsolm, *runoffm, *srad0um, *sradsum, *tradsum;
  double *albedom, *tsurfm, *tsnm, *td3m, *td4m, *td5m, *tsnm1m;
  double *td3m1m, *td4m1m, *td5m1m, *tdclm1m, *seaicem, *sicedm;
  double *ustar3m, *teffm, *glacm, *aclcm, *aclcacm, *snmelm, *runtocm;
  double *tslinm, *dsnacm, *t2maxm, *t2minm, *tsmaxm, *tsminm, *wimaxm;
  double *topmaxm, *aclcvm, *qvim, *alwcvim, *runlndm, *rgcgnm, *sodifm;
  double *srafsm, *trafsm, *sraf0m, *traf0m, *emtefm, *trsofm, *tclfsm;
  double *sclfsm, *tclf0m, *sclf0m, *tkem, *tkem1m, *auxil1m, *auxil2m;

  size_t bfsvo, bfavo;
  size_t bfsu, bfau;
  size_t bfsv, bfav;

  size_t bq, bx, bxt, blammp, bphimp, bsigmp;

  size_t bfsd, bfad;

  size_t bfstp, bfatp;
  size_t bfstpm, bfatpm;
  size_t bfsu0, bfau0;
  size_t bfsdu0, bfadu0;

  size_t bdm1, btm1, balpsm1, bqm1, bxm1, bxtm1;

  size_t bum1, bvm1;

  size_t bgeospm, btsm, bwsm, bwlm, bsnm, bslmm, baz0m, balbm, bvarpm;
  size_t bvarorm, bforestm, bvgratm, bvltm, bwsmxm, bfaom, btdm, btdclm;
  size_t bapsm, baprlm, baprcm, baprsm, bustrgwm, bvstrgwm, bvdisgwm;
  size_t baclcovm, btemp2m, bdew2m, bwind10m, bu10m, bv10m, bsradsm;
  size_t btradsm, bsrad0m, btrad0m, bvdism, bustrm, bvstrm, bahfsm;
  size_t bevapm, bahflm, bwlm1m, btsm1m, btdm1m, bwsm1m, bwdm1m, bsnm1m;
  size_t bemterm, btrsolm, brunoffm, bsrad0um, bsradsum, btradsum;
  size_t balbedom, btsurfm, btsnm, btd3m, btd4m, btd5m, btsnm1m;
  size_t btd3m1m, btd4m1m, btd5m1m, btdclm1m, bseaicem, bsicedm;
  size_t bustar3m, bteffm, bglacm, baclcm, baclcacm, bsnmelm, bruntocm;
  size_t btslinm, bdsnacm, bt2maxm, bt2minm, btsmaxm, btsminm, bwimaxm;
  size_t btopmaxm, baclcvm, bqvim, balwcvim, brunlndm, brgcgnm, bsodifm;
  size_t bsrafsm, btrafsm, bsraf0m, btraf0m, bemtefm, btrsofm, btclfsm;
  size_t bsclfsm, btclf0m, bsclf0m, btkem, btkem1m, bauxil1m, bauxil2m;

  long data_offset;
  int i, nrecl, nrecl1, nrecs;

  char *iname = NULL, *oname =
    NULL, ifile[7][FILENAME_MAX], ofile[7][FILENAME_MAX];


  if (argc <= 2)
    {
      fprintf (stderr, "Usage: %s  restart-basename   Netcdf-basename \n",
	       argv[0]);
      fprintf (stderr, "convert ECHAM4 Restart data from CRAY to Netcdf \n");
      exit (1);
    }
  /* check output filename */
  if (ofile == NULL || ifile == NULL)
    {
      fprintf (stderr, "Usage: %s  restart-basename   Netcdf-basename \n",
	       argv[0]);
      fprintf (stderr, " don't forget the filename \n\a");
      exit (1);
    }
  iname = argv[1];
  oname = argv[2];

  fprintf (stdout, "  < %s.3[1-7] \n", iname);
  fprintf (stdout, "  > %s.3[1,2,5,6,7].nc \n\n", oname);

  /* create filenames */

  for (i = 0; i < 7; i++)
    {
      sprintf (ofile[i], "%s.3%d.nc", oname, i + 1);
      sprintf (ifile[i], "%s.3%d", iname, i + 1);
    }

  for (i = 0; i < 7; i++)
    {
      if ((file[i] = fopen (ifile[i], "r")) == NULL)
	{
	  fprintf (stderr, " could not open file: %s  \n", ifile[i]);
	  exit (1);
	}
    }

  /* the files grouped together should be the right one to use togehter */

  for (i = 0; i < 7; i++)
    {

      printf ("\nExtract information from %s:\n\n", ifile[i]);

      /* check type of file and total header length */

      rewind (file[i]);
      nrecl1 = header_info (file[i], &header);

      if (i == 0)
	{

	  /* extract header records and print information */

	  nrecs = read_header_records (file[i], &header, nrecl1);

	  nrecs = (int) header.ddr1[14];
	  printf ("\nNumber of records     : %d\n", nrecs);

	  /* extract information from record 1 */

	  extract_block_1 (header.ddr1, &all, &run);

	  /* extract information from record 3 */

	  nrecl = header.ddr2[2];
	  nlon31 = header.ddr3[header.ddr3[14]];
	  header31.ddr3 = (unsigned long long *) malloc (nrecl * sizeof (unsigned long long));
	  for (ib = 0; ib < nrecl; ib++)
	    header31.ddr3[ib] = header.ddr3[ib];

	  extract_block_3 (header.ddr3, &model, &grid, header.def);

	  /* extract information from record 4 */

	  extract_block_4 (header.ddr4, &model, header.def);

	  data_offset = ftell (file[i]);

	  {
	    off_t file_bytes;
	    long header_bytes, data_words;
	    struct stat ifile_info;

	    nrecl1 = header_info (file[4], &header);
	    nrecs = read_header_records (file[4], &header, nrecl1);
	    header_bytes = ftell (file[4]);

	    stat (ifile[4], &ifile_info);
	    file_bytes = ifile_info.st_size;

	    data_words = (file_bytes - header_bytes) / 8;
	    model.ntrac = (int) (((data_words / model.ngl)
				  - 5 * model.nlp2 * model.nlev - model.nlp2 -
				  1) / (model.nlon * model.nlev));
	  }

	}
      else if (i == 1)
	{			/* Special case due to non ddr header - gl buffer */

	  printf
	    ("\n********************************************************\n");

	  printf ("\nNon ddr header in %s.\n\n", ifile[i]);

	  nrecs = read_special_header (file[i], &header, &model, nrecl1);

	  data_offset = ftell (file[i]);

	  q = read_sgl (1, data_offset, file[i], &header, &model, &bq);
	  x = read_sgl (2, data_offset, file[i], &header, &model, &bx);
	  if (model.nhtrac > 0)
	    {
	      bxt = 0;
	      xt = read_sgl (3, data_offset, file[i], &header, &model, &bxt);
	    }
	  lammp =
	    read_sgl (4, data_offset, file[i], &header, &model, &blammp);
	  phimp =
	    read_sgl (5, data_offset, file[i], &header, &model, &bphimp);
	  sigmp =
	    read_sgl (6, data_offset, file[i], &header, &model, &bsigmp);

	  make_netcdf_32 (argv, &restart[i], &header, &model, &grid,
			  &run, &all, ofile[i], q, x, xt, lammp, phimp,
			  sigmp);

	  free (q);
	  free (x);
	  free (xt);
	  free (lammp);
	  free (phimp);
	  free (sigmp);

	  printf
	    ("\n********************************************************\n");

	}
      else
	{

	  /* extract header records and print information */

	  nrecs = read_header_records (file[i], &header, nrecl1);

	  nrecs = (int) header.ddr1[14];
	  printf ("\nNumber of records     : %d\n", nrecs);

	  /* extract information from record 1 */

	  extract_block_1 (header.ddr1, &all, &run);

	  /* extract information from record 3 */

	  nlon = header.ddr3[header.ddr3[14]];

	  if (nlon != nlon31)
	    {
	      fprintf (stderr, "warning: gridsize different from unit:31\n");
	      fprintf (stderr, "         nlon = %d nlon31 = %d\n", nlon, nlon31);
	      fprintf (stderr, "         using gridsize from unit:31 !!!!\n");
	      extract_block_3 (header31.ddr3, &model, &grid, header.def);
	    }
	  else
	      extract_block_3 (header.ddr3, &model, &grid, header.def);

	  /* extract information from record 4 */

	  extract_block_4 (header.ddr4, &model, header.def);

	  data_offset = ftell (file[i]);
	}

      switch (i)
	{			/* Handle all cases separatly ... */
	  /* nmp1 = nm + 1 */
	case 0:		/* f1 buffer */

	  fsvo = read_f1 (1, data_offset, file[i], &model, &header, &bfsvo);
	  favo = read_f1 (2, data_offset, file[i], &model, &header, &bfavo);
	  fsu = read_f1 (3, data_offset, file[i], &model, &header, &bfsu);
	  fau = read_f1 (4, data_offset, file[i], &model, &header, &bfau);
	  fsv = read_f1 (5, data_offset, file[i], &model, &header, &bfsv);
	  fav = read_f1 (6, data_offset, file[i], &model, &header, &bfav);

	  break;
	case 1:		/* handled above */
	  break;
	case 2:		/* f3 buffer */

	  fsd = read_f3 (1, data_offset, file[i], &model, &header, &bfsd);
	  fad = read_f3 (2, data_offset, file[i], &model, &header, &bfad);

	  break;
	case 3:		/* f4 buffer */

	  fstp = read_f4 (1, data_offset, file[i], &model, &header, &bfstp);	/* nlevp1,2,nmp1 */
	  fatp = read_f4 (2, data_offset, file[i], &model, &header, &bfatp);	/* nlevp1,2,nmp1 */
	  fstpm = read_f4 (3, data_offset, file[i], &model, &header, &bfstpm);	/* nlevp1,2,nmp1 */
	  fatpm = read_f4 (4, data_offset, file[i], &model, &header, &bfatpm);	/* nlevp1,2,nmp1 */
	  fsu0 = read_f4 (5, data_offset, file[i], &model, &header, &bfsu0);	/* nlev */
	  fau0 = read_f4 (6, data_offset, file[i], &model, &header, &bfau0);	/* nlev */
	  fsdu0 = read_f4 (7, data_offset, file[i], &model, &header, &bfsdu0);	/* nlev */
	  fadu0 = read_f4 (8, data_offset, file[i], &model, &header, &bfadu0);	/* nlev */


	  /* write data from file 31, 33, 34 to file 31  */
	  make_netcdf_31 (argv, &restart[0], &header, &model, &grid,
			  &run, &all, ofile[0], fsvo, favo, fsu, fau, fsv,
			  fav, fsd, fad, fstp, fatp, fstpm, fatpm, fsu0, fau0,
			  fsdu0, fadu0);

	  /* file 31 */
	  free (fsvo);
	  free (favo);
	  free (fsu);
	  free (fau);
	  free (fsv);
	  free (fav);
	  /* file 31 */
	  free (fsd);
	  free (fad);
	  /* file 34 */
	  free (fstp);
	  free (fatp);
	  free (fstpm);
	  free (fatpm);
	  free (fsu0);
	  free (fau0);
	  free (fsdu0);
	  free (fadu0);

	  break;
	case 4:		/* g1 buffer */

	  /* Due to limited information in file header ntrac has to be determined by the file size */

	  ilenb = (5 * model.nlp2 * model.nlev + model.nlp2 +
		   model.nlon * model.nlev * model.ntrac + 1);

	  start = 0;
	  ilen = model.nlp2 * model.nlev;
	  vom1 =
	    read_g1r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bvm1);

	  start += model.nlp2 * model.nlev;
	  ilen = model.nlp2 * model.nlev;
	  dm1 =
	    read_g1r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bdm1);

	  start += model.nlp2 * model.nlev;
	  ilen = model.nlp2 * model.nlev;
	  tm1 =
	    read_g1r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btm1);

	  start += model.nlp2 * model.nlev;
	  ilen = model.nlp2;
	  alpsm1 =
	    read_g1r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &balpsm1);

	  start += model.nlp2;
	  ilen = model.nlp2 * model.nlev;
	  qm1 =
	    read_g1r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bqm1);

	  start += model.nlp2 * model.nlev;
	  ilen = model.nlp2 * model.nlev;
	  xm1 =
	    read_g1r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bxm1);

	  if (model.ntrac > 0)
	    {
	      start += model.nlp2 * model.nlev;
	      ilen = model.nlon * model.nlev * model.ntrac;
	      xtm1 =
		read_g1r (data_offset, file[i], ilenb, ilen, start, &header,
			  &model, &bxtm1);
	    }

	  make_netcdf_35 (argv, &restart[i], &header, &model, &grid,
			  &run, &all, ofile[i], vom1, dm1, tm1, alpsm1,
			  qm1, xm1, xtm1);

	  free (vom1);
	  free (dm1);
	  free (tm1);
	  free (alpsm1);
	  free (qm1);
	  free (xm1);
	  free (xtm1);

	  break;
	case 5:		/* g2 buffer */
	  ilen = model.nlp2 * model.nlev * model.ngl;

	  um1 =
	    read_g2r (1, data_offset, file[i], ilen, &header, &model, &bum1);
	  vm1 =
	    read_g2r (2, data_offset, file[i], ilen, &header, &model, &bvm1);

	  make_netcdf_36 (argv, &restart[i], &header, &model, &grid,
			  &run, &all, ofile[i], um1, vm1);
	  free (um1);
	  free (vm1);

	  break;
	case 6:		/* g3 buffer */

	  ilenb =
	    (4 * model.nlp2 * model.nlev + 96 * model.nlp2 +
	     2 * model.nlp2 * model.nlevp1);

	  convtoieee = 1;
	  start = 0;
	  ilen = model.nlp2;
	  geospm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bgeospm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  wsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bwsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  wlm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bwlm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  snm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsnm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  slmm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bslmm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  az0m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &baz0m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  albm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &balbm, convtoieee);

	  convtoieee = 0;
	  start += model.nlp2;
	  ilen = model.nlp2;
	  varpm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bvarpm, convtoieee);

	  convtoieee = 1;
	  start += model.nlp2;
	  ilen = model.nlp2;
	  varorm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bvarorm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  forestm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bforestm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  vgratm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bvgratm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  vltm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bvltm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  wsmxm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bwsmxm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  faom =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bfaom, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tdm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btdm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tdclm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btdclm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  apsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bapsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  aprlm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &baprlm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  aprcm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &baprcm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  aprsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &baprsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  ustrgwm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bustrgwm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  vstrgwm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bvstrgwm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  vdisgwm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bvdisgwm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  aclcovm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &baclcovm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  temp2m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btemp2m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  dew2m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bdew2m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  wind10m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bwind10m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  u10m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bu10m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  v10m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bv10m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  sradsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsradsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tradsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btradsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  srad0m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsrad0m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  trad0m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btrad0m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  vdism =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bvdism, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  ustrm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bustrm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  vstrm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bvstrm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  ahfsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bahfsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  evapm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bevapm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  ahflm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bahflm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  wlm1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bwlm1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tsm1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btsm1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tdm1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btdm1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  wsm1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bwsm1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  wdm1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bwdm1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  snm1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsnm1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2 * model.nlevp1;
	  emterm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bemterm, convtoieee);

	  start += model.nlp2 * model.nlevp1;
	  ilen = model.nlp2 * model.nlevp1;
	  trsolm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btrsolm, convtoieee);

	  start += model.nlp2 * model.nlevp1;
	  ilen = model.nlp2;
	  runoffm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &brunoffm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  srad0um =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsrad0um, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  sradsum =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsradsum, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tradsum =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btradsum, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  albedom =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &balbedom, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tsurfm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btsurfm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tsnm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btsnm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  td3m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btd3m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  td4m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btd4m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  td5m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btd5m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tsnm1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btsnm1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  td3m1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btd3m1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  td4m1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btd4m1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  td5m1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btd5m1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tdclm1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btdclm1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  seaicem =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bseaicem, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  sicedm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsicedm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  ustar3m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bustar3m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  teffm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bteffm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  glacm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bglacm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2 * model.nlev;
	  aclcm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &baclcm, convtoieee);

	  start += model.nlp2 * model.nlev;
	  ilen = model.nlp2 * model.nlev;
	  aclcacm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &baclcacm, convtoieee);

	  start += model.nlp2 * model.nlev;
	  ilen = model.nlp2;
	  snmelm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsnmelm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  runtocm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bruntocm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tslinm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btslinm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  dsnacm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bdsnacm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  t2maxm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bt2maxm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  t2minm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bt2minm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tsmaxm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btsmaxm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tsminm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btsminm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  wimaxm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bwimaxm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  topmaxm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btopmaxm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  aclcvm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &baclcvm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  qvim =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bqvim, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  alwcvim =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &balwcvim, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  runlndm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &brunlndm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  rgcgnm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &brgcgnm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  sodifm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsodifm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  srafsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsrafsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  trafsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btrafsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  sraf0m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsraf0m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  traf0m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btraf0m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2 * 2;
	  emtefm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bemtefm, convtoieee);

	  start += model.nlp2 * 2;
	  ilen = model.nlp2 * 2;
	  trsofm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btrsofm, convtoieee);

	  start += model.nlp2 * 2;
	  ilen = model.nlp2;
	  tclfsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btclfsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  sclfsm =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsclfsm, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  tclf0m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btclf0m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  sclf0m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bsclf0m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2 * model.nlev;
	  tkem =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btkem, convtoieee);

	  start += model.nlp2 * model.nlev;
	  ilen = model.nlp2 * model.nlev;
	  tkem1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &btkem1m, convtoieee);

	  start += model.nlp2 * model.nlev;
	  ilen = model.nlp2;
	  auxil1m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bauxil1m, convtoieee);

	  start += model.nlp2;
	  ilen = model.nlp2;
	  auxil2m =
	    read_g3r (data_offset, file[i], ilenb, ilen, start, &header,
		      &model, &bauxil2m, convtoieee);

	  make_netcdf_37 (argv, &restart[i], &header, &model, &grid,
			  &run, &all, ofile[i],
			  geospm, tsm, wsm, wlm, snm, slmm, az0m, albm, varpm,
			  varorm, forestm, vgratm, vltm, wsmxm, faom, tdm,
			  tdclm, apsm, aprlm, aprcm, aprsm, ustrgwm, vstrgwm,
			  vdisgwm, aclcovm, temp2m, dew2m, wind10m, u10m,
			  v10m, sradsm, tradsm, srad0m, trad0m, vdism, ustrm,
			  vstrm, ahfsm, evapm, ahflm, wlm1m, tsm1m, tdm1m,
			  wsm1m, wdm1m, snm1m, emterm, trsolm, runoffm,
			  srad0um, sradsum, tradsum, albedom, tsurfm, tsnm,
			  td3m, td4m, td5m, tsnm1m, td3m1m, td4m1m, td5m1m,
			  tdclm1m, seaicem, sicedm, ustar3m, teffm, glacm,
			  aclcm, aclcacm, snmelm, runtocm, tslinm, dsnacm,
			  t2maxm, t2minm, tsmaxm, tsminm, wimaxm, topmaxm,
			  aclcvm, qvim, alwcvim, runlndm, rgcgnm, sodifm,
			  srafsm, trafsm, sraf0m, traf0m, emtefm, trsofm,
			  tclfsm, sclfsm, tclf0m, sclf0m, tkem, tkem1m,
			  auxil1m, auxil2m);

	  free (geospm);
	  free (tsm);
	  free (wsm);
	  free (wlm);
	  free (snm);
	  free (slmm);
	  free (az0m);
	  free (albm);
	  free (varpm);
	  free (varorm);
	  free (forestm);
	  free (vgratm);
	  free (vltm);
	  free (wsmxm);
	  free (faom);
	  free (tdm);
	  free (tdclm);
	  free (apsm);
	  free (aprlm);
	  free (aprcm);
	  free (aprsm);
	  free (ustrgwm);
	  free (vstrgwm);
	  free (vdisgwm);
	  free (aclcovm);
	  free (temp2m);
	  free (dew2m);
	  free (wind10m);
	  free (u10m);
	  free (v10m);
	  free (sradsm);
	  free (tradsm);
	  free (srad0m);
	  free (trad0m);
	  free (vdism);
	  free (ustrm);
	  free (vstrm);
	  free (ahfsm);
	  free (evapm);
	  free (ahflm);
	  free (wlm1m);
	  free (tsm1m);
	  free (tdm1m);
	  free (wsm1m);
	  free (wdm1m);
	  free (snm1m);
	  free (emterm);
	  free (trsolm);
	  free (runoffm);
	  free (srad0um);
	  free (sradsum);
	  free (tradsum);
	  free (albedom);
	  free (tsurfm);
	  free (tsnm);
	  free (td3m);
	  free (td4m);
	  free (td5m);
	  free (tsnm1m);
	  free (td3m1m);
	  free (td4m1m);
	  free (td5m1m);
	  free (tdclm1m);
	  free (seaicem);
	  free (sicedm);
	  free (ustar3m);
	  free (teffm);
	  free (glacm);
	  free (aclcm);
	  free (aclcacm);
	  free (snmelm);
	  free (runtocm);
	  free (tslinm);
	  free (dsnacm);
	  free (t2maxm);
	  free (t2minm);
	  free (tsmaxm);
	  free (tsminm);
	  free (wimaxm);
	  free (topmaxm);
	  free (aclcvm);
	  free (qvim);
	  free (alwcvim);
	  free (runlndm);
	  free (rgcgnm);
	  free (sodifm);
	  free (srafsm);
	  free (trafsm);
	  free (sraf0m);
	  free (traf0m);
	  free (emtefm);
	  free (trsofm);
	  free (tclfsm);
	  free (sclfsm);
	  free (tclf0m);
	  free (sclf0m);
	  free (tkem);
	  free (tkem1m);
	  free (auxil1m);
	  free (auxil2m);
	  break;
	default:
	  break;
	}
    }

  return (0);
}
