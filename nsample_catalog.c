#include "catalog.h"
#include <getopt.h>		/* for getopt_long */
/* 

   perform a kostrov summation and / or Michael (1984) type stress
   inversion based on AKI or gCMT style earthquake focal
   mechanism/moment tensor catalogs for spatial selection criteria

   there is also bin_catalog which uses simple binning

   nmin, dist_max > 0: use bins with at least nmin entries within dist_max [km]
   nmin < 0, dist_max > 0: use the -nmin nearest entries as long as they are
   within dist_max [km]

   (c) Thorsten Becker, thwbecker@post.harvard.edu, see README

*/

static void usage(char *name, struct kostrov_sum *k, int use_aki,
		  int use_weights, int is_xy, char *out_istring, int fric)
{
  fprintf(stderr,"usage: %s [options] catalog_file\n\n",name);
  fprintf(stderr,"  reads an AKI (default) or CMT focal mechanism catalog, selects\n");
  fprintf(stderr,"  events by nearest neighbor distance criteria, performs a Kostrov\n");
  fprintf(stderr,"  summation, and computes Michael (1984) / Vavrycuk style stress\n");
  fprintf(stderr,"  tensors per sample point. see also bin_catalog for simple binning.\n\n");
  fprintf(stderr,"  the catalog file is the single required argument. for AKI format\n");
  fprintf(stderr,"  the last column is expected to be UNIX time.\n\n");
  fprintf(stderr,"  selection: with --nmin > 0, keep sample points with at least nmin\n");
  fprintf(stderr,"  events within --dist-max km. with --nmin < 0, keep the -nmin\n");
  fprintf(stderr,"  nearest events as long as they are within --dist-max km.\n\n");
  fprintf(stderr,"options (defaults in brackets):\n");
  fprintf(stderr,"  -d, --dx val             sample spacing dx, deg (or km with -x)   [%g]\n",k->dx);
  fprintf(stderr,"  -y, --dy val             sample spacing dy, if different from dx  [dx]\n");
  fprintf(stderr,"  -m, --min-mag val        minimum magnitude                       [%g]\n",k->minmag);
  fprintf(stderr,"  -M, --max-mag val        maximum magnitude                       [%g]\n",k->maxmag);
  fprintf(stderr,"  -l, --min-lon val        minimum longitude                       [%g]\n",k->dlonmin);
  fprintf(stderr,"  -r, --max-lon val        maximum longitude                       [%g]\n",k->dlonmax);
  fprintf(stderr,"  -b, --min-lat val        minimum latitude                        [%g]\n",k->dlatmin);
  fprintf(stderr,"  -t, --max-lat val        maximum latitude                        [%g]\n",k->dlatmax);
  fprintf(stderr,"  -z, --max-depth val      maximum depth                           [%g]\n",k->maxdepth);
  fprintf(stderr,"  -Z, --min-depth val      minimum depth                           [%g]\n",k->mindepth);
  fprintf(stderr,"  -p, --nmin val           min events (>0) or -nmin nearest (<0)   [%i]\n",k->nmin);
  fprintf(stderr,"  -D, --dist-max val       maximum selection distance, km          [%g]\n",k->dist_max);
  fprintf(stderr,"  -w, --weighting val      weighting: 0 none, 1 by distance        [%i]\n",use_weights);
  fprintf(stderr,"  -F, --friction-solve val 1: Vavrycuk stress only, 2/3 optimize friction range,\n");
  fprintf(stderr,"                           4 optimize friction and esti.uncertaint. [%i]\n",fric);
  fprintf(stderr,"  -o, --out-prefix str     output filename prefix                   [%s]\n",out_istring);
  fprintf(stderr,"  -x, --xy                 treat coordinates as Cartesian x y       [%s]\n",is_xy?"on":"off");
  fprintf(stderr,"  -c, --cmt                read CMT format instead of AKI            [%s]\n",use_aki?"off":"on");
  fprintf(stderr,"  -h, --help               print this help and exit\n");
}

int main(int argc, char **argv)
{
  struct cat *catalog;
  struct kostrov_sum *kostrov;
  int itmp,c;
  char out_filename[BC_CHAR_LEN],out_filename2[BC_CHAR_LEN],out_istring[BC_CHAR_LEN];
  char *catalog_file;
  BC_BOOLEAN use_aki = BC_TRUE;
  BC_BOOLEAN remove_trace =  BC_TRUE;	/* remove trace from summations */
  BC_BOOLEAN calc_stress =   BC_TRUE;
  BC_BOOLEAN compute_dtree = BC_TRUE;
  BC_BOOLEAN has_dy = BC_FALSE;		/* was dy set separately from dx? */
  int  use_weights = 0;			/* 0: none 1: by distance from center of bin */
  const int monte_carlo = 0;

  catalog=(struct cat *)calloc(1,sizeof(struct cat));
  if(!catalog)
    BC_MEMERROR(argv[0]);
  /* default parameters */
  kostrov = catalog->sum;
  kostrov_set_defaults(kostrov); /* set binning defaults */

  /* 1: additional stress inversion, 2: optimize friction 0...1, 3: 0.2...0.8, 4: 0...1 and find uncertainties */
  catalog->use_friction_solve = 4;
  snprintf(out_istring,sizeof(out_istring),"nsample");

  /*

     parse the command line flags (short and long forms)

  */
  static struct option long_options[] = {
    {"dx",             required_argument, NULL, 'd'},
    {"dy",             required_argument, NULL, 'y'},
    {"min-mag",        required_argument, NULL, 'm'},
    {"max-mag",        required_argument, NULL, 'M'},
    {"min-lon",        required_argument, NULL, 'l'},
    {"max-lon",        required_argument, NULL, 'r'},
    {"min-lat",        required_argument, NULL, 'b'},
    {"max-lat",        required_argument, NULL, 't'},
    {"max-depth",      required_argument, NULL, 'z'},
    {"min-depth",      required_argument, NULL, 'Z'},
    {"nmin",           required_argument, NULL, 'p'},
    {"dist-max",       required_argument, NULL, 'D'},
    {"weighting",      required_argument, NULL, 'w'},
    {"friction-solve", required_argument, NULL, 'F'},
    {"out-prefix",     required_argument, NULL, 'o'},
    {"xy",             no_argument,       NULL, 'x'},
    {"cmt",            no_argument,       NULL, 'c'},
    {"help",           no_argument,       NULL, 'h'},
    {NULL, 0, NULL, 0}
  };
  while((c = getopt_long(argc,argv,"d:y:m:M:l:r:b:t:z:Z:p:D:w:F:o:xch",
			 long_options,NULL)) != -1){
    switch(c){
    case 'd': sscanf(optarg,BC_PREC_FMT,&kostrov->dx); break;
    case 'y': sscanf(optarg,BC_PREC_FMT,&kostrov->dy); has_dy = BC_TRUE; break;
    case 'm': sscanf(optarg,BC_PREC_FMT,&kostrov->minmag); break;
    case 'M': sscanf(optarg,BC_PREC_FMT,&kostrov->maxmag); break;
    case 'l': sscanf(optarg,BC_PREC_FMT,&kostrov->dlonmin); break;
    case 'r': sscanf(optarg,BC_PREC_FMT,&kostrov->dlonmax); break;
    case 'b': sscanf(optarg,BC_PREC_FMT,&kostrov->dlatmin); break;
    case 't': sscanf(optarg,BC_PREC_FMT,&kostrov->dlatmax); break;
    case 'z': sscanf(optarg,BC_PREC_FMT,&kostrov->maxdepth); break;
    case 'Z': sscanf(optarg,BC_PREC_FMT,&kostrov->mindepth); break;
    case 'p': sscanf(optarg,"%i",&kostrov->nmin); break;
    case 'D': sscanf(optarg,BC_PREC_FMT,&kostrov->dist_max); break;
    case 'w': sscanf(optarg,"%i",&use_weights); break;
    case 'F': sscanf(optarg,"%i",&itmp); catalog->use_friction_solve = itmp; break;
    case 'o': snprintf(out_istring,sizeof(out_istring),"%s",optarg); break;
    case 'x': catalog->is_xy = BC_TRUE; break;
    case 'c': use_aki = BC_FALSE; break;
    case 'h':
    default:
      usage(argv[0],kostrov,(int)use_aki,use_weights,(int)catalog->is_xy,
	    out_istring,catalog->use_friction_solve);
      exit(-1);
    }
  }
  /* the catalog file is the single remaining positional argument */
  if(optind >= argc){
    usage(argv[0],kostrov,(int)use_aki,use_weights,(int)catalog->is_xy,
	  out_istring,catalog->use_friction_solve);
    exit(-1);
  }
  catalog_file = argv[optind];

  if(!has_dy)			/* default dy to dx unless set with -y */
    kostrov->dy = kostrov->dx;

  snprintf(out_filename,sizeof(out_filename),"%s.%g.%g.%i",out_istring,kostrov->dx,kostrov->dy,monte_carlo); /* backward compat */

  fprintf(stderr,"%s: catalog: %s dx: %g dy: %g min_mag: %g max_mag: %g min_lon: %g max_lon: %g min_lat: %g: max_lat: %g\n",
	  argv[0],catalog_file,kostrov->dx,kostrov->dy,kostrov->minmag,kostrov->maxmag,kostrov->dlonmin, kostrov->dlonmax, kostrov->dlatmin, kostrov->dlatmax);
  fprintf(stderr,"%s: min_depth: %g max_depth: %g use_aki: %i use_weights: %i nmin: %i dist_max: %g is_xy: %i out_name: %s\n",
	  argv[0],kostrov->mindepth,kostrov->maxdepth,(int)use_aki,use_weights,
	  kostrov->nmin,kostrov->dist_max,(int)catalog->is_xy,out_filename);

  if(use_aki){			/* aki with last column time */
    fprintf(stderr,"%s: assuming AKI format (last column is UNIX time)\n",argv[0]);
    read_catalog(catalog_file,catalog,BC_AKI,compute_dtree);
  }else{
    fprintf(stderr,"%s: assuming CMT format\n",argv[0]);
    read_catalog(catalog_file,catalog,BC_CMT,compute_dtree);
  }
  if(!catalog->n){
    fprintf(stderr,"%s: error: zero events read\n",argv[0]);
    exit(-1);
  }

  /*

     setup bins

  */
  setup_kostrov(catalog,use_weights);
  /*
     assemble based on nmin and dist_max criteria
  */
  assemble_bins_based_on_distance(catalog,remove_trace,monte_carlo,BC_TRUE);
  /*
     print non-zero
  */
  print_kostrov_bins(catalog,out_filename,monte_carlo);

  if(calc_stress){
    /* compute Andy Michael style stress tensors */
    snprintf(out_filename2,sizeof(out_filename2),"%s.%g.%g",out_istring,kostrov->dx,kostrov->dy);
    calc_stress_tensor_for_kbins(catalog);
    print_stress_tensors(catalog,out_filename2);
  }
  geo_search_destroy(catalog->tree);
  return 0;
}
