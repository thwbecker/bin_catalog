#include "catalog.h"
#include <getopt.h>		/* for getopt_long */
/* 

   perform a kostrov summation and / or Michael (1984) type stress
   inversion based on AKI or gCMT style earthquake focal
   mechanism/moment tensor catalogs for simple binning

   there is also nsample_catalog which tries to find the closest
   events and has more flexibility in terms of sampling
   
   since June 2026, some Claude code additions

   (c) 2020 - 2026, Thorsten Becker, thwbecker@post.harvard.edu, see
   README
*/

static void usage(char *name, struct kostrov_sum *k, int monte_carlo,
		  int use_aki, int weighting_method, int is_xy,
		  char *out_istring, int nmin, int fric)
{
  fprintf(stderr,"usage: %s [options] catalog_file\n\n",name);
  fprintf(stderr,"  reads an AKI (default) or CMT focal mechanism catalog, bins it,\n");
  fprintf(stderr,"  performs a Kostrov summation, and computes Michael (1984) /\n");
  fprintf(stderr,"  Vavrycuk style stress tensors per bin.\n\n");
  fprintf(stderr,"  the catalog file is the single required argument. For AKI format\n");
  fprintf(stderr,"  the last column is expected to be UNIX time (non standard).\n\n");
  fprintf(stderr,"options (defaults in brackets):\n");
  fprintf(stderr,"  -d, --dx val             bin spacing dx, deg (or km with -x)      [%g]\n",k->dx);
  fprintf(stderr,"  -y, --dy val             bin spacing dy, if different from dx     [dx]\n");
  fprintf(stderr,"  -m, --min-mag val        minimum magnitude                        [%g]\n",k->minmag);
  fprintf(stderr,"  -M, --max-mag val        maximum magnitude                        [%g]\n",k->maxmag);
  fprintf(stderr,"  -l, --min-lon val        minimum longitude                        [%g]\n",k->dlonmin);
  fprintf(stderr,"  -r, --max-lon val        maximum longitude                        [%g]\n",k->dlonmax);
  fprintf(stderr,"  -b, --min-lat val        minimum latitude                         [%g]\n",k->dlatmin);
  fprintf(stderr,"  -t, --max-lat val        maximum latitude                         [%g]\n",k->dlatmax);
  fprintf(stderr,"  -z, --max-depth val      maximum depth                            [%g]\n",k->maxdepth);
  fprintf(stderr,"  -Z, --min-depth val      minimum depth                            [%g]\n",k->mindepth);
  fprintf(stderr,"  -n, --monte-carlo val    number of Monte Carlo realizations       [%i]\n",monte_carlo);
  fprintf(stderr,"  -w, --weighting val      weighting: 0 none, 1 time, 2 frequency   [%i]\n",weighting_method);
  fprintf(stderr,"  -p, --min-events val     minimum events per bin for stress        [%i]\n",nmin);
  fprintf(stderr,"  -F, --friction-solve val 1: Vavrycuk stress only, 2/3 optimize friction range,\n");
  fprintf(stderr,"                           4 optimize friction and esti.uncertaint. [%i]\n",fric);
  fprintf(stderr,"  -o, --out-prefix str     output filename prefix                   [%s]\n",out_istring);
  fprintf(stderr,"  -x, --xy                 treat coordinates as Cartesian x y       [%s]\n",is_xy?"on":"off");
  fprintf(stderr,"  -c, --cmt                read CMT format instead of AKI           [%s]\n",use_aki?"off":"on");
  fprintf(stderr,"  -h, --help               print this help and exit\n");
}

int main(int argc, char **argv)
{
  struct cat *catalog;
  struct kostrov_sum *kostrov;
  int itmp,c;
  char out_filename[BC_CHAR_LEN],out_filename2[BC_CHAR_LEN],out_istring[BC_CHAR_LEN];
  char *catalog_file;
  int monte_carlo = 0;
  BC_BOOLEAN use_aki     =   BC_TRUE;
  BC_BOOLEAN remove_trace =  BC_TRUE;	/* remove trace from summations */
  BC_BOOLEAN calc_stress =   BC_TRUE;
  BC_BOOLEAN compute_dtree = BC_FALSE;	/* not needed */
  BC_BOOLEAN has_dy = BC_FALSE;		/* was dy set separately from dx? */
  /* default */
  int min_events_for_stress = 5; /* at least so many events to attempt a stress inversion */
  int  weighting_method = 0; /* 0: none
				1: distance from center of time
				2: normalized by numbers over time
			     */
  catalog=(struct cat *)calloc(1,sizeof(struct cat));
  if(!catalog)
    BC_MEMERROR(argv[0]);

  /* 1: additional stress inversion, 2: optimize friction 0...1, 
     3: 0.2...0.8,                   4: 0...1 and find uncertainties */
  catalog->use_friction_solve = 2;
  sprintf(out_istring,"kostrov");

  kostrov = catalog->sum;
  kostrov_set_defaults(kostrov); /* set binning defaults */

  /*

     parse the command line flags (short and long forms)

  */
  static struct option long_options[] = {
    {"dx",             required_argument, NULL, 'd'},
    {"dy",             required_argument, NULL, 'y'},
    {"min-mag",        required_argument, NULL, 'm'},
    {"max-mag",        required_argument, NULL, 'M'},
    {"monte-carlo",    required_argument, NULL, 'n'},
    {"min-lon",        required_argument, NULL, 'l'},
    {"max-lon",        required_argument, NULL, 'r'},
    {"min-lat",        required_argument, NULL, 'b'},
    {"max-lat",        required_argument, NULL, 't'},
    {"max-depth",      required_argument, NULL, 'z'},
    {"min-depth",      required_argument, NULL, 'Z'},
    {"weighting",      required_argument, NULL, 'w'},
    {"min-events",     required_argument, NULL, 'p'},
    {"friction-solve", required_argument, NULL, 'F'},
    {"out-prefix",     required_argument, NULL, 'o'},
    {"xy",             no_argument,       NULL, 'x'},
    {"cmt",            no_argument,       NULL, 'c'},
    {"help",           no_argument,       NULL, 'h'},
    {NULL, 0, NULL, 0}
  };
  while((c = getopt_long(argc,argv,"d:y:m:M:n:l:r:b:t:z:Z:w:xco:p:F:h",
			 long_options,NULL)) != -1){
    switch(c){
    case 'd': sscanf(optarg,BC_PREC_FMT,&kostrov->dx); break;
    case 'y': sscanf(optarg,BC_PREC_FMT,&kostrov->dy); has_dy = BC_TRUE; break;
    case 'm': sscanf(optarg,BC_PREC_FMT,&kostrov->minmag); break;
    case 'M': sscanf(optarg,BC_PREC_FMT,&kostrov->maxmag); break;
    case 'n': sscanf(optarg,"%i",&monte_carlo); break;
    case 'l': sscanf(optarg,BC_PREC_FMT,&kostrov->dlonmin); break;
    case 'r': sscanf(optarg,BC_PREC_FMT,&kostrov->dlonmax); break;
    case 'b': sscanf(optarg,BC_PREC_FMT,&kostrov->dlatmin); break;
    case 't': sscanf(optarg,BC_PREC_FMT,&kostrov->dlatmax); break;
    case 'z': sscanf(optarg,BC_PREC_FMT,&kostrov->maxdepth); break;
    case 'Z': sscanf(optarg,BC_PREC_FMT,&kostrov->mindepth); break;
    case 'w': sscanf(optarg,"%i",&weighting_method); break;
    case 'p': sscanf(optarg,"%i",&min_events_for_stress); break;
    case 'F': sscanf(optarg,"%i",&itmp); catalog->use_friction_solve = itmp; break;
    case 'o': sprintf(out_istring,"%s",optarg); break;
    case 'x': catalog->is_xy = BC_TRUE; break;
    case 'c': use_aki = BC_FALSE; break;
    case 'h':
    default:
      usage(argv[0],kostrov,monte_carlo,(int)use_aki,weighting_method,
	    (int)catalog->is_xy,out_istring,min_events_for_stress,
	    catalog->use_friction_solve);
      exit(-1);
    }
  }
  /* the catalog file is the single remaining positional argument */
  if(optind >= argc){
    usage(argv[0],kostrov,monte_carlo,(int)use_aki,weighting_method,
	  (int)catalog->is_xy,out_istring,min_events_for_stress,
	  catalog->use_friction_solve);
    exit(-1);
  }
  catalog_file = argv[optind];

  if(!has_dy)			/* default dy to dx unless set with -y */
    kostrov->dy = kostrov->dx;
  kostrov->nmin = min_events_for_stress;

  /* output filename for Kostrov summations */
  snprintf(out_filename,sizeof(out_filename),"%s.%g.%g.%i",out_istring,kostrov->dx,kostrov->dy,monte_carlo);

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
  switch(weighting_method){
  case 0:
    break;
  case 1:
    fprintf(stderr,"%s: weighing sum by temporal distance from tsec %g\n",
	    argv[0],catalog->tcenter);
    break;
  case 2:
    fprintf(stderr,"%s: frequency normalized weight\n",argv[0]);
    break;
  default:
    fprintf(stderr,"%s: normalizing method %i undefined\n",argv[0],weighting_method);
    exit(-1);
    break;
  }

  /*

     setup bins

  */
  setup_kostrov(catalog,weighting_method);

  /*
     sum
  */
  sum_kostrov_bins(catalog,remove_trace,monte_carlo,BC_TRUE);
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
  return 0;
}
