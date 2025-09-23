/* bin_catalog.c */
/* bvalue.c */
void calc_b_value_bpos_mc(double *, long int, double, double, double *);
void calc_b_value_bpos(double *, long int, double, double *);
void calc_b_value_marzocci(double *, long int, double, double, double *, double *);
void calc_b_value_thomas(double *, long int, double, double, double *, double *);
void calc_b_value_ml(double *, long int, double, double *, double *);
/* calc_gr.c */
/* calc_gr_time.c */
void calc_gr_switch(double *, long, double, double, double *, double *, int);
/* eigen.c */
void calc_eigensystem_vec6(double *, double *, double *, unsigned short, unsigned short);
void calc_eigensystem_sym_3x3(double [3][3], double *, double *, unsigned short, unsigned short);
void calc_eigensystem_sym_9(double *, double *, double *, unsigned short, unsigned short);
void indexx(int, double *, int *);
/* eigen_3ds_simple.c */
/* eigen_driver.c */
/* fault_eq.c */
void stridip(double, double, double, double *, double *);
void find_alt_plane(double, double, double, double *, double *, double *);
void aki2mom(double, double, double, double, double *, double *);
double mean_hor_strain(double *);
double mag2mom(double);
double mom2mag(double);
double mag2pot(double);
double scalar_mom(double *);
/* geo_kdtree.c */
geo_tree_t *geo_tree_create(int);
void geo_tree_destroy(geo_tree_t *);
int allocate_node(geo_tree_t *);
int geo_tree_add_point(geo_tree_t *, double, double, int);
void geo_tree_build(geo_tree_t *);
int build_tree_recursive(geo_tree_t *, int *, int, int);
result_array_t *geo_tree_query_radius(geo_tree_t *, double, double, double, double, int);
void search_radius_recursive(geo_tree_t *, int, double, double, double, double, result_array_t *, int);
result_array_t *geo_tree_query_k_nearest(geo_tree_t *, double, double, int);
void search_knn_recursive(geo_tree_t *, int, double, double, query_result_t *, int *, int);
result_array_t *result_array_create(int);
void result_array_add(result_array_t *, kd_tree_point, double);
void result_array_destroy(result_array_t *);
int compare_by_val(const void *, const void *);
int compare_by_distance(const void *, const void *);
void print_results(result_array_t *, const char *);
/* handle_catalog.c */
void relocate_catalog(struct cat *, struct cat *, struct cat *);
void kostrov_set_defaults(struct kostrov_sum *);
void sum_kostrov_bins(struct cat *, unsigned short, unsigned short, unsigned short);
void assemble_bins_based_on_distance(struct cat *, unsigned short, unsigned short);
void sum_smoothed_seismicity(struct cat *, int);
void make_histogram(double *, int, double, double *, double *, int **, double **, int *);
void print_histogram(int *, double *, int, FILE *);
void setup_kostrov(struct cat *, int);
void clear_bins(struct cat *);
void print_kostrov_bins(struct cat *, char *, unsigned short);
void print_stress_tensors(struct cat *, char *);
double kostrov_bdlon(int, struct kostrov_sum *);
double kostrov_bdlat(int, struct kostrov_sum *);
void print_summed_moment(struct cat *, char *);
void merge_catalog(struct cat *, struct cat *, struct cat *, double, double);
int print_catalog(char *, struct cat *, int);
int read_catalog(char *, struct cat *, int, unsigned short);
void copy_quake(struct qke *, struct qke *);
char *mode_name(int);
int read_quake(FILE *, struct qke *, int);
void print_quake(FILE *, struct qke, int);
int read_quake_aki(FILE *, struct qke *);
int read_quake_cmt(FILE *, struct qke *);
int read_quake_eng(FILE *, struct qke *);
void print_quake_aki(FILE *, struct qke);
void print_quake_cmt(FILE *, struct qke);
void print_quake_cmt_fp(FILE *, struct qke);
void create_catalog(struct cat *, long int);
void make_room_for_quake(struct cat *);
double quake_weight(double, double, double, int);
double quake_scale_lkm(double);
void add_quake_to_bin_list(unsigned int, struct bn *, double);
void assign_quake_angles(struct qke *, double *);
void swap_angles(double *);
void swap(double *, double *);
unsigned short quake_qualified(double, double, double, double, double, double);
/* handle_catalog_gmt.c */
void tensor2fpangle(double *, double *, double *, double *, double *, double *, double *);
void GMT_momten2axe(struct M_TENSOR, struct AXIS *, struct AXIS *, struct AXIS *);
void axe2dc(struct AXIS, struct AXIS, struct nodal_plane *, struct nodal_plane *);
double computed_rake2(double, double, double, double, double);
int GMT_jacobi(double *, int *, int *, double *, double *, double *, double *, int *);
/* linalg_misc_geo.c */
void ranger(double *);
double max_x_from_int_vector(double *, int *, int);
void rotate_vec6(double *, double *, double, double, double);
void sixsymtomat(double *, double [3][3]);
void mattosixsym(double [3][3], double *);
void get_gen_rot(double [3][3], double, double, double);
void rotate_mat(double [3][3], double [3][3], double [3][3]);
void remove_trace(double *);
double distance_cart(double, double, double, double);
double distance_geo(double, double, double, double, double, double);
double distance(struct cat *, struct cat *, int, int);
void get_index_vector(int **, int, int, long *);
void normalize_tens6(double *);
double tensor6_norm(double *);
void tens6to3by3(double *, double [3][3]);
double std_quick(int, double, double);
double ran2(long int *);
double gauss_ran(long int *, double);
double gasdev(long int *);
FILE *myopen(char *, char *, char *);
void sincos(double, double *, double *);
/* m02dcfp.c */
/* m02mag.c */
/* merge_catalog.c */
/* michael_leasq.c */
void michael_leasq(double *, int, int, double *, double *, double *, double *, double *);
void michael_gaus(double *, int, double *, double *);
double dabs(double);
void michael_atransa(double *, int, int, double *);
void michael_atransb(double *, int, int, double *, double *);
void michael_sigsq(double *, int, int, double *, double *, double *);
/* nsample_catalog.c */
/* slip_deviation.c */
void calc_misfits_from_single_angle_set(double *, double *, int, double *);
void slip_deviation_svec_single(double *, double *, double *, double *);
void slip_deviation_mmat_single(double [3][3], double *, double *, double *);
double slip_deviation_dotp(double *, double *, double [3][3]);
/* solve_stress_one_bin.c */
/* stability_criterion.c */
void optimize_angles_via_instability(int, double *, double *, double, double *, double *, int *);
void calc_average_instability(int, double *, double *, double, double *, double *);
void stability_criterion_eig(double *, double *, double, double *, unsigned short, double *);
/* stress_inversion.c */
void calc_stress_tensor_for_kbins(struct cat *);
void solve_stress_michael_random_sweep(int, double *, double *, double *, double *, long int *);
void adjust_stress_for_friction(int, double *, double *, double *, double *, double *, double *, double *, double *, unsigned short, unsigned short, int *, int *);
void solve_stress_michael_specified_plane(int, double *, double *, double *);
void michael_solve_lsq(int, int, int, double *, double *, double *, double *);
void my6stress2m3x3(double *, double [3][3]);
void michael_assign_to_matrix(double *, int *, double **, double **);
/* test_eigen.c */
/* test_kdtree.c */
void demo(void);
