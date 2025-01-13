/* tn.c */
extern int i4_uniform_ab(int a, int b, int *seed);
extern double normal_01_cdf(double x);
extern double normal_01_cdf_inv(double p);
extern void normal_01_cdf_values(int *n_data, double *x, double *fx);
extern double normal_01_mean(void);
extern double normal_01_moment(int order);
extern double normal_01_pdf(double x);
extern double normal_01_sample(int *seed);
extern double normal_01_variance(void);
extern double normal_ms_cdf(double x, double mu, double sigma);
extern double normal_ms_cdf_inv(double cdf, double mu, double sigma);
extern double normal_ms_mean(double mu, double sigma);
extern double normal_ms_moment(int order, double mu, double sigma);
extern double normal_ms_moment_central(int order, double mu, double sigma);
extern double normal_ms_moment_central_values(int order, double mu, double sigma);
extern double normal_ms_moment_values(int order, double mu, double sigma);
extern double normal_ms_pdf(double x, double mu, double sigma);
extern double normal_ms_sample(double mu, double sigma, int *seed);
extern double normal_ms_variance(double mu, double sigma);
extern double r8_abs(double x);
extern double r8_choose(int n, int k);
extern double r8_factorial2(int n);
extern void r8_factorial2_values(int *n_data, int *n, double *f);
extern double r8_huge(void);
extern double r8_log_2(double x);
extern double r8_mop(int i);
extern double r8_uniform_01(int *seed);
extern void r8poly_print(int n, double a[], char *title);
extern double r8poly_value_horner(int m, double c[], double x);
extern double *r8vec_linspace_new(int n, double a, double b);
extern double r8vec_max(int n, double *dvec);
extern double r8vec_mean(int n, double x[]);
extern double r8vec_min(int n, double *dvec);
extern void r8vec_print(int n, double a[], char *title);
extern double r8vec_variance(int n, double x[]);
extern void timestamp(void);
extern double truncated_normal_ab_cdf(double x, double mu, double sigma, double a, double b);
extern void truncated_normal_ab_cdf_values(int *n_data, double *mu, double *sigma, double *a, double *b, double *x, double *fx);
extern double truncated_normal_ab_cdf_inv(double cdf, double mu, double sigma, double a, double b);
extern double truncated_normal_ab_mean(double mu, double sigma, double a, double b);
extern double truncated_normal_ab_moment(int order, double mu, double sigma, double a, double b);
extern double truncated_normal_ab_pdf(double x, double mu, double sigma, double a, double b);
extern void truncated_normal_ab_pdf_values(int *n_data, double *mu, double *sigma, double *a, double *b, double *x, double *fx);
extern double truncated_normal_ab_sample(double mu, double sigma, double a, double b, int *seed);
extern double truncated_normal_ab_variance(double mu, double sigma, double a, double b);
extern double truncated_normal_a_cdf(double x, double mu, double sigma, double a);
extern void truncated_normal_a_cdf_values(int *n_data, double *mu, double *sigma, double *a, double *x, double *fx);
extern double truncated_normal_a_cdf_inv(double cdf, double mu, double sigma, double a);
extern double truncated_normal_a_mean(double mu, double sigma, double a);
extern double truncated_normal_a_moment(int order, double mu, double sigma, double a);
extern double truncated_normal_a_pdf(double x, double mu, double sigma, double a);
extern void truncated_normal_a_pdf_values(int *n_data, double *mu, double *sigma, double *a, double *x, double *fx);
extern double truncated_normal_a_sample(double mu, double sigma, double a, int *seed);
extern double truncated_normal_a_variance(double mu, double sigma, double a);
extern double truncated_normal_b_cdf(double x, double mu, double sigma, double b);
extern void truncated_normal_b_cdf_values(int *n_data, double *mu, double *sigma, double *b, double *x, double *fx);
extern double truncated_normal_b_cdf_inv(double cdf, double mu, double sigma, double b);
extern double truncated_normal_b_mean(double mu, double sigma, double b);
extern double truncated_normal_b_moment(int order, double mu, double sigma, double b);
extern double truncated_normal_b_pdf(double x, double mu, double sigma, double b);
extern void truncated_normal_b_pdf_values(int *n_data, double *mu, double *sigma, double *b, double *x, double *fx);
extern double truncated_normal_b_sample(double mu, double sigma, double b, int *seed);
extern double truncated_normal_b_variance(double mu, double sigma, double b);