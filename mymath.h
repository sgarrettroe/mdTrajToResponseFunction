
void nrerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
float *****f5tensor(long n1l,long n1h,long n2l,long n2h,long n3l,long n3h,long n4l,long n4h,long n5l,long n5h);
float **submatrix(float **a,long oldrl,long oldrh,long oldcl, long oldch,long newrl,long newcl);
void free_submatrix(float **b,long nrl,long nrh,long ncl,long nch);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m,long nrl,long nrh,long ncl,long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh);//
void free_f5tensor(float *****t,long n1l,long n1h,long n2l,long n2h,long n3l,long n3h,long n4l,long n4h,long n5l,long n5h);

//cut the rest of these?
float ran2(long *idum);
void four1(float data[], unsigned long nn, int isign);
void realft(float data[], unsigned long n, int isign);
void twofft(float data1[], float data2[], float fft1[], float fft2[], unsigned long n);
void correl(float data1[],float data2[], unsigned long n, float ans[]);
void autocorrel(float data[], unsigned long n, float ans[]);
void fft(float *P_re, float *P_im, int n);
void fft3D(float ***P_re, float ***P_im, int nt);
void fft2D(float **P_re, float **P_im, int nt);
int sign(float x);

float ** hist(float seq[],unsigned long n_steps,unsigned long n_edges);
void histEdges(float seq[],unsigned long n_steps,unsigned long n_edges,float x[]);
void histc(float seq[],unsigned long n,float x_edges[],float h[],unsigned long n_edges,unsigned long bin[]);
void jointProbabilityDensity3D(float seq[],unsigned long n_steps,float edges[],unsigned long n_edges,unsigned long lag1,unsigned long lag2,float ***rho_3);
void threePointCorrelation(float seq[],unsigned long n_steps,unsigned long lags[],unsigned long n_lags,float **c3);
void histTri(float x[],float y[],float z[],unsigned long n_steps,float edges[],unsigned long n_edges,float ***hist_tensor);

