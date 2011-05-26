#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h> //added for long options

// 14.9.05: Filled Fouriertransform with linear interpolation, reduces dramatically the offset problem
/*  
*  Original program fifthorder.c by P. Hamm calculated response functions from 
*  Langevin dynamics trajectories. 
*  8.1.2007 Modified by S. Garrett-Roe to read an MD trajectory (one file of coordinates and 
                                                                 *  one file of forces) and calculate response functions from that.
*
*/

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
void normalizeResults(int nt,unsigned long isample,
                      float S_re[],float S_im[],
                      float **P1_re,float **P1_im,float **P2_re,float **P2_im,
                      float ***R1_re,float ***R1_im,float ***R2_re,float ***R2_im,
                      float ***R3_re,float ***R3_im,float ***R4_re,float ***R4_im,
                      float Sf_re[],float Sf_im[],
                      float **P1f_re,float **P1f_im,float **P2f_re,float **P2f_im,
                      float ***R1f_re,float ***R1f_im,float ***R2f_re,float ***R2f_im,
                      float ***R3f_re,float ***R3f_im,float ***R4f_re,float ***R4f_im);
void fourierTransformResults(int nt,
                             float Sf_re[],float Sf_im[],
                             float **P1f_re,float **P1f_im,float **P2f_re,float **P2f_im,
                             float ***R1f_re,float ***R1f_im,float ***R2f_re,float ***R2f_im,
                             float ***R3f_re,float ***R3f_im,float ***R4f_re,float ***R4f_im);
void writeResults(char base_name[],int nt,char extension[],
                  float Sf_re[],float Sf_im[],
                  float **P1f_re,float **P1f_im,float **P2f_re,float **P2f_im,
                  float ***R1f_re,float ***R1f_im,float ***R2f_re,float ***R2f_im,
                  float ***R3f_re,float ***R3f_im,float ***R4f_re,float ***R4f_im);



int DEBUG_LEVEL=0;
int LOG_LEVEL=1;

//#define PI 3.1415927 //use M_PI instead!
#define order 5
#define wavenumbersToInvPs 2.99792458e-2
//#define DT 0.01 //ps, time-steps output by GROMACS //read mdp file

int 
main(int iarg, char **input) 
{
  time_t my_time=time(0); // time process started
  
  FILE *param_fid,*fid;
  char name[100],freq_file_name[100],time_file_name[100],mdp_file_name[100];
  char base_name[100],string[100];
  unsigned long i,j,l,k,m;
  int have_mdp_file=0;

  //vars for force trajectory
  unsigned long nsteps=0, natoms, nmols, nprotons=0; //these set how the
						   //file is read
  unsigned long nstepsmax=0,nprotonsmax=0; //these set how much of the
					   //trajectory is averaged
  unsigned long nstxout=0;
  float dt=0;//stepsize read from mdp or commandline
  
  //stuff for the calculating
  double sum;
  float **dw_matrix;
  float *dw; //a pointer into dw_matrix representing the frequency
             //trajectory of one proton
  float *dw_initial; //the initial frequency of each proton
  float mean,sigma,skewness; //characterize the frequency traj
  
  //vars for histograms
  unsigned long n_bins = 64; //bins in the histogram and joint prob dens
  unsigned long *bin_traj; //which histogram bin does each trajectory step fall into
  float *P_w; //the probability distribution
  float *x_edges; //defines the bins into which frequencies are grouped
  
  //vars for c2 and c3
  unsigned long n_c2;
  unsigned long n_zeropad; //length of the zeropad
  unsigned long n_c2_output; //number of points to actually output 
  unsigned long n_c3 = 64; //number of points per axis of c3
  unsigned long *lags; //vector of delays used in calculating c3
  float *c2;//2-point correlation function of one proton
  float *c2_accum; //averaged over protons
  float **c3; //3-point c.f. of one proton (not used)
  float **c3_accum; // 3-point c.f. averaged over protons
  float c2_max = 1; // limit of c2 in ps!
  float c3_max = 0.640; // limit of c3 in ps!
  int c3_lag_stepsize; // in steps


  //vars for the time axes
  //const unsigned long n_t2_array = 3; //number of time steps
  //unsigned long *t2_array,*tdiag_array; //arrays of time steps used
  unsigned long n_t2_t4_pairs; //number of times we will calculate joint prob dens
  float **t2_t4_pairs;//matrix of times (t2,t4)
  float *****rho_5tensor; //joint prob densities

  /*start processing command line*/
  opterr = 0;
  int c; /* for getopt */
  int option_index = 0;
  static struct option long_options[] =
    {
      /* any options that set a flag go here */
      //{"help",no_argument,&help_flag,1},
      //{"brief",no_argument,&flag,0},
      /* these options don't set a flag, we distinguish them by index*/
      {"n_c3",required_argument,0,'a'},
      {"n_c2_output",required_argument,0,'b'},
      {"nbins",required_argument,0,'c'},
      {"debug_level",required_argument,0,'d'},
      {"c3_max",required_argument,0,'e'},
      {"help",no_argument,0,'h'},
      {"log_level",required_argument,0,'l'},
      {"mdp",required_argument,0,'m'},
      {"nmols",required_argument,0,'n'},
      {"nprotons",required_argument,0,'p'},
      {"max_protons_to_use",required_argument,0,'q'},
      {"nsteps",required_argument,0,'s'},
      {"max_steps_to_use",required_argument,0,'t'},
      {0,0,0,0}
    };

  while ((c = getopt_long(iarg, input, "a:b:c:d:h:m:n:p:s:t:q:l",
			  long_options, &option_index))!=-1)
    switch(c){
    case 'a': /*n_c3*/
      if (sscanf(optarg,"%lu",&n_c3)!=1) nrerror("failed reading commandline option -a or --n_c3");
      break;
    case 'b': /*n_c2_output*/
      if (sscanf(optarg,"%lu",&n_c2_output)!=1) nrerror("failed reading commandline option -b or --n_c2_output");
      break;
    case 'c': /*bins*/
      if (sscanf(optarg,"%lu",&n_bins)!=1) nrerror("failed reading commandline option -c or --nbins");
      break;      
    case 'd': /* set debug level */
      if(sscanf(optarg,"%i",&DEBUG_LEVEL)!=1) nrerror("failed reading commandline option");
      break;
    case 'h': /* help and exit */
      printf("USAGE:\n   %s freq_file time_file output_file_prefix\n\n",input[0]);
      printf("optional arguments:\n");
      printf("\t-a NUM, --n_c3 NUM\n");
      printf("\t-b NUM, --n_c2_output NUM\n");
      printf("\t-c NUM, --nbins NUM\n");
      printf("\t-d NUM, --debug_level NUM\n");
      printf("\t-e NUM, --c3_max NUM\n");
      printf("\t-h NUM, --help\n");
      printf("\t-l NUM, --log_level NUM\n");
      printf("\t-m FILE, --mdp FILE\n");
      printf("\t-n NUM, --nmols NUM\n");
      printf("\t-p NUM, --nprotons NUM\n");
      printf("\t-q NUM, --max_protons_to_use NUM\n");
      printf("\t-s NUM, --nsteps NUM\n");
      printf("\t-t NUM, --max_steps_to_use NUM\n");
      exit(EXIT_SUCCESS);
    case 'l': /* set log level */
      if(sscanf(optarg,"%i",&LOG_LEVEL)!=1) nrerror("failed reading commandline option");
      break;
    case 'm': /* have mdp file */
      have_mdp_file = 1; //true
      strcpy(mdp_file_name,optarg);
      break;
    case 'n': /* set nproton number of protons in trajectory */
      if(sscanf(optarg,"%lu",&nmols)!=1) nrerror("failed reading commandline option");
      nprotons = 2*nmols;
      natoms = nmols*3;
      break;
    case 'p': /* set nproton number of protons in trajectory */
      if(sscanf(optarg,"%lu",&nprotons)!=1) nrerror("failed reading commandline option");
      nmols = nprotons/2;
      natoms = nmols*3;
      break;
    case 'q': /* max protons to use */
      if(sscanf(optarg,"%lu",&nprotonsmax)!=1) nrerror("failed reading commandline option");
      break;
    case 's': /* steps to read */
      if(sscanf(optarg,"%lu",&nsteps)!=1) nrerror("failed reading commandline option");
      break;
    case 't': /* max steps to use */
      if(sscanf(optarg,"%lu",&nstepsmax)!=1) nrerror("failed reading commandline option");
      break;
    case '?':
      nrerror("unknown command line option.");
    default:
      abort();
    }
  
  if(iarg-optind<3)
    {
      printf("Too few input arguments. Please give 3 file names.\n");
      printf("USAGE:\n   %s freq_file time_file output_file_prefix\n\n",input[0]);
      exit(EXIT_FAILURE);
    }

  //get the file names
  strcpy(freq_file_name,input[0+optind]);
  strcpy(time_file_name,input[1+optind]);
  strcpy(base_name,input[2+optind]);
  
  /* try to read the mdp file if the option was passed */
  if (have_mdp_file){
    if (LOG_LEVEL>=1)
      printf("Reading mdp file ...\n");

    //test for file existence
    fid = fopen(mdp_file_name,"rt");
    if (fid==NULL)
      nrerror("Can't open mdp file");
    fclose(fid);
    if (LOG_LEVEL>=1)
      printf("\t ... mdp file exists\n");

    //now start to look for variables
    if (LOG_LEVEL>=1)
      printf("\t ... nsteps ");
    strcpy(string,"grep nsteps ");
    strcat(string,mdp_file_name);
    fid = popen(string,"r");
    if(fscanf(fid,"nsteps =  %lu",&nsteps)!=1)
      nrerror("Can't read total nsteps from mdp file!");
    pclose(fid);
    if (LOG_LEVEL>=1)
      printf("%lu\n",nsteps);

    if (LOG_LEVEL>=1)
      printf("\t ... nstxout ");
    strcpy(string,"grep nstxout ");
    strcat(string,mdp_file_name);
    fid = popen(string,"r");
    if(fscanf(fid,"nstxout =  %lu",&nstxout)!=1)
      nrerror("Can't read nstxout from mdp file!");
    pclose(fid);
    if (LOG_LEVEL>=1)
      printf("%lu\n",nstxout);

    //calculate the number of steps that should be in the dw file
    nsteps = nsteps/nstxout+1; //don't forget the 0
    if (LOG_LEVEL>=1)
	printf("\t ... using %lu\n",nsteps);
    
    if (LOG_LEVEL>=1)
      printf("\t ... dt ");
    strcpy(string,"grep dt ");
    strcat(string,mdp_file_name);
    fid = popen(string,"r");
    if(fscanf(fid,"dt =  %f",&dt)!=1)
      nrerror("Can't read nsteps from mdp file!");
    pclose(fid);
    if (LOG_LEVEL>=1)
      printf("%f\n",dt);
    
  }//end read mdp file

  /* if not set at the command line, or the mdp file use the shell command "wc -l"
   * to scan the length of the coord file to get the number of steps
   * and read one line with "head" and then "wc -w" to get the number of molecules
   */
  if (nsteps == 0){
    if (LOG_LEVEL>0)
      printf("Determine number of steps from frequency file.\n");
    strcpy(string,"wc -l ");
    strcat(string,freq_file_name);
    fid = popen(string,"r");
    if(fid==NULL) nrerror("Opening 'wc -l freq_file_name' pipe failed.");
    fscanf(fid,"%ld",&nsteps);
    pclose(fid);
    if (LOG_LEVEL>0)
      printf("Found %ld steps in file %s\n", nsteps, freq_file_name);
  }
    
  if (nprotons==0){
    if (LOG_LEVEL>0)
      printf("Determine number of molecules from frequency file.\n");
    strcpy(string,"head -n 1 ");
    strcat(string,freq_file_name);
    strcat(string," | wc -w");
    fid = popen(string,"r");
    if(fid==NULL) nrerror("Opening 'head -n 1 freq_file_name | wc -w' pipe failed.");
    fscanf(fid,"%ld",&nprotons);
    pclose(fid);
    nmols = nprotons/2;
    natoms = nmols*3;
    if (LOG_LEVEL>0)
      printf("Found %ld protons in file %s\n", nprotons, freq_file_name);
  }

  /*
   * if not otherwise specified at the command line, set nstepsmax
   * and nprotonsmax to nsteps and nprotons 
   */
  if (nstepsmax==0) nstepsmax = nsteps;
  if (nprotonsmax==0) nprotonsmax = nprotons;
    
  //calculate lengths of the two-point corr fxn which must be a
  //power of two because of the FFT
  n_c2 = pow(2,ceil(log(nstepsmax)/log(2))); //the next largest power of 2
  n_zeropad = n_c2 - nstepsmax;
  n_c2_output = round(c2_max/dt);

  /* setup the c3 axis */
  c3_lag_stepsize = round(c3_max / dt / n_c3 );

  /* try to figure out how many time points for joint probability
     density */
  strcpy(string,"wc -l ");
  strcat(string,time_file_name);
  fid=popen(string,"r");
  if(fid==NULL) nrerror("Failed opening pipe wc -l time_file_name");
  fscanf(fid,"%ld",&n_t2_t4_pairs);
  pclose(fid);
  
  // write parameter file for what we are going to calculate
  if (LOG_LEVEL>0)
    printf("Writing parameter file\n");
  strcpy(name,base_name);
  strcat(name,"_parameters_stats.txt");
  param_fid = fopen(name,"wt");
  fprintf(param_fid,"# This is %s called as\n",input[0]);
  fprintf(param_fid,"# ");
  for (i=0;i<iarg;i++)
    fprintf(param_fid,"%s ",input[i]);
  fprintf(param_fid,"\n");
  fprintf(param_fid,"# On: %s\n",ctime(&my_time));
  fprintf(param_fid,"\n"); 
  //write more parameters
  fprintf(param_fid,"# Parameters:\n");
  fprintf(param_fid,"dt = %f\n",dt);
  //fprintf(param_fid,"NSKIP = %i\n", 0);
  fprintf(param_fid,"n_bins = %li\n", n_bins);
  fprintf(param_fid,"n_c2_output = %li\n", n_c2_output);
  fprintf(param_fid,"nsteps = %i\n", (int) nsteps);
  fprintf(param_fid,"nmols = %i\n", (int) nmols);
  fprintf(param_fid,"nstepsmax = %i\n", (int) nstepsmax);
  fprintf(param_fid,"nprotonsmax = %i\n", (int) nprotonsmax);
  fprintf(param_fid,"n_c2 = %i\n",(int) n_c2);
  fprintf(param_fid,"n_c3 = %i\n",(int) n_c3);
  fprintf(param_fid,"c3_stepsize = %f\n",(float)c3_lag_stepsize*dt);
  fprintf(param_fid,"c3_max = %f\n",c3_max);
  fprintf(param_fid,"n_t2_t4_pairs = %ld\n",n_t2_t4_pairs);
  //fprintf(param_fid,"n_t2_array = %i\n",(int) n_t2_array);
  fprintf(param_fid,"\n"); 
  fclose(param_fid);


  /*allocate arrays*/
  if (LOG_LEVEL>0)
    printf("Initialize Arrays\n");
  x_edges = vector(1,n_bins);
  P_w = vector(1,n_bins);
  bin_traj = lvector(1,nsteps);
  dw_matrix = matrix(1,nprotons,1,n_c2);
  dw_initial = vector(1,nprotons);
  c2 = vector(1,n_c2);
  c2_accum = vector(1,n_c2);
  c3 = matrix(1,n_c3,1,n_c3);
  c3_accum = matrix(1,n_c3,1,n_c3);
  lags = lvector(1,n_c3);
  t2_t4_pairs = matrix(1,n_t2_t4_pairs,1,2);
  //rho_5tensor is really only 4D not fully 5D after I changed how I
  //generate the t2 and t4 times. The first dimension is a singleton
  //dimension.
  rho_5tensor = f5tensor(1,1,1,n_t2_t4_pairs,1,n_bins,1,n_bins,1,n_bins);

    
  /*
   * Setup the time axes for the joint probability density
   */  
  strcpy(string,time_file_name);
  fid=fopen(string,"rt");
  printf("%s\n",string);
  if (fid==NULL) nrerror("Failed to open time_file for reading.");
  for (i=1;i<=n_t2_t4_pairs;i++)
    {
      if(fscanf(fid,"%f %f",&t2_t4_pairs[i][1],&t2_t4_pairs[i][2])!=2)
	nrerror("Failed reading time_file before end of file.");
    }
  fclose(fid);

    /*
     * Read frequency trajectory file
     *
     * The trajectory files should be a command-line input. Read these
     * files line by line and for each time-step calculate the resulting
     * frequency of each OH. Save this in a matrix to be used later for
     * collecting statistics. Finally, make sure to subtract the mean
     * frequency. For a 1 ns trajectory of 2000 protons molecules with
     * zeropadding, the resulting dw_matrix is ~ 1 Gi, so be careful.
     *
     */
    strcpy(name,freq_file_name);
    fid = fopen(name,"rt");
    if(fid==NULL) nrerror("opening file failed.");
    for(j=1;j<=nsteps;j++)
      for(i=1;i<=nprotons;i++)
	  fscanf(fid,"%f",&dw_matrix[i][j]);
    fclose(fid);
    
    //add the zero-pad 
    for(i=1;i<=nprotonsmax;i++)
      for(j=nstepsmax+1;j<n_c2;j++)
        dw_matrix[i][j]=0;
    
    // calculate the mean and subtract it
    sum = 0;
    for(i=1;i<=nprotonsmax;i++)
      for(j=1;j<=nstepsmax;j++)
        sum = sum + dw_matrix[i][j];

    mean = sum/(nstepsmax*nprotonsmax);
    for(i=1;i<=nprotonsmax;i++)
      for(j=1;j<=nstepsmax;j++)
        dw_matrix[i][j] = dw_matrix[i][j] - mean;
    if (LOG_LEVEL>0)
      printf("The mean frequency of the trajectory, counting all protons is %f cm-1\n",mean);
    
    //calculate the standard deviation for diagonstic purposes
    sum = 0;
    for(i=1;i<=nprotonsmax;i++)
    {
      dw = dw_matrix[i];
      for(j=1;j<=nstepsmax;j++)
        sum = sum + dw[j]*dw[j];
    }
    sigma = sqrt(sum/(nstepsmax*nprotonsmax));
    if (LOG_LEVEL>0)
      printf("The standard deviation of the trajectory, counting all protons is %f cm-1\n",sigma);
    
    
    //calculate the skewness for diagonstic purposes
    sum = 0;
    for(i=1;i<=nprotonsmax;i++)
      for(j=1;j<=nstepsmax;j++)
        sum = sum + dw_matrix[i][j]*dw_matrix[i][j]*dw_matrix[i][j];
    skewness = sum/(nstepsmax*nprotonsmax*pow(sigma,3));
    if (LOG_LEVEL>0)
      printf("The skewness of the trajectory, counting all protons is %f.\n",skewness);
    
    /* print these results to param file*/
    strcpy(name,base_name);
    strcat(name,"_parameters_stats.txt");
    param_fid = fopen(name,"at");
    fprintf(param_fid,"# Frequency trajectory statistics:\n");
    fprintf(param_fid,"# Mean     \t%f cm-1\n",mean);
    fprintf(param_fid,"# Std dev  \t%f cm-1\n",sigma);
    fprintf(param_fid,"# Skewness \t%f\n",skewness);
    fclose(param_fid);
    
    
    /*
     *
     * Calculate histogram, 2-, and 3-point correlation functions of the data set
     *
     */
    
    //first calculate a frequency histogram
    printf("Calculate histograms\n");
    for(i=1;i<=nprotonsmax;i++) 
      dw_initial[i]=dw_matrix[i][1]; //intial frequency of each proton
    
    //calculate the "edges" of the bins for the histogram based on these initial values
    histEdges(dw_initial,nprotonsmax,n_bins,x_edges);
    
    //initialize P_w
    for(i=1;i<=n_bins;i++)
      P_w[i]=0;
    
    for(i=1;i<=nprotonsmax;i++)
    {
      //point to each proton's trajectory
      dw = dw_matrix[i];
      
      //accumulate a frequency histogram
      histc(dw,(unsigned long) nstepsmax,x_edges,P_w,n_bins,bin_traj);
    }
    
    /* 
     * Calculate 2-point correlation:
     */
    if (LOG_LEVEL>0)
      printf("Calculate 2-point correlation.\n");
    //initialize c2_accum
    for(i=1;i<=n_c2;i++) c2_accum[i]=0;
    
    //loop over protons taking the autocorrelation of n_c2 points
    for(i=1;i<=nprotonsmax;i++)
    {
      //initialize dw
      dw=dw_matrix[i];
      
      //n_c2 MUST be a power of 2!
       autocorrel(dw,n_c2,c2); 
      for(j=1;j<=n_c2;j++) c2_accum[j] = c2_accum[j] + c2[j];
    }
    //normalize by the number of protons
    for(i=1;i<=n_c2;i++) c2_accum[i]/=nprotonsmax;
    
    
    /* 
     * Calculate 3-point correlations:
     */
    
    //calculate c3 on a grid of "lags" 0:c3_max/dt/n_c3:c3_max/dt
    for(i=1;i<=n_c3;i++)
      {
	lags[i] = (i-1)*c3_lag_stepsize;
      }

    for(j=1;j<=n_c3;j++) 
      for(k=1;k<=n_c3;k++) 
        c3_accum[j][k] = 0.;
        
    if (LOG_LEVEL>0)
      printf("Calculate 3-point correlation.\n");
    for(i=1;i<=nprotonsmax;i++){
      threePointCorrelation(dw_matrix[i],nstepsmax,lags,n_c3,c3);
      for(j=1;j<=n_c3;j++)
	for(k=1;k<=n_c3;k++)
	  c3_accum[j][k] += c3[j][k];
    }

    //normalize c3
    for(j=1;j<=n_c3;j++) 
      for(k=1;k<=n_c3;k++) 
        c3_accum[j][k]/=nprotonsmax;
    
    
    /*
     *    Joint probabilities
     */
    
    // calculate joint probability densities  
    if (LOG_LEVEL>0)
      printf("Calculate joint probablity densities...\n");
    if (rho_5tensor==NULL) nrerror("f5tensor failed");
    for(i=1;i<=1;i++)
      for(j=1;j<=n_t2_t4_pairs;j++)
        for(k=1;k<=n_bins;k++)
          for(l=1;l<=n_bins;l++)
	    for(m=1;m<=n_bins;m++)
              rho_5tensor[i][j][k][l][m]=0;
    
    for(j=1;j<=n_t2_t4_pairs;j++)
      {
	if (LOG_LEVEL>0)
	  printf("...for t2 = %f t4 = %f\n",t2_t4_pairs[j][1],t2_t4_pairs[j][2]);
	for(k=1;k<=nprotonsmax;k++)
	    jointProbabilityDensity3D(dw_matrix[k],nstepsmax,x_edges,n_bins,
				      t2_t4_pairs[j][1],t2_t4_pairs[j][2],
				      rho_5tensor[1][j]);
      }
    
    /*
     * output correlation functions / histograms
     */
    if (LOG_LEVEL>0)
      printf("Output correlation functions and histograms.\n");
    strcpy(name,base_name);
    strcat(name,"_hist.dat");
    fid = fopen(name,"wt");
    if(fid==NULL) nrerror("opening file failed.");
    for(i=1;i<=n_bins;i++) fprintf(fid,"%f %f\n",x_edges[i],P_w[i]);
    fclose(fid);
    
    strcpy(name,base_name);
    strcat(name,"_c2.dat");
    fid = fopen(name,"wt");
    if(fid==NULL) nrerror("opening file failed.");
    for(i=1;i<=n_c2_output;i++)
      fprintf(fid,"%f ",c2_accum[i]);
    fprintf(fid,"\n");
    fclose(fid);
    
    strcpy(name,base_name);
    strcat(name,"_c3.dat");
    fid = fopen(name,"wt");
    if(fid==NULL) nrerror("opening file failed.");
    for(j=1;j<=n_c3;j++)
    {
      for(i=1;i<=n_c3;i++) fprintf(fid,"%f ",c3_accum[j][i]);
      fprintf(fid,"\n");
    }
    fclose(fid);
    
    strcpy(name,base_name);
    strcat(name,"_rho.dat");
    fid = fopen(name,"wt");
    if(fid==NULL) nrerror("opening file failed.");
    for(i=1;i<=1;i++)
      for(j=1;j<=n_t2_t4_pairs;j++)
        for(k=1;k<=n_bins;k++)
          for(l=1;l<=n_bins;l++)
            for(m=1;m<=n_bins;m++)
              fprintf(fid,"%f\n",rho_5tensor[i][j][k][l][m]);
    fclose(fid);
    
    fflush(NULL); //write the files now
    
    if (LOG_LEVEL>0){
      printf("Free memory...\n");
      printf("...dw_matrix\n");
    }
    free_matrix(dw_matrix,1,nprotons,1,nsteps);
    
    if (LOG_LEVEL>0)
      printf("...histogram and c2 and c3 vars\n");
    free_lvector(lags,1,n_c3);
    free_vector(c2,1,n_c2);
    free_matrix(c3,1,n_c3,1,n_c3);
    free_vector(c2_accum,1,n_c2);
    free_matrix(c3_accum,1,n_c3,1,n_c3);
    free_vector(x_edges,1,n_bins);
    
    if (LOG_LEVEL>0)
      printf("...rho_5tensor and t2 arrays\n");
    free_f5tensor(rho_5tensor,1,1,1,n_t2_t4_pairs,1,n_bins,1,n_bins,1,n_bins);
    //free_lvector(t2_array,1,n_t2_array);
    //free_lvector(tdiag_array,1,n_t2_array);
    free_matrix(t2_t4_pairs,1,n_t2_t4_pairs,1,n_t2_t4_pairs);
    
    exit(EXIT_SUCCESS);
}//end main()

/*************************************************************************************/
/*******************************************************************************/
void fft2D(float **P3_re, float **P3_im, int nt)
{
  int it1,it3;
  float tmp;
  
  for(it3=1;it3<=nt;it3++)
    fft(P3_re[it3],P3_im[it3],nt);
  
  for(it3=1;it3<=nt;it3++)
    for(it1=1;it1<it3;it1++)
    {
      tmp=P3_re[it3][it1];
      P3_re[it3][it1]=P3_re[it1][it3];
      P3_re[it1][it3]=tmp;
      tmp=P3_im[it3][it1];
      P3_im[it3][it1]=P3_im[it1][it3];
      P3_im[it1][it3]=tmp;
    }
      for(it3=1;it3<=nt;it3++)
        fft(P3_re[it3],P3_im[it3],nt);
  
  for(it3=1;it3<=nt;it3++)
    for(it1=1;it1<it3;it1++)
    {
      tmp=P3_re[it3][it1];
      P3_re[it3][it1]=P3_re[it1][it3];
      P3_re[it1][it3]=tmp;
      tmp=P3_im[it3][it1];
      P3_im[it3][it1]=P3_im[it1][it3];
      P3_im[it1][it3]=tmp;
    }
}
/*******************************************************************************/
void fft3D(float ***P3_re, float ***P3_im, int nt)

{
  int it1,it3,it5;
  float tmp;
  
  for(it5=1;it5<=nt;it5++)
    for(it3=1;it3<=nt;it3++)
      fft(P3_re[it5][it3],P3_im[it5][it3],nt);
  
  
  for(it5=1;it5<=nt;it5++)
    for(it3=1;it3<=nt;it3++)
      for(it1=1;it1<it3;it1++)
      {
        tmp=P3_re[it5][it3][it1];
        P3_re[it5][it3][it1]=P3_re[it5][it1][it3];
        P3_re[it5][it1][it3]=tmp;
        tmp=P3_im[it5][it3][it1];
        P3_im[it5][it3][it1]=P3_im[it5][it1][it3];
        P3_im[it5][it1][it3]=tmp;
      }
        
        
        for(it5=1;it5<=nt;it5++)
          for(it3=1;it3<=nt;it3++)
            fft(P3_re[it5][it3],P3_im[it5][it3],nt);
  
  for(it5=1;it5<=nt;it5++)
    for(it3=1;it3<=nt;it3++)
      for(it1=1;it1<it3;it1++)
      {
        tmp=P3_re[it5][it3][it1];
        P3_re[it5][it3][it1]=P3_re[it5][it1][it3];
        P3_re[it5][it1][it3]=tmp;
        tmp=P3_im[it5][it3][it1];
        P3_im[it5][it3][it1]=P3_im[it5][it1][it3];
        P3_im[it5][it1][it3]=tmp;
      }
        
        for(it3=1;it3<=nt;it3++)
          for(it5=1;it5<=nt;it5++)
            for(it1=1;it1<it5;it1++)
            {
              tmp=P3_re[it5][it3][it1];
              P3_re[it5][it3][it1]=P3_re[it1][it3][it5];
              P3_re[it1][it3][it5]=tmp;
              tmp=P3_im[it5][it3][it1];
              P3_im[it5][it3][it1]=P3_im[it1][it3][it5];
              P3_im[it1][it3][it5]=tmp;
            }
              
              
              for(it5=1;it5<=nt;it5++)
                for(it3=1;it3<=nt;it3++)
                  fft(P3_re[it5][it3],P3_im[it5][it3],nt);
  
  for(it3=1;it3<=nt;it3++)
    for(it5=1;it5<=nt;it5++)
      for(it1=1;it1<it5;it1++)
      {
        tmp=P3_re[it5][it3][it1];
        P3_re[it5][it3][it1]=P3_re[it1][it3][it5];
        P3_re[it1][it3][it5]=tmp;
        tmp=P3_im[it5][it3][it1];
        P3_im[it5][it3][it1]=P3_im[it1][it3][it5];
        P3_im[it1][it3][it5]=tmp;
      }
        
}
/*****************************************************************************/
void fft(float *P_re,float *P_im, int n)
{
  int i,j,ir,nfill=16;
  float *a;
  
  a=vector(1,2*n*nfill);
  for (i=1;i<n;i++)
  {
    a[nfill*2*(i-1)+1]=P_re[i];
    a[nfill*2*(i-1)+2]=P_im[i];
    for(j=1;j<nfill;j++)
    {
	  a[nfill*2*(i-1)+2*j+1]=P_re[i]+(P_re[i+1]-P_re[i])*j/nfill;
	  a[nfill*2*(i-1)+2*j+2]=P_im[i]+(P_im[i+1]-P_im[i])*j/nfill;
    }
  }
  a[nfill*2*(n-1)+1]=P_re[n];
  a[nfill*2*(n-1)+2]=P_im[n];
  for(j=1;j<nfill;j++)
  {
    a[nfill*2*(n-1)+2*j+1]=P_re[n]-P_re[n]*j/nfill;
    a[nfill*2*(n-1)+2*j+2]=P_im[n]-P_im[n]*j/nfill;
  }
  
  
  four1(a,n*nfill,+1);
  for (i=1;i<=n;i++)
  {
    ir=i+(nfill-1)*n+n/2;
    if(i>n/2)ir-=n*nfill;
    P_re[i]=a[2*ir-1];
    P_im[i]=a[2*ir];
  }
  
  
  free_vector(a,1,2*n);
}
/*************************************************************************************/

int sign(float x)
{
  if(x>0) return(1);
  if(x<0) return(-1);
  return(0);
}

/********************Numerical Recipiers Routines****************************/

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(float data[], unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software ,2:. */




/*****************************************************************************/
#define NR_END 1
#ifndef FREE_ARG 
#define FREE_ARG char*
#endif

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  printf("Numerical Recipes run-time error...\n");
  printf("%s\n",error_text);
  printf("...now exiting to system...\n");
  exit(EXIT_FAILURE);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;
  
  v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  if (!v) nrerror("allocation failure in vector()");
  return (v-nl+NR_END);
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;
  
  /* allocate pointers to rows */
  m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return (m);
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;
  
  /* allocate pointers to pointers to rows */
  t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;
  
  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  
  /* return pointer to array of pointers to rows */
  return t;
}

float *****f5tensor(long n1l,long n1h,long n2l,long n2h,long n3l,long n3h,long n4l,long n4h,long n5l,long n5h)
//allocate a float 5Dtensor with range t[n1l,n1h][n2l,n2h][n3l,n3h][n4l,n4h][n5l,n5h]
{
  long i,j,k,l,n1=n1h-n1l+1,n2=n2h-n2l+1,n3=n3h-n3l+1,n4=n4h-n4l+1,n5=n5h-n5l+1;
  float *****t;
  
  //allocate pointers to pointers to pointers to pointers (!!) to rows
  t=(float *****)malloc((size_t)((n1+NR_END)*sizeof(float ****)));
  if(!t) nrerror("allocation failure 1 in f5tensor()" );
  t +=NR_END;
  t -=n1l;
  
  //allocate pointers to pointers to pointers to rows
  t[n1l]=(float ****)malloc((size_t)((n1*n2+NR_END)*sizeof(float ***)));
  if(!t[n1l]) nrerror("allocation failure 2 in f5tensor()");
  t[n1l] +=NR_END;
  t[n1l] -=n2l;
  
  //allocate pointers to pointers to rows
  t[n1l][n2l]=(float ***) malloc((size_t)((n1*n2*n3+NR_END)*sizeof(float **)));
  if (!t[n1l][n2l]) nrerror("allocation failure 3 in f5tensor()");
  t[n1l][n2l] += NR_END;
  t[n1l][n2l] -= n3l;
  
  /* allocate pointers to rows and set pointers to them */
  t[n1l][n2l][n3l]=(float **) malloc((size_t)((n1*n2*n3*n4+NR_END)*sizeof(float*)));
  if (!t[n1l][n2l][n3l]) nrerror("allocation failure 4 in f5tensor()");
  t[n1l][n2l][n3l] += NR_END;
  t[n1l][n2l][n3l] -= n4l;
  
  /* allocate rows and set pointers to them */
  t[n1l][n2l][n3l][n4l]=(float *) malloc((size_t)((n1*n2*n3*n4*n5+NR_END)*sizeof(float)));
  if (!t[n1l][n2l][n3l][n4l]) nrerror("allocation failure 5 in f5tensor()");
  t[n1l][n2l][n3l][n4l] += NR_END;
  t[n1l][n2l][n3l][n4l] -= n5l;
  
  for(l=n4l+1;l<=n4h;l++) t[n1l][n2l][n3l][l]=t[n1l][n2l][n3l][l-1]+n5;
  for(k=n3l+1;k<=n3h;k++){
    t[n1l][n2l][k]=t[n1l][n2l][k-1]+n4;
    t[n1l][n2l][k][n4l]=t[n1l][n2l][k-1][n4l]+n4*n5;
    for(l=n4l+1;l<=n4h;l++) t[n1l][n2l][k][l]= t[n1l][n2l][k][l-1]+n5;
  }
  for(j=n2l+1;j<=n2h;j++){
    t[n1l][j]=t[n1l][j-1]+n3;
    t[n1l][j][n3l]=t[n1l][j-1][n3l]+n3*n4;
    t[n1l][j][n3l][n4l]=t[n1l][j-1][n3l][n4l]+n3*n4*n5;
    for(k=n3l+1;k<=n3h;k++){
      t[n1l][j][k]= t[n1l][j][k-1]+n4;
      t[n1l][j][k][n4l]=t[n1l][j][k-1][n4l]+n4*n5;
      for(l=n4l+1;l<=n4h;l++){
        t[n1l][j][k][l]=t[n1l][j][k][l-1]+n5;
      }
      for(l=n4l+1;l<=n4h;l++) t[n1l][j][n3l][l]=t[n1l][j][n3l][l-1]+n5;
    }
  }
  for(i=n1l+1;i<=n1h;i++){
    t[i]=t[i-1]+n2;
    t[i][n2l]=t[i-1][n2l]+n2*n3;
    t[i][n2l][n3l]=t[i-1][n2l][n3l]+n2*n3*n4;
    t[i][n2l][n3l][n4l]=t[i-1][n2l][n3l][n4l]+n2*n3*n4*n5;
    for(j=n2l+1;j<=n2h;j++){
      t[i][j]=t[i][j-1]+n3;
      t[i][j][n3l]=t[i][j-1][n3l]+n3*n4;
      t[i][j][n3l][n4l]=t[i][j-1][n3l][n4l]+n3*n4*n5;
      for(k=n3l+1;k<=n3h;k++){
        t[i][j][k]=t[i][j][k-1]+n4;
        t[i][j][k][n4l]=t[i][j][k-1][n4l]+n4*n5;
        for(l=n4l+1;l<=n4h;l++) t[i][j][k][l]=t[i][j][k][l-1]+n5;
        for(l=n4l+1;l<=n4h;l++) t[i][j][n3l][l]=t[i][j][n3l][l-1]+n5;
      }
    }
    for(k=n3l+1;k<=n3h;k++){
      t[i][n2l][k]=t[i][n2l][k-1]+n4;
      t[i][n2l][k][n4l]=t[i][n2l][k-1][n4l]+n4*n5;
      for(l=n4l+1;l<=n4h;l++) t[i][n2l][k][l]=t[i][n2l][k][l-1]+n5;
    }
    for(l=n4l+1;l<=n4h;l++) t[i][n2l][n3l][l]=t[i][n2l][n3l][l-1]+n5;
  }  
  return t;
}

void free_f5tensor(float *****t,long n1l,long n1h,long n2l,long n2h,long n3l,long n3h,long n4l,long n4h,long n5l,long n5h)
//free a tensor allocated with f5tensor()
{
  free((FREE_ARG) (t[n1l][n2l][n3l][n4l]+n5l-NR_END));
  free((FREE_ARG) (t[n1l][n2l][n3l]+n4l-NR_END));
  free((FREE_ARG) (t[n1l][n2l]+n3l-NR_END));
  free((FREE_ARG) (t[n1l]+n2l-NR_END));
  free((FREE_ARG) (t+n1l-NR_END));
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
  nh++;
  free((char*) (v+nl-NR_END));
}


void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
  nch++;nrh++;
  free((char*) (m[nrl]+ncl-NR_END));
  free((char*) (m+nrl-NR_END));
}


int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;
  
  /* allocate pointers to rows */
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  
  /* allocate rows and set pointers to them */
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;
  
  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((char*)(v+nl-NR_END));
}


void free_imatrix(int **m,long nrl,long nrh,long ncl,long nch)
/* free an int matrix allocated by imatrix() */
{
  nch++;nrh++;
  free((char*) (m[nrl]+ncl-NR_END));
  free((char*) (m+nrl-NR_END));
}


void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
  free((char*) (t[nrl][ncl]+ndl-NR_END));
  free((char*) (t[nrl]+ncl-NR_END));
  free((char*) (t+nrl-NR_END));
}


/************************************************************************************/


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


void twofft(float data1[], float data2[], float fft1[], float fft2[], 
            unsigned long n) 
/*Given two real input arrays data1[1..n] and data2[1..n], this
routine calls four1and returns two complex output arrays,
fft1[1..2n] and fft2[1..2n], each of complex length n (i.e., real
                                                       length2*n), which contain the discrete Fourier transforms of the
respective data arrays. n MUST be an integer power of 2. */
{ 
  void four1(float data[], unsigned long nn, int isign); 
  unsigned long nn3,nn2,jj,j; 
  float rep,rem,aip,aim; 
  nn3=1+(nn2=2+n+n); 
  for (j=1,jj=2;j<=n;j++,jj+=2) { //Pack the two real arrays into one
                                  //complexarray.
    fft1[jj-1]=data1[j]; 
    fft1[jj]=data2[j]; 
  } 
  four1(fft1,n,1); //Transform the complex array. 
  fft2[1]=fft1[2]; 
  fft1[2]=fft2[2]=0.0; 
  for (j=3;j<=n+1;j+=2) { 
    rep=0.5*(fft1[j]+fft1[nn2-j]); //Use symmetries to separate the
                                   //two transforms.
    rem=0.5*(fft1[j]-fft1[nn2-j]); 
    aip=0.5*(fft1[j+1]+fft1[nn3-j]); 
    aim=0.5*(fft1[j+1]-fft1[nn3-j]); 
    fft1[j]=rep; //Ship them out in two complexarrays. 
    fft1[j+1]=aim; 
    fft1[nn2-j]=rep; 
    fft1[nn3-j] = -aim; 
    fft2[j]=aip; 
    fft2[j+1] =-rem; 
    fft2[nn2-j]=aip; 
    fft2[nn3-j]=rem; 
  } 
}

void realft(float data[], unsigned long n, int isign) 
/*Calculates the Fourier transform of a set of n real-valued data
points. Replaces this data (which is stored in array data[1..n]) by
the positive frequency half of its complex Fourier transform. The
real-valued first and last components of the complex transform are
returned as elements data[1] and data[2], respectively. n must be a
power of 2. This routine also calculates the inverse transform of
acomplex data array if it is the transform of real data. (Result in
                                                          this case must be multiplied by 2/n.) */
{ 
  void four1(float data[], unsigned long nn, int isign); 
  unsigned long i,i1,i2,i3,i4,np3; 
  float c1=0.5,c2,h1r,h1i,h2r,h2i; 
  double wr,wi,wpr,wpi,wtemp,theta; //Double precision for the
                                    //trigonometric recurrences. 
  theta=3.141592653589793/(double) (n>>1); //Initialize the recurrence. 
  if (isign == 1){ 
    c2 =-0.5; 
    four1(data,n>>1,1); //The forward transform is here. 
  }else { 
    c2=0.5; //Otherwise set up for an inverse transform.
    theta = -theta; 
  } 
  wtemp=sin(0.5*theta); 
  wpr = -2.0*wtemp*wtemp; 
  wpi=sin(theta); 
  wr=1.0+wpr; 
  wi=wpi; 
  np3=n+3; 
  for (i=2;i<=(n>>2);i++) {// Case i=1 done separately below. 
    i4=1+(i3=np3-(i2=1+(i1=i+i-1))); 
    h1r=c1*(data[i1]+data[i3]); //The two separate transforms are
                                //separated out of data.
    h1i=c1*(data[i2]-data[i4]); 
    h2r =-c2*(data[i2]+data[i4]); 
    h2i=c2*(data[i1]-data[i3]); 
    data[i1]=h1r+wr*h2r-wi*h2i;// Here they are recombined to form the
                               // true transform of the original real data. 
      data[i2]=h1i+wr*h2i+wi*h2r; 
      data[i3]=h1r-wr*h2r+wi*h2i; 
      data[i4] = -h1i+wr*h2i+wi*h2r; 
      wr=(wtemp=wr)*wpr-wi*wpi+wr; //The recurrence. 
      wi=wi*wpr+wtemp*wpi+wi; 
  } 
  if (isign == 1){
    data[1] = (h1r=data[1])+data[2]; //Squeeze the first and last data
                                     //together to get them all within
                                     //the original array. 
    data[2] = h1r-data[2]; 
  }else { 
    data[1]=c1*((h1r=data[1])+data[2]); 
    data[2]=c1*(h1r-data[2]); 
    four1(data,n>>1,-1); //This is the inverse transform for the
                         //case isign = -1.
  } 
} 

void correl(float data1[], float data2[], unsigned long n, float ans[]) 
/*Computes the correlation of two real data sets data1[1..n] and data2[1..n]
(including any user-supplied zeropadding). n MUST be aninteger power
of two. The answer is returned as the first n points in ans[1..2*n]
stored in wrap-aroundorder, i.e., correlations at increasingly
negative lags are in ans[n] on down to ans[n/2+1], while correlations
at increasingly positive lags are in ans[1] (zerolag) on up to
ans[n/2]. Note that ans must be supplied in the calling program with
length at least 2*n, since it is also used as workingspace. Sign
convention of this routine: if data1 lags data2, i.e., is shifted to
the right of it, then ans will show a peak at positive lags. */
{ 
  void realft(float data[], unsigned long n, int isign); 
  void twofft(float data1[], float data2[], float fft1[], float fft2[], 
              unsigned long n); 
  unsigned long no2,i; 
  float dum,*fft; 
  fft=vector(1,n<<1); 
  twofft(data1,data2,fft,ans,n); //Transform both data vectors at once. 
  no2=n>>1; //Normalization for inverseFFT. 
  for (i=2;i<=n+2;i+=2) { 
    //Multiply to find FFT of their correlation. 
    ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2; 
    ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2; 
  } 
  ans[2]=ans[n+1]; //Pack first and last into one element. 
  realft(ans,n,-1); //Inverse transform gives correlation. 
  free_vector(fft,1,n<<1); 
} 

void autocorrel(float data[], unsigned long n, float ans[]) 
/*Computes the autocorrelation of one real data set data[1..n]
(including any user-supplied zeropadding). n MUST be aninteger power
of two. The answer is returned ans[1..n] (note difference with correl)
stored in wrap-aroundorder, i.e., correlations at increasingly
negative lags are in ans[n] on down to ans[n/2+1], while correlations
at increasingly positive lags are in ans[1] (zerolag) on up to
ans[n/2]. */
{ 
  void realft(float data[], unsigned long n, int isign); 
  unsigned long no2,i; 
  float *data2;
  
  //protect data because f*C@!g realft overwrites its inputs
  data2 = vector(1,n);
  for(i=1;i<=n;i++) //optimize with memcopy?
    data2[i]=data[i];
  
  realft(data2,n,1);
  no2=n>>1; //Normalization for inverse FFT. 
  for (i=2;i<=n;i+=2) { 
    ans[i-1]=(data2[i-1]*data2[i-1]+data2[i]*data2[i])/no2; 
    ans[i]=0;
  } 
  ans[2]=data2[2]*data2[2]/no2;
  realft(ans,n,-1); //Inverse transform gives correlation. 
  free_vector(data2,1,n); 
} 

void histc(float seq[],unsigned long n,float x_edges[],float h[],unsigned long n_edges,unsigned long bin[])
/*take a sequence seq[1..n] and histogram it into bins given by
x_edges returning the histogram in h[1..nedges] and the bin number of
each point in the sequence as bin[1..n]. N.B. h is not reinitialized
by the routine, so it can be used cumulatively (and you had better
initialize it the first time you call histc). Values of seq outside
the range of x_edges are folded into the first and last bins,
respectively. It is probably best to ensure that x_edges contains all
the data.*/
{
  long i;
  float min_x,max_x,width,norm;
  min_x = x_edges[1];
  max_x = x_edges[n_edges];
  width = max_x - min_x;
  norm=(n_edges-1)/width;
  
  for(i=1;i<=n;i++)
  {
    bin[i] = floor( (seq[i]-min_x)*norm ) + 1 ;
    if (bin[i]<1) bin[i]=1;
    if (bin[i]>n_edges) bin[i]=n_edges;
    //printf("%f\t%d\t%f\n",seq[i],(int) bin[i],x_edges[bin[i]]);
    
    h[bin[i]]++;
  }
}

#ifndef NR_END
#define NR_END 1
#endif

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;
  
  v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
  if (!v) nrerror("allocation failure in lvector()");
  return v-nl+NR_END;
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}


void histEdges(float seq[],unsigned long n_steps,unsigned long n_edges,float x[])
{
  float max_x,min_x,range;
  unsigned long i;
  
  max_x=seq[1];
  min_x=seq[1];
  for(i=2;i<=n_steps;i++)
  {
    if(seq[i]>max_x) max_x = seq[i];
    if(seq[i]<min_x) min_x = seq[i];
  }
  
  range = max_x - min_x;
  
  //add a little cushion at the extremes
  max_x = max_x + 0.1*range;
  min_x = min_x - 0.1*range;
  range = max_x - min_x;
  
  //initialize x and y
  for(i=1;i<=n_edges;i++)
  {
    x[i] = min_x + (i-1)*range/(n_edges-1);
  }
}

float ** hist(float seq[],unsigned long n_steps,unsigned long n_edges)
/*calculate the histogram of a sequence seq[1..n_steps] divided into
n_bins, returning the x_edges and histogram values packed into
histostruct[1..2][1..n_bins]. x_edges = histostruct[1][..],
histogram = [2][1..nbins]*/
{
  float **histostruct,*x,*y;
  unsigned long i,*dummy;
  
  dummy = lvector(1,n_steps);
  
  histostruct = matrix(1,2,1,n_edges);
  //point x and y at the rows of histostruct
  x = histostruct[1];
  y = histostruct[2];
  
  histEdges(seq,n_steps,n_edges,x);
  
  //initialize the histogram to 0
  for(i=1;i<=n_edges;i++)
  {
    y[i] = 0;
  }
  
  histc(seq,n_steps,x,y,n_edges,dummy);
  
  free_lvector(dummy,1,n_steps);
  return histostruct;
}

void jointProbabilityDensity3D(float seq[],unsigned long n_steps,float edges[],unsigned long n_edges,unsigned long lag1,unsigned long lag2,float ***rho_3)
{
  unsigned long n_steps_reduced;
  
  n_steps_reduced = n_steps-lag1-lag2;
  //  printf("jointProbabilityDensity3D calling histTri\n");
  histTri(&seq[1],&seq[lag1+1],&seq[lag1+lag2+1],n_steps_reduced,edges,n_edges,rho_3);
  
}
void histTri(float x[],float y[],float z[],unsigned long n_steps,float edges[],unsigned long n_edges,float ***hist_tensor)
{
  unsigned long *xbin,*ybin,*zbin,n_edges2,n_edges3,i;
  float *pf,*dummy,*xyz,*xyz_edges;
  
  n_edges2 = n_edges*n_edges;
  n_edges3 = n_edges*n_edges2;
  
  //allocate vars
  dummy = vector(1,n_edges);
  xbin  = lvector(1,n_steps);
  ybin  = lvector(1,n_steps);
  zbin  = lvector(1,n_steps);
  
  //  printf("histTri calling histc (1)\n");
  //calculate independent histograms (stored in dummy), but it is the
  //trajectory part that we want (xbin, ybin, zbin)
  histc(x,n_steps,edges,dummy,n_edges,xbin);
  //  printf("histTri calling histc (2)\n");
  histc(y,n_steps,edges,dummy,n_edges,ybin);
  //  printf("histTri calling histc (3)\n");
  histc(z,n_steps,edges,dummy,n_edges,zbin);
  
  
  //loop over steps and calculate the 1D index into the 3D histogram
  //tensor from xbin, ybin and zbin
  //  printf("allocate xyz\n");
  xyz = vector(1,n_steps);
  xyz_edges = vector(1,n_edges3);
  
  for(i=1;i<=n_steps;i++)
    xyz[i] = (float) zbin[i] + (ybin[i]-1)*n_edges + (xbin[i]-1)*n_edges2;
  
  for(i=1;i<=n_edges3;i++)
    xyz_edges[i] = (float) i;
  
  //  printf("set pointer\n");
  pf = &hist_tensor[1][1][1]-1; //point to the data block of hist_tensor (offset so it is h[1..n^3])
  
  //  printf("histTri calling histc (4)\n");
  histc(xyz,n_steps,xyz_edges,pf,n_edges3,xbin); //xbin now serves as
                                                 //a dummy var
  
  //free memory
  free_vector(dummy,1,n_edges);
  free_lvector(xbin,1,n_steps);
  free_lvector(ybin,1,n_steps);
  free_lvector(zbin,1,n_steps);
  free_vector(xyz,1,n_steps);
  free_vector(xyz_edges,1,n_edges);
}

void threePointCorrelation(float seq[],unsigned long n_steps,unsigned long lags[],unsigned long n_lags,float **c3)
{
  //float **c3;
  unsigned long i,j,k,t1,t2,n_steps_reduced;
  
  //c3 = matrix(1,n_lags,1,n_lags);
  
  for(i=1;i<=n_lags;i++)
  {
    t1=lags[i];
    for(j=1;j<=n_lags;j++)
	{
	  t2=lags[j];
	  
	  c3[i][j]=0;
	  n_steps_reduced = n_steps - t1 - t2;
	  //loop over sequence calculating the mean
	  for(k=1;k<=n_steps_reduced;k++)
	    c3[i][j]+=seq[k]*seq[k+t1]*seq[k+t1+t2];
	  c3[i][j]/=n_steps_reduced;
	}
  } 
  //return c3;
}
void normalizeResults(int nt,unsigned long isample,
                      float S_re[],float S_im[],
                      float **P1_re,float **P1_im,float **P2_re,float **P2_im,
                      float ***R1_re,float ***R1_im,float ***R2_re,float ***R2_im,
                      float ***R3_re,float ***R3_im,float ***R4_re,float ***R4_im,
                      float Sf_re[],float Sf_im[],
                      float **P1f_re,float **P1f_im,float **P2f_re,float **P2f_im,
                      float ***R1f_re,float ***R1f_im,float ***R2f_re,float ***R2f_im,
                      float ***R3f_re,float ***R3f_im,float ***R4f_re,float ***R4f_im)
{
  unsigned long it1,it3,it5;
  
  for(it1=1;it1<=nt;it1++)
  {
    Sf_re[it1]=S_re[it1]/isample;
    Sf_im[it1]=S_im[it1]/isample;
    for(it3=1;it3<=nt;it3++)
	{
	  P1f_re[it3][it1]=P1_re[it3][it1]/isample;
	  P1f_im[it3][it1]=P1_im[it3][it1]/isample;
	  P2f_re[it3][it1]=P2_re[it3][it1]/isample;
	  P2f_im[it3][it1]=P2_im[it3][it1]/isample;
	  for(it5=1;it5<=nt;it5++)
      {
        R1f_re[it5][it3][it1]=R1_re[it5][it3][it1]/isample;
        R1f_im[it5][it3][it1]=R1_im[it5][it3][it1]/isample;
        R2f_re[it5][it3][it1]=R2_re[it5][it3][it1]/isample;
        R2f_im[it5][it3][it1]=R2_im[it5][it3][it1]/isample;
        R3f_re[it5][it3][it1]=R3_re[it5][it3][it1]/isample;
        R3f_im[it5][it3][it1]=R3_im[it5][it3][it1]/isample;
        R4f_re[it5][it3][it1]=R4_re[it5][it3][it1]/isample;
        R4f_im[it5][it3][it1]=R4_im[it5][it3][it1]/isample;
      }
	}
  }
}

void fourierTransformResults(int nt,
                             float Sf_re[],float Sf_im[],
                             float **P1f_re,float **P1f_im,float **P2f_re,float **P2f_im,
                             float ***R1f_re,float ***R1f_im,float ***R2f_re,float ***R2f_im,
                             float ***R3f_re,float ***R3f_im,float ***R4f_re,float ***R4f_im)
{
  /*Calculate 1D-Fouriertrafo*/
  printf("Calculate 1D-Fouriertransform\n");
  fft(Sf_re,Sf_im,nt);
  /*Calculate 2D-Fouriertrafo*/
  printf("Calculate 2D-Fouriertransform\n");
  fft2D(P1f_re,P1f_im,nt);
  fft2D(P2f_re,P2f_im,nt);
  /*Calculate 3D-Fouriertrafo*/
  printf("Calculate 3D-Fouriertransform\n");
  fft3D(R1f_re,R1f_im,nt);
  fft3D(R2f_re,R2f_im,nt);
  fft3D(R3f_re,R3f_im,nt);
  fft3D(R4f_re,R4f_im,nt);
}

void writeResults(char base_name[],int nt,char extension[],
                  float Sf_re[],float Sf_im[],
                  float **P1f_re,float **P1f_im,float **P2f_re,float **P2f_im,
                  float ***R1f_re,float ***R1f_im,float ***R2f_re,float ***R2f_im,
                  float ***R3f_re,float ***R3f_im,float ***R4f_re,float ***R4f_im)
{
  char name[100];
  unsigned long it1,it3,it5;
  FILE *out1D,*out2D,*out3D;
  
  /*write results*/
  printf("Write Results\n");
  strcpy(name,base_name);
  strcat(name,"_spec1D");
  strcat(name,extension);
  out1D=fopen(name,"wt");
  //	    for(it1=ncut*nt/8+2;it1<=(8-ncut)*nt/8;it1++)
  for(it1=1;it1<=nt;it1++)
  {
    fprintf(out1D,"%8.5f ",Sf_re[it1]);
    fprintf(out1D,"\n");
  }
  fclose(out1D);
  
  strcpy(name,base_name);
  strcat(name,"_spec2D");
  strcat(name,extension);
  out2D=fopen(name,"wt");
  //	    for(it1=ncut*nt/8+2;it1<=(8-ncut)*nt/8;it1++)
  //	      for(it3=ncut*nt/8+2;it3<=(8-ncut)*nt/8;it3++)
  for(it1=1;it1<=nt;it1++)
    for(it3=1;it3<=nt;it3++)
    {
      //		  fprintf(out2D,"%8.5f %8.5f",P1f_re[it3][it1],P2f_re[it3][nt+2-it1]);
      fprintf(out2D,"%8.5f %8.5f",P1f_re[it3][it1],P2f_re[it3][nt+1-it1]);
      fprintf(out2D,"\n");
    }
      fclose(out2D);
  
  strcpy(name,base_name);
  strcat(name,"_spec3D");
  strcat(name,extension);
  out3D=fopen(name,"wt");
  //	    for(it1=ncut*nt/8+2;it1<=(8-ncut)*nt/8;it1++)
  //	      for(it3=ncut*nt/8+2;it3<=(8-ncut)*nt/8;it3++)
  //		for(it5=ncut*nt/8+2;it5<=(8-ncut)*nt/8;it5++)
  for(it1=1;it1<=nt;it1++)
    for(it3=1;it3<=nt;it3++)
      for(it5=1;it5<=nt;it5++)
      {
        //		    fprintf(out3D,"%8.5f %8.5f %8.5f %8.5f",R1f_re[it5][it3][it1],			    R2f_re[it5][it3][nt+2-it1],R3f_re[it5][nt+2-it3][it1],R4f_re[it5][nt+2-it3][nt+2-it1]);
        fprintf(out3D,"%8.5f %8.5f %8.5f %8.5f",
                R1f_re[it5][it3][it1],
                R2f_re[it5][it3][nt+1-it1],
                R3f_re[it5][nt+1-it3][it1],
                R4f_re[it5][nt+1-it3][nt+1-it1]);
        fprintf(out3D,"\n");
      }
        fclose(out3D);
}

float **submatrix(float **a,long oldrl,long oldrh,long oldcl, long oldch,long newrl,long newcl)
{
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  float **m;
  
  //allocate array of pointers to rows
  m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("allocation failure in submatrix()");
  m += NR_END;
  m -= newrl;
  
  //set pointers to rows
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
  
  //return pointer to array of pointers to rows
  return m;
}

void free_submatrix(float **b,long nrl,long nrh,long ncl,long nch)
{
  free((FREE_ARG) (b+nrl-NR_END));
}

