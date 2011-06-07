#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <fnmatch.h>

/*
 *  7.2.2007 S. Garrett-Roe to read an MD trajectory (one file of coordinates and 
 *  one file of forces) and calculate a frequnecy trajectory.
 *
 *  26 May 2011 start to rewrite with better input of options, the
 *  possibility to use other mappings of coords and forces to
 *  frequency, the possibility to output dipole and anharmonicity
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

#define PI 3.1415927
#define wavenumbersToInvPs 2.99792458e-2
int DEBUG_LEVEL = 0;

void display_usage( void )
{
  puts("mdTrajToResponse - convert a gromacs md trajectory into a frequency and dipole trajectory.\n" 
       "\n"
       "Usage:\n"
       "mdTrajToFreq [options] -p input_parameter_file -f forces.xvg -c coords.xvg -o base_name ...\n"
       "\n"
       "Options:\n"
       "-p parameter file \n"
       "-c coords file \n"
       "-f forces file \n"
       "-o output base name"
       "-v -verbose\n"
       "\tThe verbosity level (debug information). Multiple v's (e.g. -vv) increase the verbosity level\n"
       "\n"
       "Input file format:\n"
       "Files should be plain text. Comments lines start with a #.\n"
       "Subsequent lines should have one parameter assignment per line.\n"
       "Allowed paramters are:\n"
       " dt step size from the gromacs trajectory (in ps)\n"
       " order the highest order of response function to calculate (1, 3, or 5)\n"
       " fit_order the level of taylor expansion from field to frequency and/or dipole moment. 0 is just a constant. 1 is linear (most parameterizations). 2 is quadratic, etc. This must be assigned before the parameters are passed. Parameters a are for frequency, b for x_{ij}, and mu_mug is the ratio of condensed phase ot gas phase dipole moment. These two-dimensional matrices are zero-based. a and b have size [0..(order-1)/2)][0..fit_order]. mu_mug has size [0..fit_order] Subscript 0 is for v = 0-1, 1 is for v = 1-2, 2 is for v = 2-3.\n"
       " so the \\omega_{01} = a[0][0] + a[0][1] E + a[0][2] E^2 + ... similarly for b and mu."
       " These are assigned in the parameter file as a_i_j  = val. (e.g. a_0_1 = 0.1 => a[0][1]= 0.1).\n"
       "q_H charge on the proton (in units of e) to convert force to field\n"
       "\n"
       "note on units! To convert from Skinner's paramters into gromacs units you divide by a factor of 49614. So the OD (SPC/E) w_10 from Schmidt J Chem Phys 123, 044513 (2005) is\n"
"\tw_10 = 2792.7  + 7559.5 / q_H / 49614 Force_H\n"
"where q_H is the proton charge (0.4238e for SPC/E) and F_H is the force on the proton (in kJ mol^{-1} nm^{-1}) projected onto the OD axis.\n"
       "flag_massweightedforces = {0 | 1}. Skinner type parameters this is 0 (false). For geissler type parameters this is 1 (true).\n"
       "flag_noncondon 0 (false) use condon approx and do not output dipole trajectory. 1 (true) calculate dipole moment and output it as a trajectory also.\n"
       "flag_twolevelsystem 1 (true) only calculate 01 transitions. 0 (false) calculate n_levels = (order+1)/2 excited state signals (eg order 3 -> n_levels = 2 v=0-1 and v=1-2)\n."
       "--compress_output use gz on output files\n"
       "\n");
  exit( EXIT_FAILURE );
}

//the short command line options, the letter is followed by a : if the
// option takes an argument
static const char *optString = "c:f:o:p:vh?";

static const struct option longOpts[] = {
  // all the long options go here with the letter that they map do (or
  //  0 if they have no short version)
  //  { "render",no_argument, NULL, 0},
  { "coords", required_argument,NULL,'c'},
  { "forces", required_argument,NULL,'f'},
  { "output", required_argument,NULL,'o'},
  { "parameters", required_argument,NULL,'p'},
  { "compress_output",no_argument,NULL,0},
  { "verbose",no_argument,NULL,'v'},
  { "help",no_argument,NULL, 'h'},
  { NULL, no_argument, NULL,0}
};

struct globalArgs_t {
  //the parameters that need to be globally visible can be stord in
  //this data structure
  int test; 
  float dt;
  int order;
  int fit_order; //order of expansion
  float **a; //zero based matrix size (0,(order-1)/2,0,fit_order)
  float **b; //zero based matrix size (0,(order-1)/2,0,fit_order)
  float *mu_mug; 
  float q_H; //charge on the proton (to convert force to field)
  int flag_massweightedforces;
  int flag_noncondon;
  int flag_twolevelsystem;
  int flag_compressoutput;
} globalArgs;

void initialize_globalArgs( void ){
  // generate the default values of the global arguments above
  globalArgs.test = 1;
  globalArgs.flag_compressoutput = 0;
}

void display_options( void ){
  // show the user what the current global arguments (options) are
  time_t my_time=time(0); // time process started
  int i,j;

  printf("\n");
  printf("Global options (globalArgs)\n");
  printf("\n");

  printf("test \t%d\n",globalArgs.test); 
  printf("dt \t%f ps\n",globalArgs.dt); 
  for (i=0;i<=(globalArgs.order-1)/2;i++){
    for (j=0;j<=globalArgs.fit_order;j++){
      printf("a_%i_%i = %f\n",i,j,globalArgs.a[i][j]);
    }
  }
  for (i=0;i<=(globalArgs.order-1)/2;i++){
    for (j=0;j<=globalArgs.fit_order;j++){
      printf("b_%i_%i = %f\n",i,j,globalArgs.b[i][j]);
    }
  }
  for (j=0;j<=globalArgs.fit_order;j++){
    printf("mu_mug_%i = \t%f \n",j,globalArgs.mu_mug[j]); 
  }
  
  printf("flag_massweightedforces \t%d\n",globalArgs.flag_massweightedforces); 
  printf("flag_noncondon \t%d\n",globalArgs.flag_noncondon); 
  printf("flag_twolevelsystem \t%d\n",globalArgs.flag_twolevelsystem); 
  printf("flag_compressoutput \t%d\n",globalArgs.flag_compressoutput); 

  printf("\n\n");
}


void mdTrajToFreq(const char* coord_file_name, const char* force_file_name, const char* base_name)
{
  FILE *coord_fid,*force_fid;
  FILE *fid; 
  FILE **w_fid_array;// = malloc(3*sizeof(FILE*));//TEST TEST
  FILE **x_fid_array;// = malloc(3*sizeof(FILE*));//TEST TEST!!!
  char *name,*string;
  //,coord_file_name[100],force_file_name[100],base_name[100],
  
  const float NA = 6.023e23; // mol^-1
  const float mass = 1.67262158e-27; //kg
  const float mH = 1.; //a.u.
  const float mO = 16.; //a.u.
  const float c = 2.99792458e10; // cm/s
  const float h = 6.626e-34; // Js

  float mu; //reduced mass
  
  //vars for force trajectory
  float dummy;
  unsigned long nsteps, natoms, nmols, nprotons;
  
  float *posO,*posH1,*posH2,*forceO,*forceH1,*forceH2,*force,*bond,bond_length;
  unsigned long i,j,k;
  float w; //frequency
  float x; //position matrix element
  float force_proj; //forces projected on OH bond
  float field;
  int i_level;
  float dipole;

  //int dummy;

  const int flag_massweightedforces = globalArgs.flag_massweightedforces;
  const int flag_noncondon = globalArgs.flag_noncondon;
  const int flag_twolevelsystem = globalArgs.flag_twolevelsystem;
  //  const float **a = globalArgs.a;
  //  const float **b = globalArgs.b;
  float **a, **b; //why doesn't the above work???
  const float *mu_mug = globalArgs.mu_mug;
  const int n_levels = (globalArgs.flag_twolevelsystem == 0 ? (globalArgs.order+1)/2 : 1); //TEST TEST!!!
  const int fit_order = globalArgs.fit_order;
  const int flag_compressoutput = globalArgs.flag_compressoutput;
  const float q_H = globalArgs.q_H;


  a = globalArgs.a;
  b = globalArgs.b;

  if (DEBUG_LEVEL>=1){
    printf("entering mdTrajtoFreq()\n");
    printf("nlevels = %d\n",n_levels);
    for (i=0;i<=(globalArgs.order-1)/2;i++){
      for (j=0;j<=globalArgs.fit_order;j++){
	printf("g.a_%i_%i = %f\ta_%i_%i = %f\n",(int)i,(int)j,globalArgs.a[i][j],(int)i,(int)j,a[i][j]);
      }
    }
    for (i=0;i<=(globalArgs.order-1)/2;i++){
      for (j=0;j<=globalArgs.fit_order;j++){
	printf("g.b_%i_%i = %f\tb_%i_%i = %f\n",(int)i,(int)j,globalArgs.b[i][j],(int)i,(int)j,b[i][j]);
      }
    }
    for (j=0;j<=globalArgs.fit_order;j++){
      printf("g.mu_mug_%i = \t%f\tmu_mu_g_%i = \t%f\n",(int)j,globalArgs.mu_mug[j],(int)j,mu_mug[j]); 
    }
  }//end debug printf's

  mu = mH*mO/(mH + mO);
  //fact = (float) sqrt(9./4.*(Anh/w0)*1e18*1e6/pow(NA,2)/(w0*c*h)/(mu*mass)*1e-24);
  //fact = sqrt(9/4*(Anh/w0)*10^18*10^6/NA^2/(w0*c*h)/(mu*mass)*10^-24*wavenumbersToInvPs);
  
  // use the shell command "wc -l" to scan the length of the coord
  // file to get the number of steps
  // and read one line with "head" and then "wc -w" to get the
  // number of molecules
  printf("Determine number of steps and number of molecules from coordinate file.\n");
  if (asprintf(&string,"wc -l %s",coord_file_name)<0){
    fprintf(stderr,"failed to write string");
    exit(EXIT_FAILURE);
  }
  fid = popen(string,"r");
  free(string);
  if(fid==NULL) nrerror("Opening 'wc -l coord_file_name' pipe failed.");
  fscanf(fid,"%ld",&nsteps);
  pclose(fid);
  printf("Found %ld steps in file %s\n", nsteps, coord_file_name);
  
  if (asprintf(&string,"head -n 1 %s | wc -w",coord_file_name) < 0){
    fprintf(stderr,"failed to write string");
    exit(EXIT_FAILURE);
  }
  fid = popen(string,"r");
  free(string);
  if(fid==NULL) nrerror("Opening 'head -n 1 coord_file_name | wc -w' pipe failed.");
  fscanf(fid,"%ld",&natoms);
  pclose(fid);
  
  natoms = (natoms-1)/3; //substract 1 for the time stamp 
  //and div by 3 for x y z coords
  nmols = natoms/3; //3 atoms per H2O molecule
  nprotons = nmols*2;//2 protons
  printf("Found %ld atoms, %ld molecules in file %s\n", natoms, nmols, coord_file_name);
  
  posO    = vector(1,3);
  posH1   = vector(1,3);
  posH2   = vector(1,3);
  forceO  = vector(1,3);
  forceH1 = vector(1,3);
  forceH2 = vector(1,3);
  bond    = vector(1,3);
  force   = vector(1,3);
    
  /*
   * Read trajectory files (coordinate and forces)
   *
   * The trajectory files should be a command-line input. Read these
   * files line by line and for each time-step calculate the resulting
   * frequency of each OH. Save this in a matrix to be used later for
   * collecting statistics. Finally, make sure to subtract the mean
   * frequency. For a 1 ns trajectory of 2000 protons molecules with
   * zeropadding, the resulting dw_matrix is ~ 1 Gi, so be careful.
   *
   */
  
  
  // open coordinate and force files
  coord_fid = fopen(coord_file_name,"rt");
  if(coord_fid==NULL) nrerror("opening coordinate file failed.");
  force_fid = fopen(force_file_name,"rt");
  if(force_fid==NULL) nrerror("opening force file failed.");
  printf("Open trajectory files\n");  
  
  //set these up as /dev/NULL if not desired and then just calculate everything
  w_fid_array = malloc(n_levels*sizeof(FILE*));
  x_fid_array = malloc(n_levels*sizeof(FILE*));

  //open output files
  if (flag_compressoutput==0){
    //open text files
    if (DEBUG_LEVEL>=0) printf("open files for text output\n");
    for (i_level=0;i_level<n_levels;i_level++){
      //freq files
      //sprintf(name,"%s_dw%i.dat",base_name,(int)i_level);
      if (asprintf(&name,"%s_dw%i.dat",base_name,(int)i_level) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
      w_fid_array[i_level]=fopen(name,"wt"); 
      if (w_fid_array[i_level]==NULL){
	fprintf(stderr,"Error opening file %s for output\n",name);
	exit(EXIT_FAILURE);
      }
      free(name);

      //mu files
      if (flag_noncondon==1){
  	    //sprintf(name,"%s_mu%i.dat",base_name,(int)i_level);
	if (asprintf(&name,"%s_mu%i.dat",base_name,(int)i_level) < 0)
	  {
	    fprintf(stderr,"failed to write string");
	    exit(EXIT_FAILURE);
	  }
      } else {
	    //sprintf(name,"/dev/null");
	if (asprintf(&name,"/dev/null") < 0)
	  {
	    fprintf(stderr,"failed to write string");
	    exit(EXIT_FAILURE);
	  }
      }      
      x_fid_array[i_level]=fopen(name,"wt"); 
      if (x_fid_array[i_level]==NULL){
	fprintf(stderr,"Error opening file %s for output\n",name);
	exit(EXIT_FAILURE);
      }
      free(name);
    }
  } else {
    //open a pipe to gz to compress the output
    if (DEBUG_LEVEL>=0) printf("open pipe to gzip files for compressed output\n");
    for (i_level=0;i_level<n_levels;i_level++){
      //freq files
      //    sprintf(name,"gzip > %s_dw%i.dat.gz",base_name,(int)i_level);
      if (asprintf(&name,"gzip > %s_dw%i.dat.gz",base_name,(int)i_level) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
      if (DEBUG_LEVEL>=2) printf("%s\n",name);
      w_fid_array[i_level] = popen(name,"w"); //TEST TEST
      if (w_fid_array[i_level]==NULL){
	fprintf(stderr,"Error opening pipe %s for output\n",name);
	exit(EXIT_FAILURE);
      }
      free(name);

      //mu files
      if (flag_noncondon==1){
	//sprintf(name,"gzip > %s_mu%i.dat.gz",base_name,(int)i_level);
	if (asprintf(&name,"gzip > %s_mu%i.dat.gz",base_name,(int)i_level) < 0)
	  {
	    fprintf(stderr,"failed to write string");
	    exit(EXIT_FAILURE);
	  }
      } else {
	//    sprintf(name,"/dev/null");
	if (asprintf(&name,"/dev/null") < 0)
	  {
	    fprintf(stderr,"failed to write string");
	    exit(EXIT_FAILURE);
	  }
      }        
      x_fid_array[i_level]=popen(name,"w"); //TEST TEST
      if (x_fid_array[i_level]==NULL){
	fprintf(stderr,"Error opening pipe %s for output\n",name);
	exit(EXIT_FAILURE);
      }
      free(name);

    }
  }

  if(DEBUG_LEVEL>=3)
    for (i_level=0;i_level<n_levels;i_level++){
      printf("%p\n",w_fid_array[i_level]);
      printf("%p\n",x_fid_array[i_level]);
    }

  /*
   * loop over timesteps 
   */
  for(i=1;i<=nsteps;i++)
    {
      // the first number is the time, which I throw away
      // if reading the file doesn't correctly read a number then nrerror
      if(fscanf(coord_fid,"%f",&dummy) !=1)
        nrerror("Reading coord file failed before end of file (time).");
      if(fscanf(force_fid,"%f",&dummy)!=1)
        nrerror("Reading force file failed before end of file (time).");
      for(j=1;j<=nmols;j++)
	{
	  // read configuration for one time-step
	  if(fscanf(coord_fid,"%f %f %f",&posO[1],&posO[2],&posO[3])!=3)
	    nrerror("Reading coord file failed before end of file (O atom).");
	  if(fscanf(coord_fid,"%f %f %f",&posH1[1],&posH1[2],&posH1[3])!=3)
	    nrerror("Reading coord file failed before end of file (H1 atom).");
	  if(fscanf(coord_fid,"%f %f %f",&posH2[1],&posH2[2],&posH2[3])!=3)
	    nrerror("Reading coord file failed before end of file (H2 atom).");
	  
	  // read in forces for one time-step
	  if(fscanf(force_fid,"%f %f %f",&forceO[1],&forceO[2],&forceO[3])!=3)
	    nrerror("Reading force file failed before end of file (O atom).");
	  if(fscanf(force_fid,"%f %f %f",&forceH1[1],&forceH1[2],&forceH1[3])!=3)
	    nrerror("Reading force file failed before end of file (H1 atom).");
	  if(fscanf(force_fid,"%f %f %f",&forceH2[1],&forceH2[2],&forceH2[3])!=3)
	    nrerror("Reading force file failed before end of file (H2 atom).");
        
        

	  // calculate the bond vector and bond length for H1
	  bond_length=0;
	  for(k=1;k<=3;k++)
	    bond[k] = posH1[k] - posO[k];
	  bond_length=sqrt(bond[1]*bond[1] + bond[2]*bond[2] + bond[3]*bond[3]);
	  
	  if (flag_massweightedforces==1){
	    // calculate mass weighted forces
	    for(k=1;k<=3;k++)
	      force[k] = (forceH1[k]/mH - forceO[k]/mO)*mu;
	  } else {
	    // just use the proton force
	    for(k=1;k<=3;k++)
	      force[k] = forceH1[k];
	  }
       
	  //calculate projection of forces on H1
	  force_proj = 0;
	  for(k=1;k<=3;k++)
	    force_proj += force[k]*bond[k];
	  field = force_proj/q_H;


	  //calculate frequencies and anharmonicities
	  for (i_level=0;i_level<n_levels;i_level++){
	    // scale field to effective OH stretch frequency
	    w = 0;
	    for (k=0;k<=fit_order;k++)
	      //w+=0.001;
	      w += a[i_level][k]*pow(field,k); //TEST TEST 
	    
	    //dw_matrix[2*j-1][i] *= fact/bond_length/FUDGE;
	    
	    // non-condon effects here
	    x = 0;
	    //calculate the matrix elements here
	    for (k=0;k<=fit_order;k++)
	      x += b[i_level][k]*pow(w,k); //TEST TEST 
	    //calculate the harmonic dipole ???
	    dipole = 0;
	    for (k=0;k<=fit_order;k++)
	      dipole = mu_mug[k]*pow(field,k);
	    //combine
	    dipole = dipole*x;

	    //output freq 
	    if(DEBUG_LEVEL>=3){
	      printf("print w = %f step %lu mol %lu i_level %i pointer %p\n",w,i,j,i_level,w_fid_array[i_level]);
	      printf("%p\n",w_fid_array[i_level]);
	      printf("%p\n",x_fid_array[i_level]);
	    }
	    fprintf(w_fid_array[i_level],"%12.5f ",w);
	    if (DEBUG_LEVEL>=3) fflush(w_fid_array[i_level]);

		  
	    //output dipole
	    if (DEBUG_LEVEL>=3) printf("print x = %f step %lu mol %lu\n",dipole,i,j);
	    fprintf(x_fid_array[i_level],"%12.5f ",dipole);
	    if (DEBUG_LEVEL>=3) fflush(x_fid_array[i_level]);
	    
	  }//end loop over levels
	
	
	}//end for j=1:nmols
      
      //print a newline to each freq and mu file
      for (i_level=0;i_level<n_levels;i_level++){
	fprintf(w_fid_array[i_level],"\n");
	fprintf(x_fid_array[i_level],"\n");
      }

    } //end for i = 1:nsteps

  // close coordinate and force files
  fclose(coord_fid);
  fclose(force_fid);
  
  // close output files
  for (i_level=0;i_level<n_levels;i_level++){
    fclose(w_fid_array[i_level]);
    fclose(x_fid_array[i_level]);
  }
  
    
  if (DEBUG_LEVEL>=1) printf("clean up memory\n");
  if (DEBUG_LEVEL>=1) printf("...positions and forces\n");
  free_vector(posO,1,3);
  free_vector(posH1,1,3);
  free_vector(posH2,1,3);
  free_vector(forceO,1,3);
  free_vector(forceH1,1,3);
  free_vector(forceH2,1,3);
  free_vector(bond,1,3);
  free_vector(force,1,3);
  
  // free fid_arrays
  free_matrix(a,0,n_levels,0,fit_order);
  free_matrix(b,0,n_levels,0,fit_order);
  free(w_fid_array);
  free(x_fid_array);

} //end mdTrajToFreq()

void read_input_parameters(char *parameter_file_name){
  //read the input paramter file
  FILE *fid;
  char *unparsed_line;
  char *line;
  char name[100];
  //char *name;
  //int bytes_read;
  float val;
  int count=0;
  size_t len;
  int i=0,j=0; //indices into matrices

  printf("read parameters\n");

  if(access(parameter_file_name,R_OK)==-1){
    fprintf(stderr,"Unable to read %s\n",parameter_file_name);
    exit(EXIT_FAILURE);
  }
  fid = fopen(parameter_file_name,"rt");
  while(feof(fid)==0){
    count += 1;

    
    unparsed_line = fgetln(fid,&len);
    if (DEBUG_LEVEL>=3) printf("len  %i\n",(int)len);
    line = calloc(len,sizeof(char));
    strncpy(line,unparsed_line,(int)len);
    if (DEBUG_LEVEL>=2)
      printf("line %i: %s",count,line);
    
    /* //getline is not in FreeBSD's libc yet ... 
    line = (char *) malloc(100 * sizeof(char));
    bytes_read= getline(&line, &len,fid);
    if (bytes_read == -1){
      fprintf(stderr,"Unable to read line of parameter file:\n%s",line);
      exit(EXIT_FAILURE);
    }
    */
    if (len == 0){
      continue;
    } //end if len>0

    if(strncmp(line,"#",1)==0) {
      // if the first character is a comment character do nothing on this loop
      if (DEBUG_LEVEL>=1) printf("comment: %s",line);
      continue;
    }

    // if it is not a comment character try to process it
    //if (sscanf(line,"%as = %f",&name,&val) < 2){ //I don't know why this doesn't work!!! (also uncomment name declaration above
    if (sscanf(line,"%s = %f",name,&val) < 2){
      printf("skip line: \n%s",line);
      continue;
    }
    
    if (DEBUG_LEVEL>=1) printf("assignment: %s = %f\n",name,val);
    
    //now that we have name and value try to do assignments if we
    //recognize the name
    //step size
    if (strcmp("dt",name) == 0){
      globalArgs.dt = val;
    }
    //order of response functions to calculate
    if (strcmp("order",name) ==0){
      globalArgs.order = (int) val;
    }
    //order of taylor expasion of field to frequency calculations
    if (strcmp("fit_order",name) ==0){
      if (val>1){
	fprintf(stderr,"You requested fit_order=%i, but fit order >1 not yet implemented! \n",(int) val);
	exit(EXIT_FAILURE);
      }
      
      globalArgs.fit_order = (int) val;
      //initialize the parameter arrays
      globalArgs.a = matrix(0,(globalArgs.order-1)/2,0,globalArgs.fit_order);
      globalArgs.b = matrix(0,(globalArgs.order-1)/2,0,globalArgs.fit_order);
      globalArgs.mu_mug = vector(1,globalArgs.fit_order);
    }
    //parameters for frequency \omega_{ji}
    if (fnmatch("a_*_*",name,FNM_CASEFOLD)==0){
      sscanf(name,"a_%d_%d",&i,&j);
      if ( DEBUG_LEVEL >= 1 ) printf("assigning a_%d_%d = %f\n",i,j,val);
      if ( i <= ( globalArgs.order - 1 ) / 2 ) 
	if ( j <= globalArgs.fit_order )
	  globalArgs.a[i][j] = val;
    }
    //parameters for coordinate x_{ji}
    if (fnmatch("b_*_*",name,FNM_CASEFOLD)==0){
      sscanf(name,"b_%d_%d",&i,&j);
      if (DEBUG_LEVEL >=1) printf("assigning b_%d_%d = %f\n",i,j,val);
      if ( i <= ( globalArgs.order - 1 ) / 2 ) 
	if ( j <= globalArgs.fit_order )
	  globalArgs.b[i][j] = val;
    }
    //parameters for dipole moment
    if (fnmatch("mu_mug_*",name,FNM_CASEFOLD)==0){
      sscanf(name,"mu_mug_%d",&j);
      if (DEBUG_LEVEL >=1) printf("assigning mu_mug_%d = %f\n",j,val);
      if ( j <= globalArgs.fit_order )
	globalArgs.mu_mug[j] = val;
    }
    if (fnmatch("q_H",name,FNM_CASEFOLD)==0){
      globalArgs.q_H = val;
    }
    if (fnmatch("flag_twolevelsystem",name,FNM_CASEFOLD)==0){
      globalArgs.flag_twolevelsystem = (int) val;
    }
    if (fnmatch("flag_noncondon",name,FNM_CASEFOLD)==0){
      globalArgs.flag_noncondon = (int) val;
    }
    if (fnmatch("flag_massweightedforces",name,FNM_CASEFOLD)==0){
      globalArgs.flag_massweightedforces = (int) val;
    }
    if (fnmatch("flag_compressoutput",name,FNM_CASEFOLD)==0){
      globalArgs.flag_compressoutput= (int) val;
    }
    
  } //end while(feof(fid)==0)

  fclose(fid);

}

int main( int argc, char *argv[] ) {
  //This is the main function, the glue that sticks the pieces
  //together

  //declare main() variables
  int opt = 0, longIndex=0; //for getopt_long
  char *parameter_file_name,*coord_file_name,*force_file_name,*base_name;

  /* initialize global arguments */
  initialize_globalArgs();

  printf("process options\n");
  /* process command line arguments */
  while((opt = getopt_long( argc, argv, optString,longOpts, &longIndex))!=-1)
    switch(opt){
    case 'c': //coord file
      if (asprintf(&coord_file_name,"%s",optarg) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
      break;
    case 'f': //forces file
      if (asprintf(&force_file_name,"%s",optarg) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
      break;
    case 'o': //output base name
      if (asprintf(&base_name,"%s",optarg) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
      break;
    case 'p': //parameter file
      //strcpy(parameter_file_name,optarg);
      if (asprintf(&parameter_file_name,"%s",optarg) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
      break;
    case 't': //test
      sscanf(optarg,"%i", &globalArgs.test);
      break;
    case 'v': /* increase verbosity */
      DEBUG_LEVEL+=1;
      break;
    case 'h': /* fall-through is intentional */
    case '?':
      display_usage();
      break;
    case 0: /* all options that don't have a short one character version */
      if (fnmatch("compress_output",longOpts[longIndex].name, FNM_CASEFOLD ) == 0 ){
	globalArgs.flag_compressoutput=1;
      }
/*       if( strcmp( "link_radius", longOpts[longIndex].name ) == 0 ){ */
/* 	sscanf(optarg,"%f",&val); */
/* 	globalArgs.link_radius = val; */
/*       } else if( strcmp( "link_minimum_weight", longOpts[longIndex].name ) == 0 ){ */
/* 	sscanf(optarg,"%f",&val); */
/* 	globalArgs.link_minimum_weight = val;      */
/*       } else { */
/* 	fprintf(stderr,"WARNING!!! Commandline option %s not recognized!!! \n",string); */
/* 	exit(EXIT_FAILURE); */
/*       } */
      break;
    default:
      /* should never get here */
      fprintf(stderr,"WARNING!!! Commandline option %s not recognized!!! \n",longOpts[longIndex].name);
      exit(EXIT_FAILURE);
      break;
    }

  if(argc-optind!=0)
    {
      display_usage();
    }

  // input parameters
  read_input_parameters(parameter_file_name);

  // make sure options are okay
  display_options();

  //input required data
  //read_input_data(parameter_file_name);

  //do the main calculation
  mdTrajToFreq(coord_file_name,force_file_name,base_name);

  //output the results

  //cleanup

  // we're done!
  exit(EXIT_SUCCESS);

}//end main()

/*************************************************************************************
 *************************************************************************************
 * 
 *
 *    Math helper functions go after here
 *
 *
 *************************************************************************************
 *************************************************************************************/
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
x_edges returning the histogram in h[1..nedges] and the bin number
of each point in the sequence as bin[1..n]. N.B. h is not
reinitialized by the routine, so it can be used cumulatively (and
                                                              you had better initialize it the first time you call histc). Values
of seq outside the range of x_edges are folded into the first and
last bins, respectively. It is probably best to ensure that x_edges
contains all the data.*/
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
	  
	  //c3[i][j]=0;
	  n_steps_reduced = n_steps - t1 - t2;
	  //loop over sequence calculating the mean
	  for(k=1;k<=n_steps_reduced;k++)
	    c3[i][j]+=seq[k]*seq[k+t1]*seq[k+t2];
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

