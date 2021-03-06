#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

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




#define PI 3.1415927
#define order 5
#define wavenumbersToInvPs 2.99792458e-2
//#define NMOLS 1019 //number of waters in the box
//#define NPROTONS 2*NMOLS
//#define NSTEPS 10001 //1 ns trajectory (now  a variable determined at runtime)

int main(int iarg, char **input) 
{
  int nt=32;                  /* Number of time Points for t1, t3, and
  * t5 -- must be a power of 2!  */
  int ntint=6;                /* number of time-steps to pre-integrate */
  float dt=0.01; //ps, time-steps output by GROMACS
  int nskip=15; //number of steps to skip between trajectory samples
                 //(0.150 ps with 0.010 ps steps for H2O)
  //  int npoptimes = 3;          /* number of population times to evaluate (obsolete) */
  
  time_t my_time=time(0); // time process started
  
  long idum;
  FILE *param_fid,*fid;
  char extension[100],numb[100],name[100],freq_file_name[100],time_file_name[100],base_name[100],string[100];
    
  //vars for force trajectory
  unsigned long nsteps=-1, natoms, nmols, nprotons=-1;
  
  double sum;
  float **dw_matrix;
  float *dw; //a pointer into dw_matrix representing the frequency
             //trajectory of one proton
  float *dw_initial; //the initial frequency of each proton
  float mean,sigma,skewness; //characterize the frequency traj
    
  //vars for time arrays (soon obsolete)
  //  unsigned long *t2_array,*tdiag_array; //arrays of time steps used
  //  const unsigned long n_t2_array = 3; //number of time steps
    
  //vars for response function calcultions
  float **pdw1,**pdw2,*sdw,***rdw1,***rdw2,***rdw3,***rdw4;
  float ***R1_re,***R1_im,***R2_re,***R2_im,***R3_re,***R3_im,***R4_re,***R4_im;
  float **P1_re,**P2_re,**P2_im,**P1_im,*S_re,*S_im;
  float ***R1f_re,***R1f_im,***R2f_re,***R2f_im,***R3f_re,***R3f_im,***R4f_re,***R4f_im,
    **P1f_re,**P2f_re,**P2f_im,**P1f_im,*Sf_re,*Sf_im;
  unsigned long n_t2_t4_pairs; //number of times we will calculate joint prob dens
  float **t2_t4_pairs;//matrix of times (t2,t4)
  float t2,t4;
  unsigned long i,j,isample,it,it1,it3,it5,nt2,nt4,ntraject,nsamples,ipoptime,iproton;
  float *dwint;           /* trajectories after preintegrating/coarse graining  */   

  //change this logic to allow restart optional arguments

  opterr = 0;
  int c; /* for getopt */

  printf("This is %s called as\n",input[0]);
  for (i=0;i<iarg;i++)
    printf("%s ",input[i]);
  printf("\n");
  printf("On: %s\n",ctime(&my_time));
  printf("\n"); 

  /* read the command line parameters and set values as appropriate */
  while ((c = getopt(iarg,input,"g:i:n:p:s:t:"))!=-1)
    switch(c){
    case 'h': /* help line and exit */
      printf("USAGE:\n   %s [-g nt] [-i ntint] [-n nsteps] [-p nprotons] [-s nskip] [-t dt] freq_file time_file output_file_prefix\n\n",input[0]);
      exit(EXIT_FAILURE);
    case 'g': /* set nt number of grid points */
      if(sscanf(optarg,"%i",&nt)!=1) nrerror("failed reading commandline option");
      break;
    case 'i': /* set ntint number of (pre)integration steps */
      if(sscanf(optarg,"%i",&ntint)!=1) nrerror("failed reading commandline option");
      break;
    case 'n': /* set nsteps number of steps in the md simulation */
      if(sscanf(optarg,"%lu",&nsteps)!=1) nrerror("failed reading commandline option");
      break;
    case 'p': /* set nproton number of protons in trajectory */
      if(sscanf(optarg,"%lu",&nprotons)!=1) nrerror("failed reading commandline option");
      break;
    case 's': /* set nskip number of steps to skip to make independent samples */
      if(sscanf(optarg,"%i",&nskip)!=1) nrerror("failed reading commandline option");
      break;
    case 't': /* set dt size of base time-step from md simulation */
      if(sscanf(optarg,"%f",&dt)!=1) nrerror("failed reading commandline option");
      break;
    case '?':
      nrerror("unknown command line option.");
      /*      if (isprint(optopt))
	fprintf(stderr,"Unknown option `-%c'\n",optopt);
      else
	fprintf(stderr,"Unknown option character `\\x%x'.\n",optopt);
	exit(EXIT_FAILURE);*/
    default:
      abort();
    }
      
    if(iarg-optind<3)
    {
      printf("Too few input arguments. Please give 3 file names.\n");
      printf("USAGE:\n   %s [-g nt] [-i ntint] [-n nsteps] [-p nprotons] [-s nskip] [-t dt] freq_file time_file output_file_prefix\n\n",input[0]);
      exit(EXIT_FAILURE);
    }

    strcpy(freq_file_name,input[0+optind]);
    strcpy(time_file_name,input[1+optind]);
    strcpy(base_name,input[2+optind]);
    
    /* if not specified on the comand line... */
    if (nsteps==-1){
      //use the shell command "wc -l" to scan the length of the coord file to get the number of steps
      // and read one line with "head" and then "wc -w" to get the number of molecules
      printf("Determine number of steps and number of molecules from frequency file.\n");
      strcpy(string,"wc -l ");
      strcat(string,freq_file_name);
      fid = popen(string,"r");
      if(fid==NULL) nrerror("Opening 'wc -l freq_file_name' pipe failed.");
      fscanf(fid,"%ld",&nsteps);
      pclose(fid);
      printf("Found %ld steps in file %s\n", nsteps, freq_file_name);
    }
    if (nprotons==-1){
      strcpy(string,"head -n 1 ");
      strcat(string,freq_file_name);
      strcat(string," | wc -w");
      fid = popen(string,"r");
      if(fid==NULL) nrerror("Opening 'head -n 1 freq_file_name | wc -w' pipe failed.");
      fscanf(fid,"%ld",&nprotons);
      pclose(fid);
      nmols = nprotons/2;
      natoms = nmols*3;
      printf("Found %ld protons in file %s\n", nprotons, freq_file_name);
    }
    else{
      nmols = nprotons/2;
      natoms = nmols*3;
    }      
    
    /*
     * Setup the time axes for the joint probability density
     */
    //set up the time axis that will be used
    //note: because of the way the response functions are calculated,
    // it is better to have t2 AND t2/2 be multiples of ntint
    /*
      t2_array = lvector(1,n_t2_array);
    tdiag_array = lvector(1,n_t2_array);
    
    t2_array[1] = 0;
    t2_array[2] = 2*ntint; //2*6*dt = 0.120 ps
    t2_array[3] = 18*ntint; //18*6*dt = 1.08 ps
    for(i=1;i<=n_t2_array;i++)
      tdiag_array[i] = t2_array[i]/2;
    */
    strcpy(string,"wc -l ");
    strcat(string,time_file_name);
    fid=popen(string,"r");
    if(fid==NULL) nrerror("Failed opening pipe wc -l time_file_name");
    fscanf(fid,"%ld",&n_t2_t4_pairs);
    pclose(fid);
    t2_t4_pairs = matrix(1,n_t2_t4_pairs,1,2);

    strcpy(string,time_file_name);
    fid=fopen(string,"rt");
    if (fid==NULL) nrerror("Failed to open time_file for reading.");
    for (i=1;i<=n_t2_t4_pairs;i++)
      {
	if(fscanf(fid,"%f %f",&t2_t4_pairs[i][1],&t2_t4_pairs[i][2])!=2)
	  nrerror("Failed reading time_file before end of file.");
      }
    fclose(fid);

    /* write parameter file */
    printf("Writing parameter file\n");
    strcpy(name,base_name);
    strcat(name,"_parameters_R5.txt");
    param_fid = fopen(name,"wt");
    fprintf(param_fid,"# This is %s called as\n",input[0]);
    for (i=0;i<iarg;i++)
      fprintf(param_fid,"%s ",input[i]);
    fprintf(param_fid,"\n");
    fprintf(param_fid,"# On: %s\n",ctime(&my_time));
    fprintf(param_fid,"\n"); 
    fprintf(param_fid,"# Parameters:\n");
    fprintf(param_fid,"dt = %f\n",dt);
    fprintf(param_fid,"nskip = %i\n", nskip);
    fprintf(param_fid,"nsteps = %i\n", (int) nsteps);
    fprintf(param_fid,"nmols = %i\n", (int) nmols);
    //    fprintf(param_fid,"n_t2_array = %i\n",(int) n_t2_array);
    fprintf(param_fid,"nt = %i\n",(int) nt);
    fprintf(param_fid,"ntint = %i\n",(int) ntint);
    fprintf(param_fid,"n_t2_t4_pairs = %ld\n",n_t2_t4_pairs);
    fprintf(param_fid,"\n"); 
    
    /*initialize random generator*/
    idum=-10;
    
    /*allocate arrays*/
    printf("Initialize Arrays\n");
    
    R1_re=f3tensor(1,nt,1,nt,1,nt);
    R1_im=f3tensor(1,nt,1,nt,1,nt);
    R2_re=f3tensor(1,nt,1,nt,1,nt);
    R2_im=f3tensor(1,nt,1,nt,1,nt);
    R3_re=f3tensor(1,nt,1,nt,1,nt);
    R3_im=f3tensor(1,nt,1,nt,1,nt);
    R4_re=f3tensor(1,nt,1,nt,1,nt);
    R4_im=f3tensor(1,nt,1,nt,1,nt);
    P1_re=matrix(1,nt,1,nt);
    P1_im=matrix(1,nt,1,nt);
    P2_re=matrix(1,nt,1,nt);
    P2_im=matrix(1,nt,1,nt);
    S_re=vector(1,nt);
    S_im=vector(1,nt);
    
    R1f_re=f3tensor(1,nt,1,nt,1,nt);
    R1f_im=f3tensor(1,nt,1,nt,1,nt);
    R2f_re=f3tensor(1,nt,1,nt,1,nt);
    R2f_im=f3tensor(1,nt,1,nt,1,nt);
    R3f_re=f3tensor(1,nt,1,nt,1,nt);
    R3f_im=f3tensor(1,nt,1,nt,1,nt);
    R4f_re=f3tensor(1,nt,1,nt,1,nt);
    R4f_im=f3tensor(1,nt,1,nt,1,nt);
    P1f_re=matrix(1,nt,1,nt);
    P1f_im=matrix(1,nt,1,nt);
    P2f_re=matrix(1,nt,1,nt);
    P2f_im=matrix(1,nt,1,nt);
    Sf_re=vector(1,nt);
    Sf_im=vector(1,nt);
    
    
    sdw=vector(1,nt);
    pdw1=matrix(1,nt,1,nt);
    pdw2=matrix(1,nt,1,nt);
    rdw1=f3tensor(1,nt,1,nt,1,nt);
    rdw2=f3tensor(1,nt,1,nt,1,nt);
    rdw3=f3tensor(1,nt,1,nt,1,nt);
    rdw4=f3tensor(1,nt,1,nt,1,nt);
    
    dw_matrix = matrix(1,nprotons,1,nsteps);
    dw_initial = vector(1,nprotons);

    /*
     * Read frequency file
     *
     * The frequency file can be generated from a gromacs trajectory
     * by first exporting force and coordinate files, then running
     * mdTrajToFreq.c on these files.  The frequency of each OH is in
     * a matrix for collecting statistics. Finally, make sure to
     * subtract the mean frequency. For a 1 ns trajectory of 2000
     * protons molecules with zeropadding, the resulting dw_matrix is
     * ~ 1 Gi, so be careful.
     *
     */
      strcpy(name,freq_file_name);
      fid = fopen(name,"rt");
      strcpy(string,"opening ");
      strcat(string,name);
      strcat(string," failed.");
      if(fid==NULL) nrerror(string);
      for(j=1;j<=nsteps;j++)
	for(i=1;i<=nprotons;i++)
	  {
	    fscanf(fid,"%f",&dw_matrix[i][j]);
	  }
      fclose(fid);    
        
    // calculate the mean and subtract it
    sum = 0;
    for(i=1;i<=nprotons;i++)
    {
      //dw=dw_matrix[i];
      for(j=1;j<=nsteps;j++)
        //sum = sum + dw[j];
        sum = sum + dw_matrix[i][j];
    }
    mean = sum/(nsteps*nprotons);
    for(i=1;i<=nprotons;i++)
      for(j=1;j<=nsteps;j++)
        dw_matrix[i][j] = dw_matrix[i][j] - mean;
    printf("The mean frequency of the trajectory, counting all protons is %f cm-1\n",mean);
    
    //calculate the standard deviation for diagonstic purposes
    sum = 0;
    for(i=1;i<=nprotons;i++)
    {
      dw = dw_matrix[i];
      for(j=1;j<=nsteps;j++)
        sum = sum + dw[j]*dw[j];
    }
    sigma = sqrt(sum/(nsteps*nprotons));
    printf("The standard deviation of the trajectory, counting all protons is %f cm-1\n",sigma);
    
    
    //calculate the skewness for diagonstic purposes
    sum = 0;
    for(i=1;i<=nprotons;i++)
      for(j=1;j<=nsteps;j++)
        sum = sum + dw_matrix[i][j]*dw_matrix[i][j]*dw_matrix[i][j];
    skewness = sum/(nsteps*nprotons*pow(sigma,3));
    printf("The skewness of the trajectory, counting all protons is %f.\n",skewness);
    
    /* print these results to param file*/
    fprintf(param_fid,"# Frequency trajectory statistics:\n");
    fprintf(param_fid,"# Mean     \t%f cm-1\n",mean);
    fprintf(param_fid,"# Std dev  \t%f cm-1\n",sigma);
    fprintf(param_fid,"# Skewness \t%f\n",skewness);
    fclose(param_fid);
        
    
    /*
     *
     * Sample the trajectories for non-linear polarization 
     *
     */

    //if flag_restart then pick use the saved ipoptime

    //loop over population times calculating t2=t,t4=0, and t2=t/2,t4=t/2 
    for(ipoptime=1;ipoptime<=n_t2_t4_pairs;ipoptime++)
    {
      t2 = t2_t4_pairs[ipoptime][1];
      t4 = t2_t4_pairs[ipoptime][2];
      printf("Time delay %li of %li, t2 = %f t4 = %f\n",ipoptime,n_t2_t4_pairs,t2,t4);
      
      /*initialize arrays*/
      for(it1=1;it1<=nt;it1++)
      {
        S_re[it1]=0;
        S_im[it1]=0;
        for(it3=1;it3<=nt;it3++)
        {
	      P1_re[it3][it1]=0;
	      P1_im[it3][it1]=0;
	      P2_re[it3][it1]=0;
	      P2_im[it3][it1]=0;
	      for(it5=1;it5<=nt;it5++)
          {
            R1_re[it5][it3][it1]=0;
            R1_im[it5][it3][it1]=0;
            R2_re[it5][it3][it1]=0;
            R2_im[it5][it3][it1]=0;
            R3_re[it5][it3][it1]=0;
            R3_im[it5][it3][it1]=0;
            R4_re[it5][it3][it1]=0;
            R4_im[it5][it3][it1]=0;
          }
        }
      }
      //or use values from the restart file
      
      strcpy(extension,"_t2_");
      sprintf(numb,"%d",(int)(t2));
      strcat(extension,numb);
      strcat(extension,"_t4_");
      sprintf(numb,"%d",(int)(t4));
      strcat(extension,numb);
      strcat(extension,".dat");
      printf("Extension: %s\n",extension);
      
      
      //nt2 and nt4 are the number of coarse grained points
      //during population times in trajectory
      nt4=round(t4/ntint);
      nt2=round(t2/ntint);
      
      
     
 ntraject=3*nt+nt2+nt4; // number of points in the coarse-grained trajectory
      dwint=vector(1,ntraject); // that selection averaged over ntint
                                // points (for H20 about 6)
      
      //calculate the maximum number of dw trajectories in the total
      //time evolution as <nsamples>
      if (nsteps<ntraject*ntint)
      {
        nrerror("ERROR: trajectory is too short for response function calculation!");
      }
      else
      {   
        nsamples = floor( (nsteps-ntraject*ntint)/nskip );
        printf("The number of samples for this set of delays is %6.4li\n",nsamples);
      }
      
      //if flag_restart then first_iproton=iproton_restart else =1 first_isample=isample_restart else 1
      //loop over each proton <iproton>
      printf("Start processing frequency trajectories calculating response functions.\n");
      for(iproton=1;iproton<=nprotons;iproton++) for(isample=1;isample<=nsamples;isample++)
      {
        if((isample+nsamples*(iproton-1))%500==0) printf("trajectory #%ld\n",(isample+nsamples*(iproton-1)));
        
        // point dw at the appropriate row of dw_matrix
        // offset by an amount (nskip) to make this an independent,
        // uncorrlated trajectory
        dw=&dw_matrix[iproton][1+(isample-1)*nskip];
        
        /*
         * 
         * Calculate integrals in packages of ntint small intervals, use trapez method
         *
         */
        for(it=1;it<=ntraject;it++)
	    {
	      dwint[it]=(dw[(it-1)*ntint+1]+dw[it*ntint+1])/2;
	      for(j=2;j<=ntint;j++)
            dwint[it]+=dw[(it-1)*ntint+j];
	      dwint[it]*=dt/ntint;
	    }
          
        /*
         * Calculate nonlinear polarisation
         *
         */
        sdw[1]=0;
        pdw1[1][1]=0;
        pdw2[1][1]=0;
        rdw1[1][1][1]=0;
        rdw2[1][1][1]=0;
        rdw3[1][1][1]=0;
        rdw4[1][1][1]=0;
        for(it1=1;it1<=nt;it1++)
	    {
          if (it1<nt)
          {
            sdw[it1+1]=sdw[it1]+dwint[it1];
            pdw1[1][it1+1]=pdw1[1][it1]+dwint[it1];
            pdw2[1][it1+1]=pdw2[1][it1]-dwint[it1];
            rdw1[1][1][it1+1]=rdw1[1][1][it1]-dwint[it1];
            rdw2[1][1][it1+1]=rdw2[1][1][it1]+dwint[it1];
            rdw3[1][1][it1+1]=rdw3[1][1][it1]-dwint[it1];
            rdw4[1][1][it1+1]=rdw4[1][1][it1]+dwint[it1];
          }
          S_re[it1]+=cos(sdw[it1]);
          S_im[it1]+=sin(sdw[it1]);
          if (order>=3) for(it3=1;it3<=nt;it3++)
          {
            if (it3<nt)
            {
              pdw1[it3+1][it1]=pdw1[it3][it1]-dwint[it3+nt2+it1];
              pdw2[it3+1][it1]=pdw2[it3][it1]-dwint[it3+nt2+it1];
              rdw1[1][it3+1][it1]=rdw1[1][it3][it1]-dwint[it3+nt2+it1];
              rdw2[1][it3+1][it1]=rdw2[1][it3][it1]-dwint[it3+nt2+it1];
              rdw3[1][it3+1][it1]=rdw3[1][it3][it1]+dwint[it3+nt2+it1];
              rdw4[1][it3+1][it1]=rdw4[1][it3][it1]+dwint[it3+nt2+it1];
            }
            P1_re[it3][it1]+=cos(pdw1[it3][it1]);
            P1_im[it3][it1]+=sin(pdw1[it3][it1]);
            P2_re[it3][it1]+=cos(pdw2[it3][it1]);
            P2_im[it3][it1]+=sin(pdw2[it3][it1]);
            if (order>=5) for(it5=1;it5<=nt;it5++)
            {
              if (it5<nt)
              {
                rdw1[it5+1][it3][it1]=rdw1[it5][it3][it1]-dwint[it5+nt4+it3+nt2+it1];
                rdw2[it5+1][it3][it1]=rdw2[it5][it3][it1]-dwint[it5+nt4+it3+nt2+it1];
                rdw3[it5+1][it3][it1]=rdw3[it5][it3][it1]-dwint[it5+nt4+it3+nt2+it1];
                rdw4[it5+1][it3][it1]=rdw4[it5][it3][it1]-dwint[it5+nt4+it3+nt2+it1];
              }
              R1_re[it5][it3][it1]+=cos(rdw1[it5][it3][it1]);
              R1_im[it5][it3][it1]+=sin(rdw1[it5][it3][it1]);
              R2_re[it5][it3][it1]+=cos(rdw2[it5][it3][it1]);
              R2_im[it5][it3][it1]+=sin(rdw2[it5][it3][it1]);
              R3_re[it5][it3][it1]+=cos(rdw3[it5][it3][it1]);
              R3_im[it5][it3][it1]+=sin(rdw3[it5][it3][it1]);
              R4_re[it5][it3][it1]+=cos(rdw4[it5][it3][it1]);
              R4_im[it5][it3][it1]+=sin(rdw4[it5][it3][it1]);
            }
          }
        }
            if((isample+nsamples*(iproton-1))%20000==0)
            {
              /*write intermediate results*/
              normalizeResults(nt,(isample+nsamples*(iproton-1)),
                               S_re,S_im,
                               P1_re,P1_im,P2_re,P2_im,
                               R1_re,R1_im,R2_re,R2_im,R3_re,R3_im,R4_re,R4_im,
                               Sf_re,Sf_im,
                               P1f_re,P1f_im,P2f_re,P2f_im,
                               R1f_re,R1f_im,R2f_re,R2f_im,R3f_re,R3f_im,R4f_re,R4f_im);
              
              fourierTransformResults(nt,
                                      Sf_re,Sf_im,
                                      P1f_re,P1f_im,P2f_re,P2f_im,
                                      R1f_re,R1f_im,R2f_re,R2f_im,R3f_re,R3f_im,R4f_re,R4f_im);
              
              writeResults(base_name,nt,extension,
                           Sf_re,Sf_im,
                           P1f_re,P1f_im,P2f_re,P2f_im,
                           R1f_re,R1f_im,R2f_re,R2f_im,R3f_re,R3f_im,R4f_re,R4f_im); 
            } //end if (isample+nsamples*(iproton-1)%2000==0 

	    //somewhere here unset flag_restart

        } //end isample loop
          
          //final results
          normalizeResults(nt,nsamples*nprotons,
                           S_re,S_im,
                           P1_re,P1_im,P2_re,P2_im,
                           R1_re,R1_im,R2_re,R2_im,R3_re,R3_im,R4_re,R4_im,
                           Sf_re,Sf_im,
                           P1f_re,P1f_im,P2f_re,P2f_im,
                           R1f_re,R1f_im,R2f_re,R2f_im,R3f_re,R3f_im,R4f_re,R4f_im);
          
          fourierTransformResults(nt,
                                  Sf_re,Sf_im,
                                  P1f_re,P1f_im,P2f_re,P2f_im,
                                  R1f_re,R1f_im,R2f_re,R2f_im,R3f_re,R3f_im,R4f_re,R4f_im);
          
          writeResults(base_name,nt,extension,
                       Sf_re,Sf_im,
                       P1f_re,P1f_im,P2f_re,P2f_im,
                       R1f_re,R1f_im,R2f_re,R2f_im,R3f_re,R3f_im,R4f_re,R4f_im);
          
          //free dwint for the next t2,t4
          free_vector(dwint,1,ntraject);
    }//end t2, t4 loop
        
        
    printf("Free memory...\n");
    printf("...response functions R1\n");
    free_f3tensor(R1_re,1,nt,1,nt,1,nt);
    free_f3tensor(R1_im,1,nt,1,nt,1,nt);
    printf("...response functions R2\n");
    free_f3tensor(R2_re,1,nt,1,nt,1,nt);
    free_f3tensor(R2_im,1,nt,1,nt,1,nt);
    printf("...response functions R3\n");
    free_f3tensor(R3_re,1,nt,1,nt,1,nt);
    free_f3tensor(R3_im,1,nt,1,nt,1,nt);
    printf("...response functions R4\n");
    free_f3tensor(R4_re,1,nt,1,nt,1,nt);
    free_f3tensor(R4_im,1,nt,1,nt,1,nt);
    printf("...response functions P1\n");
    free_matrix(P1_re,1,nt,1,nt);
    free_matrix(P1_im,1,nt,1,nt);
    printf("...response functions P2\n");
    free_matrix(P2_re,1,nt,1,nt);
    free_matrix(P2_im,1,nt,1,nt);
    printf("...response functions S\n");
    free_vector(S_re,1,nt);
    free_vector(S_im,1,nt);
    
    printf("...response averages\n");
    free_f3tensor(R1f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R1f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R2f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R2f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R3f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R3f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R4f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R4f_im,1,nt,1,nt,1,nt);
    free_matrix(P1f_re,1,nt,1,nt);
    free_matrix(P1f_im,1,nt,1,nt);
    free_matrix(P2f_re,1,nt,1,nt);
    free_matrix(P2f_im,1,nt,1,nt);
    free_vector(Sf_re,1,nt);
    free_vector(Sf_im,1,nt);
    
    printf("...response intermediates\n");
    free_vector(sdw,1,nt);
    free_matrix(    pdw1,1,nt,1,nt);
    free_matrix(    pdw2,1,nt,1,nt);
    free_f3tensor(    rdw1,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw2,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw3,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw4,1,nt,1,nt,1,nt);
    
    printf("...dw_matrix\n");
    free_matrix(dw_matrix,1,nprotons,1,nsteps);
    
    free_matrix(t2_t4_pairs,1,n_t2_t4_pairs,1,n_t2_t4_pairs);
    //free_lvector(t2_array,1,n_t2_array);
    //free_lvector(tdiag_array,1,n_t2_array);
    
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
  
  /* do the fft */
  four1(a,n*nfill,+1);

  /* now take n/2 points from each end and shift them to the middle of the spectrum
     like fftshift of matlab*/
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

