#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>
#include <fnmatch.h>
#include "globalArgs.h"
#include "mymath.h"

// 14.9.05: Filled Fouriertransform with linear interpolation, reduces dramatically the offset problem
/*  
*  Original program fifthorder.c by P. Hamm calculated response functions from 
*  Langevin dynamics trajectories. 
*  8.1.2007 Modified by S. Garrett-Roe to read an MD trajectory (one file of coordinates and 
                                                                 *  one file of forces) and calculate response functions from that.
*
*/



#define wavenumbersToInvPs 2.99792458e-2

void display_usage( int argc, char* argv[] )
{
  puts("freqTrajToR5 - calculate a response function from an input frequency trajectory\n"
       "\n"
       "Usage:\n");
  printf("USAGE:\n   %s [-g nt] [-i ntint] [-n nsteps] [-q nprotons] [-r iproton] [-s nskip]  -t time_file -o output_file_prefix\n\n",argv[0]);
  puts(" other options:\n"
       "-e --nsteps_in_file number of time steps in freq trajectory\n"
       "-f --nprotons_in_file number of protons in freq trajectory\n"
       "-g --nt number of time points per time axis in response function\n"
       "--dt the timestep from the gromacs file\n"
       "-i --ntint number of steps preintegrated\n"
       "-n --nsteps number of steps of the trajectory to use\n"
       "-q --nprotons number of protons in the trajectory to use\n"
       "-r --proton_offset proton offset in the trajectory to use\n"
       "-s --nskip number of steps to skip between mini-trajectories\n"
       "-t --time_file name of time file\n");
  
printf("PREVIOUS USAGE:\n   %s [-g nt] [-i ntint] [-n nsteps]"
	     "[-p nprotons] [-s nskip] [-t dt] [-f nprotons_in_file]"
	     "[-o proton_offset] [-v verbosity_level ] [-z flagFFT]"
	     "freq_file time_file output_file_prefix\n\n",argv[0]);
       
  exit(EXIT_FAILURE);
}

static const char *optString = "e:f:g:i:n:o:p:r:s:t:vh?";

static const struct option longOpts[] = {
  // all the long options go here with the letter that they map do (or
  //  0 if they have no short version)
  { "nsteps_in_file", required_argument,NULL,'e'},
  { "nprotons_in_file", required_argument,NULL,'f'},
  { "nt", required_argument,NULL,'g'},
  { "nprotons", required_argument,NULL,'q'},
  { "proton_offset", required_argument,NULL,'r'},
  { "ntint", required_argument,NULL,'i'},
  { "nsteps", required_argument,NULL,'n'},
  { "nskip", required_argument,NULL,'s'},
  { "time_file",required_argument,NULL,'t'},
  { "dt",required_argument,NULL,0},
  { "output",required_argument,NULL,'o'},
  { "verbose",no_argument,NULL,'v'},
  { "help",no_argument,NULL, 'h'},
  { NULL, no_argument, NULL,0}
};

void read_time_file(const char *time_file_name,float ***t2_t4_pairs,int *n_t2_t4_pairs){
  FILE *fid;
  char *string;
  int i;

  if (DEBUG_LEVEL >= 1) printf("reading time file...\n");

  if (fnmatch(time_file_name,"",FNM_CASEFOLD)==0){
    fprintf(stderr,"no time file provided\n");
    exit(EXIT_FAILURE);
  }

    /*
     * Setup the time axes for the joint probability density
     */
    //set up the time axis that will be used
    //note: because of the way the response functions are calculated,
    // it is better to have t2 AND t2/2 be multiples of ntint
    if (asprintf(&string,"wc -l %s",time_file_name)==0)
      {
	fprintf(stderr,"unable to write to string\n");
	exit(EXIT_FAILURE);	
      }

    fid=popen(string,"r");
    if(fid==NULL) nrerror("Failed opening pipe wc -l time_file_name");
    fscanf(fid,"%d",n_t2_t4_pairs);
    pclose(fid);
    free(string);
    if (DEBUG_LEVEL>=1) printf("n_t2_t4_pairs = %i\n",*n_t2_t4_pairs);
    
    *t2_t4_pairs = matrix(1,*n_t2_t4_pairs,1,2);

    fid=fopen(time_file_name,"rt");
    if (fid==NULL) nrerror("Failed to open time_file for reading.");
    for (i=1;i<=*n_t2_t4_pairs;i++)
      {
	if(fscanf(fid,"%f %f",&((*t2_t4_pairs)[i][1]),&((*t2_t4_pairs)[i][2]))!=2)
	  nrerror("Failed reading time_file before end of file.");
	if (DEBUG_LEVEL>=1) printf("t2_t4_pairs[%i][1,2] = (%f,%f)\n",*n_t2_t4_pairs,(*t2_t4_pairs)[i][1],(*t2_t4_pairs)[i][2]);

      }
    fclose(fid);
}

void read_nsteps_and_nprotons_from_freq( char *parameter_file_name ){
  char *freq_file_name;
  FILE *fid;
  char *string;

    /* number of steps to use from the trajectory */
  int nsteps_in_file = globalArgs.nsteps_in_file; 
  /* number of steps to use from the trajectory */
  int nprotons_in_file = globalArgs.nprotons_in_file; 

  int flag_compressedinput = 0; //not the global one!

  freq_file_name = globalArgs.w_file_names[0];

  /* if not specified on the comand line... */
  if (nsteps_in_file==-1){
    //use the shell command "wc -l" to scan the length of the coord file to get the number of steps
    // and read one line with "head" and then "wc -w" to get the number of molecules
    printf("Determine number of steps and number of molecules from frequency file.\n");

    if (fnmatch(freq_file_name,"*.gz",FNM_CASEFOLD) == 0)
      flag_compressedinput = 1;


    if (flag_compressedinput == 0){
      //for straight ascii
      if (asprintf(&string,"wc -l %s",freq_file_name) < 0) 
	nrerror("failed to write string");
    } else {
      //for gzipped input
      if (asprintf(&string,"gzip -cd %s | wc -l",freq_file_name) < 0) 
	nrerror("failed to write string");
    }
    fid = popen(string,"r");
    if(fid==NULL)
      {
	fprintf(stderr,"Opening %s failed.",string);
	exit(EXIT_FAILURE);
      }
    fscanf(fid,"%d",&nsteps_in_file);
    pclose(fid);
    free(string);
    printf("Found %d steps in file %s\n", nsteps_in_file, freq_file_name);
    globalArgs.nsteps_in_file = nsteps_in_file;
  }
  if (nprotons_in_file==-1){
    if (flag_compressedinput==0){
      if (asprintf(&string,"head -n 1 %s | wc -w",freq_file_name) == 0) 
	nrerror("failed to write string");
    } else {
      if (asprintf(&string,"gzip -cd %s | head -n 1 | wc -w",freq_file_name) == 0) 
	nrerror("failed to write string");
    }
    fid = popen(string,"r");
    if(fid==NULL)
      {
	fprintf(stderr,"Opening %s failed.",string);
	exit(EXIT_FAILURE);
      }
    fscanf(fid,"%d",&nprotons_in_file);
    pclose(fid);
    free(string);
    printf("Found %d protons in file %s\n", nprotons_in_file, freq_file_name);
    globalArgs.nprotons_in_file = nprotons_in_file;
  }

  //repeat for x file to make sure they are the same???

  // add these values to the parameter file
  //  globalArgs.nmols_in_file = nmols_in_file;
  //  gaWriteInt(parameter_file_name,"nmols_in_file",nmols_in_file); //only output one freq per molecule
  globalArgs.nprotons_in_file = nprotons_in_file;
  gaWriteInt(parameter_file_name,"nprotons_in_file",nprotons_in_file); //only output one freq per molecule
  globalArgs.nsteps_in_file = nsteps_in_file;
  gaWriteInt(parameter_file_name,"nsteps_in_file",nsteps_in_file);
}

void read_x_files(int nprotons, int proton_offset, int nsteps, int step_offset,float ***x_matrices,float *mean_x){
  FILE *fid;
  char *string;
  char **mu_file_names = globalArgs.mu_file_names;
  char *freq_file_name;

  /* number of excited states */
  const int n_levels = globalArgs.n_levels;
  /* number of steps to use from the trajectory */
  int nsteps_in_file = globalArgs.nsteps_in_file; 
  /* number of steps to use from the trajectory */
  int nprotons_in_file = globalArgs.nprotons_in_file; 

  int flag_compressedinput = 0; //not the global one!

  int i,j,i_level;
  float mean,sigma,skewness; //characterize the frequency traj
  double sum;
  float **dw_matrix;
  float **dw_read;

  if (DEBUG_LEVEL >= 0) printf("\nreading dipole file...\n");

  if (fnmatch("*.gz",mu_file_names[0],FNM_CASEFOLD) == 0)
    flag_compressedinput = 1;

  // matrix to read in entire file
  dw_read = matrix(1,nprotons_in_file,1,nsteps_in_file);
  
  //for each level create a matrix of frequencies
  for (i_level=0;i_level<n_levels;i_level++){
    // smaller matrices for selected parts of trajectory
    x_matrices[i_level] = matrix(1,nprotons,1,nsteps);
    dw_matrix = x_matrices[i_level]; //TEST TEST

    freq_file_name = mu_file_names[i_level];

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

    if (flag_compressedinput == 0){
      if (asprintf(&string,"%s",freq_file_name) < 0 ) nrerror("failed to write string");
      fid = fopen(string,"rt");
    } else {
      if (asprintf(&string,"gzip -cd %s",freq_file_name) < 0 ) nrerror("failed to write string");
      fid = popen(string,"r");
    }
    if(fid==NULL) nrerror(string);
    free(string);

    // read in entire trajectory
    for(j=1;j<=nsteps;j++)
      for(i=1;i<=nprotons_in_file;i++)
	{
	  fscanf(fid,"%f",&dw_read[i][j]);
	}

    //close file or pipe
    if (flag_compressedinput == 0){
      fclose(fid);    
    } else {
      pclose(fid);    
    }      
        
    //grab the desired part of the trajectory
    for(j=1;j<=nsteps;j++)
      for(i=1;i<=nprotons;i++)
	{
	  dw_matrix[i][j]=dw_read[proton_offset+i][step_offset+j];
	}

    // calculate the mean and subtract it
    sum = 0;
    for(i=1;i<=nprotons;i++)
    {
      for(j=1;j<=nsteps;j++)
        sum = sum + dw_matrix[i][j];
    }
    mean = sum/(nsteps*nprotons);
    mean_x[i_level] = mean;

/*     for(i=1;i<=nprotons;i++) */
/*       for(j=1;j<=nsteps;j++) */
/*         dw_matrix[i][j] = dw_matrix[i][j] - mean; */
    printf("The mean dipole of the trajectory, counting all protons is %f \n",mean);
    
    //calculate the standard deviation for diagonstic purposes
    sum = 0;
    for(i=1;i<=nprotons;i++)
    {
      for(j=1;j<=nsteps;j++)
        sum = sum + dw_matrix[i][j]*dw_matrix[i][j];
    }
    sigma = sqrt(sum/(nsteps*nprotons));
    printf("The standard deviation of the trajectory, counting all protons is %f \n",sigma);
    
    
    //calculate the skewness for diagonstic purposes
    sum = 0;
    for(i=1;i<=nprotons;i++)
      for(j=1;j<=nsteps;j++)
        sum = sum + dw_matrix[i][j]*dw_matrix[i][j]*dw_matrix[i][j];
    skewness = sum/(nsteps*nprotons*pow(sigma,3));
    printf("The skewness of the trajectory, counting all protons is %f.\n",skewness);

  } //end for i_levels
  
  //free up memory
  printf("Free memory ... dw_read\n");
  free_matrix(dw_read,1,nprotons_in_file,1,nsteps);

}

void read_w_files(int nprotons, int proton_offset, int nsteps, int step_offset,float ***w_matrices,float *mean_w){

  FILE *fid;
  char *string;
  char **w_file_names = globalArgs.w_file_names;
  char *freq_file_name;
  
  /* number of excited states */
  const int n_levels = globalArgs.n_levels;
  /* number of steps to use from the trajectory */
  int nsteps_in_file = globalArgs.nsteps_in_file; 
  /* number of steps to use from the trajectory */
  int nprotons_in_file = globalArgs.nprotons_in_file; 

  int flag_compressedinput = 0; //not the global one!

  //float *dw_initial; //the initial frequency of each proton
  int i,j,i_level;
  float mean,sigma,skewness; //characterize the frequency traj
  double sum;
  float **dw_matrix;
  float **dw_read;

  if (DEBUG_LEVEL >= 0) printf("\nreading freq file...\n");
  if (DEBUG_LEVEL >= 2) 
    for (i=0;i<n_levels;i++)
      printf("g.fname =  %s fname = %s\n",globalArgs.w_file_names[i],w_file_names[i]);

  if (fnmatch("*.gz",w_file_names[0],FNM_CASEFOLD) == 0)
    flag_compressedinput = 1;
  
  // matrix to read in entire file
  dw_read = matrix(1,nprotons_in_file,1,nsteps_in_file);

  //for each level create a matrix of frequencies
  for (i_level=0;i_level<n_levels;i_level++){
    
    // smaller matrices for selected parts of trajectory
    w_matrices[i_level] = matrix(1,nprotons,1,nsteps);
    dw_matrix = w_matrices[i_level]; //TEST TEST

    freq_file_name = w_file_names[i_level];
    if (DEBUG_LEVEL >= 2) printf("... file %s\n",freq_file_name);

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
    if (flag_compressedinput == 0){
      if (asprintf(&string,"%s",freq_file_name) < 0 ) nrerror("failed to write string");
      fid = fopen(string,"rt");
    } else {
      if (asprintf(&string,"gzip -cd %s",freq_file_name) < 0 ) nrerror("failed to write string");
      fid = popen(string,"r");
    }
    if (DEBUG_LEVEL >= 2) printf("using %s = %p\n",string,fid);
    if(fid==NULL) nrerror(string);
    free(string);

    // read in entire trajectory
    for(j=1;j<=nsteps_in_file;j++)
      for(i=1;i<=nprotons_in_file;i++)
	{
	  fscanf(fid,"%f",&dw_read[i][j]);
	}

    //close file or pipe
    if (flag_compressedinput == 0){
      fclose(fid);    
    } else {
      pclose(fid);    
    }      
    
    //grab the desired part of the trajectory
    for(j=1;j<=nsteps;j++)
      for(i=1;i<=nprotons;i++)
	{
	  dw_matrix[i][j]=dw_read[proton_offset+i][step_offset+j];
	}

    // calculate the mean and subtract it
    sum = 0;
    for(i=1;i<=nprotons;i++)
    {
      for(j=1;j<=nsteps;j++)
        sum = sum + dw_matrix[i][j];
    }
    mean = sum/(nsteps*nprotons);
    mean_w[i_level] = mean;

    for(i=1;i<=nprotons;i++)
      for(j=1;j<=nsteps;j++)
        dw_matrix[i][j] = dw_matrix[i][j] - mean;
    printf("The mean frequency of the trajectory, counting all protons is %f cm-1\n",mean);
    
    //calculate the standard deviation for diagonstic purposes
    sum = 0;
    for(i=1;i<=nprotons;i++)
    {
      for(j=1;j<=nsteps;j++)
        sum = sum + dw_matrix[i][j]*dw_matrix[i][j];
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

  } //end for i_levels
  
  //free up memory
  printf("Free mememory ... dw_read\n");
  free_matrix(dw_read,1,nprotons_in_file,1,nsteps);

}


void writeResultsTime(const char base_name[],const int nt,const char extension[],
		      const int n_levels,
		      float Sf_re[],float Sf_im[],
		      float **P1f_re,float **P1f_im,float **P2f_re,float **P2f_im,
		      float **P3f_re,float **P3f_im,float **P4f_re,float **P4f_im,
		      float **P1tot_re,float **P1tot_im,float **P2tot_re,float **P2tot_im,
		      float ***R1f_re,float ***R1f_im,float ***R2f_re,float ***R2f_im,
		      float ***R3f_re,float ***R3f_im,float ***R4f_re,float ***R4f_im,
		      float ***R5f_re, float ***R5f_im, float ***R6f_re, float ***R6f_im,
		      float ***R7f_re, float ***R7f_im, float ***R8f_re, float ***R8f_im,
		      float ***R9f_re, float ***R9f_im, float ***R10f_re,float ***R10f_im,
		      float ***R11f_re,float ***R11f_im,float ***R12f_re,float ***R12f_im,
		      float ***R13f_re,float ***R13f_im,float ***R14f_re,float ***R14f_im,
		      float ***R15f_re,float ***R15f_im,float ***R16f_re,float ***R16f_im,
		      float ***R17f_re,float ***R17f_im,float ***R18f_re,float ***R18f_im,
		      float ***R19f_re,float ***R19f_im,float ***R20f_re,float ***R20f_im,
                      float ***R1tot_re,float ***R1tot_im,float ***R2tot_re,float ***R2tot_im,
                      float ***R3tot_re,float ***R3tot_im,float ***R4tot_re,float ***R4tot_im)
/* This version is for outputting the data in the time domain with no
 * flipping or FT or anything like that, just the data as it is
 * calculated with the highest time argument as the fastest changing
 * element (t5 then t3 then t1) */
{
  char *name;
  int it1,it3,it5;
  FILE *out1D,*out2D,*out3D;
  
  printf("Write Results\n");
  if (asprintf(&name,"%s_t_spec1D%s",base_name,extension) < 0) nrerror("failed to write string");

  /* 1D */
  out1D=fopen(name,"wt");
  for(it1=1;it1<=nt;it1++)
  {
    fprintf(out1D,"%10.5g %10.5g",Sf_re[it1],Sf_im[it1]);
    fprintf(out1D,"\n");
  }
  fclose(out1D);
  free(name);

  /* 2D gsb + se (which are always calculated) */
  if (asprintf(&name,"%s_t_spec2D_peak1%s",base_name,extension) < 0) nrerror("failed to write string");
  out2D=fopen(name,"wt");
  for(it1=1;it1<=nt;it1++)
    for(it3=1;it3<=nt;it3++)
    {
      fprintf(out2D,"%10.5g %10.5g %10.5g %10.5g",
	      P1f_re[it3][it1],P1f_im[it3][it1],
	      P2f_re[it3][it1],P2f_im[it3][it1]);
      fprintf(out2D,"\n");
    }
  fclose(out2D);
  free(name);

  /* if we calculated higher states print them out, too */
  if (n_levels>=2){
    /* 2D esa */
    if (asprintf(&name,"%s_t_spec2D_peak2%s",base_name,extension) < 0) nrerror("failed to write string");
    out2D=fopen(name,"wt");
    for(it1=1;it1<=nt;it1++)
      for(it3=1;it3<=nt;it3++)
	{
	  fprintf(out2D,"%10.5g %10.5g %10.5g %10.5g",
		  P3f_re[it3][it1],P3f_im[it3][it1],
		  P4f_re[it3][it1],P4f_im[it3][it1]);
	  fprintf(out2D,"\n");
	}
    fclose(out2D);
    free(name);

    /* 2D total spectrum */
    if (asprintf(&name,"%s_t_spec2D%s",base_name,extension) < 0) nrerror("failed to write string");
    out2D=fopen(name,"wt");
    for(it1=1;it1<=nt;it1++)
      for(it3=1;it3<=nt;it3++)
	{
	  fprintf(out2D,"%10.5g %10.5g %10.5g %10.5g",
		  P1tot_re[it3][it1],P1tot_im[it3][it1],
		  P2tot_re[it3][it1],P2tot_im[it3][it1]);
	  fprintf(out2D,"\n");
	}
    fclose(out2D);
    free(name);

  if (asprintf(&name,"%s_t_spec3D_peak1%s",base_name,extension) < 0) nrerror("failed to write string");
  out3D=fopen(name,"wt");
  for(it1=1;it1<=nt;it1++)
    for(it3=1;it3<=nt;it3++)
      for(it5=1;it5<=nt;it5++)
      {
        fprintf(out3D,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g",
                R1f_re[it5][it3][it1],
                R1f_im[it5][it3][it1],
                R2f_re[it5][it3][it1],
                R2f_im[it5][it3][it1],
                R3f_re[it5][it3][it1],
                R3f_im[it5][it3][it1],
                R4f_re[it5][it3][it1],
                R4f_im[it5][it3][it1]);
        fprintf(out3D,"\n");
      }
  fclose(out3D);
  free(name);

  if (asprintf(&name,"%s_t_spec3D_peak2%s",base_name,extension) < 0) nrerror("failed to write string");
  out3D=fopen(name,"wt");
  for(it1=1;it1<=nt;it1++)
    for(it3=1;it3<=nt;it3++)
      for(it5=1;it5<=nt;it5++)
      {
        fprintf(out3D,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g",
                R5f_re[it5][it3][it1],
                R5f_im[it5][it3][it1],
                R6f_re[it5][it3][it1],
                R6f_im[it5][it3][it1],
                R7f_re[it5][it3][it1],
                R7f_im[it5][it3][it1],
                R8f_re[it5][it3][it1],
                R8f_im[it5][it3][it1]);
        fprintf(out3D,"\n");
      }
  fclose(out3D);
  free(name);

  if (n_levels>=3){
    if (asprintf(&name,"%s_t_spec3D_peak3%s",base_name,extension) < 0) nrerror("failed to write string");
    out3D=fopen(name,"wt");
    for(it1=1;it1<=nt;it1++)
      for(it3=1;it3<=nt;it3++)
	for(it5=1;it5<=nt;it5++)
	  {
	    fprintf(out3D,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g",
		    R9f_re[it5][it3][it1],
		    R9f_im[it5][it3][it1],
		    R10f_re[it5][it3][it1],
		    R10f_im[it5][it3][it1],
		    R11f_re[it5][it3][it1],
		    R11f_im[it5][it3][it1],
		    R12f_re[it5][it3][it1],
		    R12f_im[it5][it3][it1]);
	    fprintf(out3D,"\n");
	  }
    fclose(out3D);
    free(name);
  } /* end n_levels>=3 */

  if (asprintf(&name,"%s_t_spec3D_peak4%s",base_name,extension) < 0) nrerror("failed to write string");
  out3D=fopen(name,"wt");
  for(it1=1;it1<=nt;it1++)
    for(it3=1;it3<=nt;it3++)
      for(it5=1;it5<=nt;it5++)
      {
        fprintf(out3D,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g",
                R13f_re[it5][it3][it1],
                R13f_im[it5][it3][it1],
                R14f_re[it5][it3][it1],
                R14f_im[it5][it3][it1],
                R15f_re[it5][it3][it1],
                R15f_im[it5][it3][it1],
                R16f_re[it5][it3][it1],
                R16f_im[it5][it3][it1]);
        fprintf(out3D,"\n");
      }
  fclose(out3D);
  free(name);
  if (asprintf(&name,"%s_t_spec3D_peak5%s",base_name,extension) < 0) nrerror("failed to write string");
  out3D=fopen(name,"wt");
  for(it1=1;it1<=nt;it1++)
    for(it3=1;it3<=nt;it3++)
      for(it5=1;it5<=nt;it5++)
      {
        fprintf(out3D,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g",
                R17f_re[it5][it3][it1],
                R17f_im[it5][it3][it1],
                R18f_re[it5][it3][it1],
                R18f_im[it5][it3][it1],
                R19f_re[it5][it3][it1],
                R19f_im[it5][it3][it1],
                R20f_re[it5][it3][it1],
                R20f_im[it5][it3][it1]);
        fprintf(out3D,"\n");
      }
  fclose(out3D);
  free(name);
  } /* end if n_levels>=2 */

  if (asprintf(&name,"%s_t_spec3D%s",base_name,extension) < 0) nrerror("failed to write string");
  out3D=fopen(name,"wt");
  for(it1=1;it1<=nt;it1++)
    for(it3=1;it3<=nt;it3++)
      for(it5=1;it5<=nt;it5++)
      {
        fprintf(out3D,"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g",
                R1tot_re[it5][it3][it1],
                R1tot_im[it5][it3][it1],
                R2tot_re[it5][it3][it1],
                R2tot_im[it5][it3][it1],
                R3tot_re[it5][it3][it1],
                R3tot_im[it5][it3][it1],
                R4tot_re[it5][it3][it1],
                R4tot_im[it5][it3][it1]);
        fprintf(out3D,"\n");
      }
  fclose(out3D);
  free(name);

}

void normalizeResults(int nt,float dt,unsigned long isample,
		      int n_levels, int flag_noncondon,float *mean_w,float *mean_mu,
                      float S_re[],float S_im[],
                      float **P1_re,float **P1_im,float **P2_re,float **P2_im,
                      float **P3_re,float **P3_im,float **P4_re,float **P4_im,
                      float ***R1_re, float ***R1_im, float ***R2_re, float ***R2_im,
                      float ***R3_re, float ***R3_im, float ***R4_re, float ***R4_im,
		      float ***R5_re, float ***R5_im, float ***R6_re, float ***R6_im,
		      float ***R7_re, float ***R7_im, float ***R8_re, float ***R8_im,
		      float ***R9_re, float ***R9_im, float ***R10_re,float ***R10_im,
		      float ***R11_re,float ***R11_im,float ***R12_re,float ***R12_im,
		      float ***R13_re,float ***R13_im,float ***R14_re,float ***R14_im,
		      float ***R15_re,float ***R15_im,float ***R16_re,float ***R16_im,
		      float ***R17_re,float ***R17_im,float ***R18_re,float ***R18_im,
		      float ***R19_re,float ***R19_im,float ***R20_re,float ***R20_im,
                      float Sf_re[],float Sf_im[],
                      float **P1f_re,float **P1f_im,float **P2f_re,float **P2f_im,
                      float **P3f_re,float **P3f_im,float **P4f_re,float **P4f_im,
                      float **P1tot_re,float **P1tot_im,float **P2tot_re,float **P2tot_im,
                      float ***R1f_re, float ***R1f_im, float ***R2f_re, float ***R2f_im,
                      float ***R3f_re, float ***R3f_im, float ***R4f_re, float ***R4f_im,
		      float ***R5f_re, float ***R5f_im, float ***R6f_re, float ***R6f_im,
		      float ***R7f_re, float ***R7f_im, float ***R8f_re, float ***R8f_im,
		      float ***R9f_re, float ***R9f_im, float ***R10f_re,float ***R10f_im,
		      float ***R11f_re,float ***R11f_im,float ***R12f_re,float ***R12f_im,
		      float ***R13f_re,float ***R13f_im,float ***R14f_re,float ***R14f_im,
		      float ***R15f_re,float ***R15f_im,float ***R16f_re,float ***R16f_im,
		      float ***R17f_re,float ***R17f_im,float ***R18f_re,float ***R18f_im,
		      float ***R19f_re,float ***R19f_im,float ***R20f_re,float ***R20f_im,
                      float ***R1tot_re,float ***R1tot_im,float ***R2tot_re,float ***R2tot_im,
                      float ***R3tot_re,float ***R3tot_im,float ***R4tot_re,float ***R4tot_im)
/* normalize for number of samples, frequency shifts from rotating
   frame, amplitudes due to dipoles if condon, and number of
   population pathways */
{
  int it1,it3,it5;
  float mu_01_2,mu_12_2,mu_23_2;
  float shift_w,shift_w_1_1,shift_w_1_2,shift_w_2_1,shift_w_2_2,shift_w_2_3;
  float a,b,phi,phi1,phi2,phi3,phi4,phi5;
  float mw;

  mu_01_2 = 1;
  mu_12_2 = 1;
  mu_23_2 = 1;
  if (flag_noncondon==0){
    if (n_levels>=1) mu_01_2 = mean_mu[0]*mean_mu[0];
    if (n_levels>=2) mu_12_2 = mean_mu[1]*mean_mu[1];
    if (n_levels>=3) mu_23_2 = mean_mu[2]*mean_mu[2];
  }
  

  /* defaults */
  shift_w = 0;
  shift_w_1_1 = 0;
  shift_w_1_2 = 0;
  shift_w_2_1 = 0;
  shift_w_2_2 = 0;
  shift_w_2_3 = 0;
  /* calculate shift frequency based on number of levels we've calculated */
  if (n_levels==1){
    /* for 2D spectra */
    shift_w = 0;
    /* for 3D spectra */
    shift_w_1_1 = 0;
    shift_w_1_2 = 0;
    shift_w_2_1 = 0;
    shift_w_2_2 = 0;
    shift_w_2_3 = 0;
  }
  if (n_levels==2){
    /* for 2D spectra */
    shift_w =(mean_w[0] - mean_w[1])/2; // 0;//  M_PI/(2*dt);// 
    /* for 3D spectra */
    shift_w_1_1 = shift_w;
    shift_w_1_2 = -shift_w;
    shift_w_2_1 = shift_w;
    shift_w_2_2 = -shift_w;
    shift_w_2_3 = 0;
  }
  if (n_levels==3){
    /* for 2D spectra */
    shift_w =(mean_w[0] - mean_w[1])/2; // 0;//  M_PI/(2*dt);// 

    /* for 3D spectra */
    mw = (mean_w[0] + mean_w[1] + mean_w[2])/3;
    shift_w_1_1 = shift_w; /*this needs to be worked out*/
    shift_w_1_2 = -shift_w; /*this needs to be worked out*/
    shift_w_2_1 = -1*(mean_w[0] - mw); /*this needs to be worked out*/
    shift_w_2_2 = -1*(mean_w[1] - mw); /*this needs to be worked out*/
    shift_w_2_3 = -1*(mean_w[2] - mw); /*this needs to be worked out*/
    /* take the average of the frequencies and shift each by the difference of the frequency from the mean */
    /* shift_w_1_1  = +Delta/2; shift_w_1_2 = -Delta/2; shift_w_2_1 = +Delta; shift_w_2_2 = 0; shift_w_2_3 = -Delta;*/
  }

  printf("shift_w = %f dt = %f\n",shift_w,dt);

  for(it1=1;it1<=nt;it1++)
    {
      Sf_re[it1] = mu_01_2 * S_re[it1] / isample;
      Sf_im[it1] = mu_01_2 * S_im[it1] / isample;
      for(it3=1;it3<=nt;it3++)
	{
	  /* the factor of 2 is for the 2 diagrams which contribute,
	     ground state bleach and excited state emission. The a and
	     b terms are are do the rotation in the complex plane to
	     account for the frequency shifts due to anharmonicity. */
	  phi = -shift_w*dt*(it3-1);
	  a = 2 * mu_01_2 * mu_01_2 * P1_re[it3][it1] / isample;
	  b = 2 * mu_01_2 * mu_01_2 * P1_im[it3][it1] / isample;
	  P1f_re[it3][it1] = a * cos(phi) - b * sin(phi);
	  P1f_im[it3][it1] = a * sin(phi) + b * cos(phi);

	  a = 2 * mu_01_2 * mu_01_2 * P2_re[it3][it1] / isample;
	  b = 2 * mu_01_2 * mu_01_2 * P2_im[it3][it1] / isample;
	  P2f_re[it3][it1] = a * cos(phi) - b * sin(phi);
	  P2f_im[it3][it1] = a * sin(phi) + b * cos(phi);

	  if (n_levels>=2){
	    /* the -1 factor is because of an odd number of arrows on
	       the left of the diagrams (ie it is excited state
	       absorption */
	    a = -1 * mu_12_2 * mu_01_2 * P3_re[it3][it1] / isample;
	    b = -1 * mu_12_2 * mu_01_2 * P3_im[it3][it1] / isample;
	    P3f_re[it3][it1] = a * cos(-phi) - b * sin(-phi);
	    P3f_im[it3][it1] = a * sin(-phi) + b * cos(-phi);
	    a = -1 * mu_12_2 * mu_01_2 * P4_re[it3][it1] / isample;
	    b = -1 * mu_12_2 * mu_01_2 * P4_im[it3][it1] / isample;
	    P4f_re[it3][it1] = a * cos(-phi) - b * sin(-phi);
	    P4f_im[it3][it1] = a * sin(-phi) + b * cos(-phi);

	    P1tot_re[it3][it1] = P1f_re[it3][it1] + P3f_re[it3][it1];
	    P1tot_im[it3][it1] = P1f_im[it3][it1] + P3f_im[it3][it1];
	    P2tot_re[it3][it1] = P2f_re[it3][it1] + P4f_re[it3][it1];
	    P2tot_im[it3][it1] = P2f_im[it3][it1] + P4f_im[it3][it1];
	  }
	  for(it5=1;it5<=nt;it5++)
	    {
	      phi1 = -shift_w_1_1*dt*(it3-1) - shift_w_2_1*dt*(it5-1);
	      a = 4 * mu_01_2 * mu_01_2 * mu_01_2 * R1_re[it5][it3][it1] / isample;
	      b = 4 * mu_01_2 * mu_01_2 * mu_01_2 * R1_im[it5][it3][it1] / isample;
	      R1f_re[it5][it3][it1]  = a * cos(phi1) - b * sin(phi1);
	      R1f_im[it5][it3][it1]  = a * sin(phi1) + b * cos(phi1);
	      a = 4 * mu_01_2 * mu_01_2 * mu_01_2 * R2_re[it5][it3][it1] / isample;
	      b = 4 * mu_01_2 * mu_01_2 * mu_01_2 * R2_im[it5][it3][it1] / isample;
	      R2f_re[it5][it3][it1]  =  a * cos(phi1) - b * sin(phi1);
	      R2f_im[it5][it3][it1]  =  a * sin(phi1) + b * cos(phi1);

	      phi1 = shift_w_1_1*dt*(it3-1) - shift_w_2_1*dt*(it5-1);
	      a = 4 * mu_01_2 * mu_01_2 * mu_01_2 * R3_re[it5][it3][it1] / isample;
	      b = 4 * mu_01_2 * mu_01_2 * mu_01_2 * R3_im[it5][it3][it1] / isample;
	      R3f_re[it5][it3][it1]  = a * cos(phi1) - b * sin(phi1);
	      R3f_im[it5][it3][it1]  = a * sin(phi1) + b * cos(phi1);
	      a = 4 * mu_01_2 * mu_01_2 * mu_01_2 * R4_re[it5][it3][it1] / isample;
	      b = 4 * mu_01_2 * mu_01_2 * mu_01_2 * R4_im[it5][it3][it1] / isample;
	      R4f_re[it5][it3][it1]  =  a * cos(phi1) - b * sin(phi1);
	      R4f_im[it5][it3][it1]  =  a * sin(phi1) + b * cos(phi1);

	      if (n_levels>=2){
		phi2 = -shift_w_1_1*dt*(it3-1) - shift_w_2_2*dt*(it5-1);
		a = -2 * mu_12_2 * mu_01_2 * mu_01_2 * R5_re[it5][it3][it1] / isample;
		b = -2 * mu_12_2 * mu_01_2 * mu_01_2 * R5_im[it5][it3][it1] / isample;
		R5f_re[it5][it3][it1]  = a * cos(phi2) - b * sin(phi2);
		R5f_im[it5][it3][it1]  = a * sin(phi2) + b * cos(phi2);
		a = -2 * mu_12_2 * mu_01_2 * mu_01_2 * R6_re[it5][it3][it1] / isample;
		b = -2 * mu_12_2 * mu_01_2 * mu_01_2 * R6_im[it5][it3][it1] / isample;
		R6f_re[it5][it3][it1]  =  a * cos(phi2) - b * sin(phi2);
		R6f_im[it5][it3][it1]  =  a * sin(phi2) + b * cos(phi2);

		phi2 = shift_w_1_1*dt*(it3-1) - shift_w_2_2*dt*(it5-1);
		a = -2 * mu_12_2 * mu_01_2 * mu_01_2 * R7_re[it5][it3][it1] / isample;
		b = -2 * mu_12_2 * mu_01_2 * mu_01_2 * R7_im[it5][it3][it1] / isample;
		R7f_re[it5][it3][it1]  = a * cos(phi2) - b * sin(phi2);
		R7f_im[it5][it3][it1]  = a * sin(phi2) + b * cos(phi2);
		a = -2 * mu_12_2 * mu_01_2 * mu_01_2 * R8_re[it5][it3][it1] / isample;
		b = -2 * mu_12_2 * mu_01_2 * mu_01_2 * R8_im[it5][it3][it1] / isample;
		R8f_re[it5][it3][it1]  =  a * cos(phi2) - b * sin(phi2);
		R8f_im[it5][it3][it1]  =  a * sin(phi2) + b * cos(phi2);

		if (n_levels>=3){
		  phi3 = -shift_w_1_2*dt*(it3-1) - shift_w_2_3*dt*(it5-1);
		  a = 1 * mu_23_2 * mu_12_2 * mu_01_2 * R9_re[it5][it3][it1] / isample;
		  b = 1 * mu_23_2 * mu_12_2 * mu_01_2 * R9_im[it5][it3][it1] / isample;
		  R9f_re[it5][it3][it1]  = a * cos(phi3) - b * sin(phi3);
		  R9f_im[it5][it3][it1]  = a * sin(phi3) + b * cos(phi3);
		  a = 1 * mu_23_2 * mu_12_2 * mu_01_2 * R10_re[it5][it3][it1] / isample;
		  b = 1 * mu_23_2 * mu_12_2 * mu_01_2 * R10_im[it5][it3][it1] / isample;
		  R10f_re[it5][it3][it1]  =  a * cos(phi3) - b * sin(phi3);
		  R10f_im[it5][it3][it1]  =  a * sin(phi3) + b * cos(phi3);

		  phi3 = shift_w_1_2*dt*(it3-1) - shift_w_2_3*dt*(it5-1);
		  a = 1 * mu_23_2 * mu_12_2 * mu_01_2 * R11_re[it5][it3][it1] / isample;
		  b = 1 * mu_23_2 * mu_12_2 * mu_01_2 * R11_im[it5][it3][it1] / isample;
		  R11f_re[it5][it3][it1]  = a * cos(phi3) - b * sin(phi3);
		  R11f_im[it5][it3][it1]  = a * sin(phi3) + b * cos(phi3);
		  a = 1 * mu_23_2 * mu_12_2 * mu_01_2 * R12_re[it5][it3][it1] / isample;
		  b = 1 * mu_23_2 * mu_12_2 * mu_01_2 * R12_im[it5][it3][it1] / isample;
		  R12f_re[it5][it3][it1]  =  a * cos(phi3) - b * sin(phi3);
		  R12f_im[it5][it3][it1]  =  a * sin(phi3) + b * cos(phi3);
		}

	      phi4 = -shift_w_1_2*dt*(it3-1) - shift_w_2_2*dt*(it5-1);
	      a = -2 * mu_12_2 * mu_12_2 * mu_01_2 * R13_re[it5][it3][it1] / isample;
	      b = -2 * mu_12_2 * mu_12_2 * mu_01_2 * R13_im[it5][it3][it1] / isample;
	      R13f_re[it5][it3][it1]  = a * cos(phi4) - b * sin(phi4);
	      R13f_im[it5][it3][it1]  = a * sin(phi4) - b * cos(phi4);
	      a = -2 * mu_12_2 * mu_12_2 * mu_01_2 * R14_re[it5][it3][it1] / isample;
	      b = -2 * mu_12_2 * mu_12_2 * mu_01_2 * R14_im[it5][it3][it1] / isample;
	      R14f_re[it5][it3][it1]  =  a * cos(phi4) - b * sin(phi4);
	      R14f_im[it5][it3][it1]  =  a * sin(phi4) + b * cos(phi4);

	      phi4 = shift_w_1_2*dt*(it3-1) - shift_w_2_2*dt*(it5-1);
	      a = -2 * mu_12_2 * mu_12_2 * mu_01_2 * R15_re[it5][it3][it1] / isample;
	      b = -2 * mu_12_2 * mu_12_2 * mu_01_2 * R15_im[it5][it3][it1] / isample;
	      R15f_re[it5][it3][it1]  = a * cos(phi4) - b * sin(phi4);
	      R15f_im[it5][it3][it1]  = a * sin(phi4) - b * cos(phi4);
	      a = -2 * mu_12_2 * mu_12_2 * mu_01_2 * R16_re[it5][it3][it1] / isample;
	      b = -2 * mu_12_2 * mu_12_2 * mu_01_2 * R16_im[it5][it3][it1] / isample;
	      R16f_re[it5][it3][it1]  =  a * cos(phi4) - b * sin(phi4);
	      R16f_im[it5][it3][it1]  =  a * sin(phi4) + b * cos(phi4);

	      phi5 = -shift_w_1_2*dt*(it3-1) - shift_w_2_1*dt*(it5-1);
	      a = 1 * mu_01_2 * mu_12_2 * mu_01_2 * R17_re[it5][it3][it1] / isample;
	      b = 1 * mu_01_2 * mu_12_2 * mu_01_2 * R17_im[it5][it3][it1] / isample;
	      R17f_re[it5][it3][it1]  = a * cos(phi5) - b * sin(phi5);
	      R17f_im[it5][it3][it1]  = a * sin(phi5) + b * cos(phi5);
	      a = 1 * mu_01_2 * mu_12_2 * mu_01_2 * R18_re[it5][it3][it1] / isample;
	      b = 1 * mu_12_2 * mu_12_2 * mu_01_2 * R18_im[it5][it3][it1] / isample;
	      R18f_re[it5][it3][it1]  =  a * cos(phi5) - b * sin(phi5);
	      R18f_im[it5][it3][it1]  =  a * sin(phi5) + b * cos(phi5);

	      phi5 = shift_w_1_2*dt*(it3-1) - shift_w_2_1*dt*(it5-1);
	      a = 1 * mu_01_2 * mu_12_2 * mu_01_2 * R19_re[it5][it3][it1] / isample;
	      b = 1 * mu_01_2 * mu_12_2 * mu_01_2 * R19_im[it5][it3][it1] / isample;
	      R19f_re[it5][it3][it1]  = a * cos(phi5) - b * sin(phi5);
	      R19f_im[it5][it3][it1]  = a * sin(phi5) + b * cos(phi5);
	      a = 1 * mu_01_2 * mu_12_2 * mu_01_2 * R20_re[it5][it3][it1] / isample;
	      b = 1 * mu_01_2 * mu_12_2 * mu_01_2 * R20_im[it5][it3][it1] / isample;
	      R20f_re[it5][it3][it1]  =  a * cos(phi5) - b * sin(phi5);
	      R20f_im[it5][it3][it1]  =  a * sin(phi5) + b * cos(phi5);
	      }

	      R1tot_re[it5][it3][it1] = R1f_re[it5][it3][it1] 
		+ R5f_re[it5][it3][it1] 
		+ R9f_re[it5][it3][it1] 
		+ R13f_re[it5][it3][it1] 
		+ R17f_re[it5][it3][it1];
	      R1tot_im[it5][it3][it1] = R1f_im[it5][it3][it1] 
		+ R5f_im[it5][it3][it1] 
		+ R9f_im[it5][it3][it1] 
		+ R13f_im[it5][it3][it1] 
		+ R17f_im[it5][it3][it1];
	      R2tot_re[it5][it3][it1] = R2f_re[it5][it3][it1] 
		+ R6f_re[it5][it3][it1] 
		+ R10f_re[it5][it3][it1] 
		+ R14f_re[it5][it3][it1] 
		+ R18f_re[it5][it3][it1];
	      R2tot_im[it5][it3][it1] = R2f_im[it5][it3][it1] 
		+ R6f_im[it5][it3][it1] 
		+ R10f_im[it5][it3][it1] 
		+ R14f_im[it5][it3][it1] 
		+ R18f_im[it5][it3][it1];
	      R3tot_re[it5][it3][it1] = R3f_re[it5][it3][it1] 
		+ R7f_re[it5][it3][it1] 
		+ R11f_re[it5][it3][it1] 
		+ R15f_re[it5][it3][it1] 
		+ R19f_re[it5][it3][it1];
	      R3tot_im[it5][it3][it1] = R3f_im[it5][it3][it1] 
		+ R7f_im[it5][it3][it1] 
		+ R11f_im[it5][it3][it1] 
		+ R15f_im[it5][it3][it1] 
		+ R19f_im[it5][it3][it1];
	      R4tot_re[it5][it3][it1] = R4f_re[it5][it3][it1] 
		+ R8f_re[it5][it3][it1] 
		+ R12f_re[it5][it3][it1] 
		+ R16f_re[it5][it3][it1] 
		+ R20f_re[it5][it3][it1];
	      R4tot_im[it5][it3][it1] = R4f_im[it5][it3][it1] 
		+ R8f_im[it5][it3][it1] 
		+ R12f_im[it5][it3][it1] 
		+ R16f_im[it5][it3][it1] 
		+ R20f_im[it5][it3][it1];
	    }
	}
    }
}

void freqTrajToR5( const char *base_name, float **t2_t4_pairs, const int n_t2_t4_pairs,const int nsteps,const int nprotons,float ***w_matrices, float ***x_matrices,float *mean_w,float *mean_mu ){
  
  //vars for force trajectory
  // unsigned long nsteps=-1, natoms, nmols, nprotons=-1;
  
  char *extension;
  
  //float **dw_matrix,**mu_matrix;
  const float wavenumbersToRadiansPerPs = 2*M_PI*2.9979e10*1e-12;
  /* Number of time Points for t1, t3, and t5 -- must be a power of 2!  */
  const int nt = globalArgs.nt; 
  /* number of time-steps to pre-integrate */
  const int ntint = globalArgs.ntint;                
  /* ps, time-steps output by GROMACS */
  const float dt = globalArgs.dt*ntint;//*wavenumbersToRadiansPerPs; 
  /* number of steps to skip between trajectory samples */
  const int nskip = globalArgs.nskip; 
  /* order of spectroscopy we will calculate */
  const int order = globalArgs.order;
  /* number of excited states */
  const int n_levels = globalArgs.n_levels;
  /* noncondon calc or not */
  const int flag_noncondon = globalArgs.flag_noncondon;

  float *dw; //a pointer into dw_matrix representing the frequency
             //trajectory of one proton
  float *mu;

  //vars for response function calcultions
  float *sdw;
  float **pdw1,**pdw2,**pdw3,**pdw4;
  float ***rdw1, ***rdw2, ***rdw3, ***rdw4;
  float ***rdw5, ***rdw6, ***rdw7, ***rdw8;
  float ***rdw9, ***rdw10,***rdw11,***rdw12;
  float ***rdw13,***rdw14,***rdw15,***rdw16;
  float ***rdw17,***rdw18,***rdw19,***rdw20;

  //3D response functions
  float ***R1_re, ***R1_im, ***R2_re, ***R2_im, ***R3_re, ***R3_im, ***R4_re, ***R4_im;
  float ***R5_re, ***R5_im, ***R6_re, ***R6_im, ***R7_re, ***R7_im, ***R8_re, ***R8_im;
  float ***R9_re, ***R9_im, ***R10_re,***R10_im,***R11_re,***R11_im,***R12_re,***R12_im;
  float ***R13_re,***R13_im,***R14_re,***R14_im,***R15_re,***R15_im,***R16_re,***R16_im;
  float ***R17_re,***R17_im,***R18_re,***R18_im,***R19_re,***R19_im,***R20_re,***R20_im;
  //2D response functions
  float **P1_re,**P1_im,**P2_re,**P2_im,**P3_re,**P3_im,**P4_re,**P4_im;
  //1D response functions
  float *S_re,*S_im;
  // for saving data
  float *Sf_re,*Sf_im;
  float **P1f_re,**P2f_re,**P2f_im,**P1f_im;
  float **P3f_re,**P4f_re,**P3f_im,**P4f_im;
  float **P1tot_re,**P1tot_im,**P2tot_re,**P2tot_im;
  float ***R1f_re,***R1f_im,***R2f_re,***R2f_im,***R3f_re,***R3f_im,***R4f_re,***R4f_im;
  float ***R5f_re, ***R5f_im, ***R6f_re, ***R6f_im, ***R7f_re, ***R7f_im, ***R8f_re, ***R8f_im;
  float ***R9f_re, ***R9f_im, ***R10f_re,***R10f_im,***R11f_re,***R11f_im,***R12f_re,***R12f_im;
  float ***R13f_re,***R13f_im,***R14f_re,***R14f_im,***R15f_re,***R15f_im,***R16f_re,***R16f_im;
  float ***R17f_re,***R17f_im,***R18f_re,***R18f_im,***R19f_re,***R19f_im,***R20f_re,***R20f_im;
  float ***R1tot_re,***R1tot_im,***R2tot_re,***R2tot_im,***R3tot_re,***R3tot_im,***R4tot_re,***R4tot_im;

  float t2,t4;
  
  int i,j,isample,i_level,it,it1,it3,it5,nt2,nt4,ntraject,nsamples,ipoptime,iproton;
  float *dwint1;          /* trajectories after preintegrating/coarse graining  */   
  float *dwint2;          /* trajectories after preintegrating/coarse graining  */   
  float *dwint3;          /* trajectories after preintegrating/coarse graining  */   
  float *muint1;          /* trajectories after preintegrating/coarse graining  */   
  float *muint2;          /* trajectories after preintegrating/coarse graining  */   
  float *muint3;          /* trajectories after preintegrating/coarse graining  */   
 
  float mean_w01;
  float mean_w12;
  float mean_w23;
  float mean_mu01;
  float mean_mu12;
  float mean_mu23;

  float mu2_mu1;
  float mu4_mu3_mu2_mu1;
  float mu6_mu5_mu4_mu3_mu2_mu1;
  double sum,sum2;

  /* DEBUG NONCONDON CALCS!!! */
  float mu_01_2=1,mu_12_2=1,mu_23_2=1;
  if (flag_noncondon==1){
    if (n_levels>=1) mu_01_2 = mean_mu[0]*mean_mu[0];
    if (n_levels>=2) mu_12_2 = mean_mu[1]*mean_mu[1];
    if (n_levels>=3) mu_23_2 = mean_mu[2]*mean_mu[2];
  }

  /*allocate arrays*/
  printf("Initialize Arrays\n");
  
  /* 3D peak 1 */
  R1_re=f3tensor(1,nt,1,nt,1,nt);
  R1_im=f3tensor(1,nt,1,nt,1,nt);
  R2_re=f3tensor(1,nt,1,nt,1,nt);
  R2_im=f3tensor(1,nt,1,nt,1,nt);
  R3_re=f3tensor(1,nt,1,nt,1,nt);
  R3_im=f3tensor(1,nt,1,nt,1,nt);
  R4_re=f3tensor(1,nt,1,nt,1,nt);
  R4_im=f3tensor(1,nt,1,nt,1,nt);
  /* 3D peak 2 */
  R5_re=f3tensor(1,nt,1,nt,1,nt);
  R5_im=f3tensor(1,nt,1,nt,1,nt);
  R6_re=f3tensor(1,nt,1,nt,1,nt);
  R6_im=f3tensor(1,nt,1,nt,1,nt);
  R7_re=f3tensor(1,nt,1,nt,1,nt);
  R7_im=f3tensor(1,nt,1,nt,1,nt);
  R8_re=f3tensor(1,nt,1,nt,1,nt);
  R8_im=f3tensor(1,nt,1,nt,1,nt);
  /* 3D peak 3 */
  R9_re=f3tensor(1,nt,1,nt,1,nt);
  R9_im=f3tensor(1,nt,1,nt,1,nt);
  R10_re=f3tensor(1,nt,1,nt,1,nt);
  R10_im=f3tensor(1,nt,1,nt,1,nt);
  R11_re=f3tensor(1,nt,1,nt,1,nt);
  R11_im=f3tensor(1,nt,1,nt,1,nt);
  R12_re=f3tensor(1,nt,1,nt,1,nt);
  R12_im=f3tensor(1,nt,1,nt,1,nt);
  /* 3D peak 4 */
  R13_re=f3tensor(1,nt,1,nt,1,nt);
  R13_im=f3tensor(1,nt,1,nt,1,nt);
  R14_re=f3tensor(1,nt,1,nt,1,nt);
  R14_im=f3tensor(1,nt,1,nt,1,nt);
  R15_re=f3tensor(1,nt,1,nt,1,nt);
  R15_im=f3tensor(1,nt,1,nt,1,nt);
  R16_re=f3tensor(1,nt,1,nt,1,nt);
  R16_im=f3tensor(1,nt,1,nt,1,nt);
  /* 3D peak 5 */
  R17_re=f3tensor(1,nt,1,nt,1,nt);
  R17_im=f3tensor(1,nt,1,nt,1,nt);
  R18_re=f3tensor(1,nt,1,nt,1,nt);
  R18_im=f3tensor(1,nt,1,nt,1,nt);
  R19_re=f3tensor(1,nt,1,nt,1,nt);
  R19_im=f3tensor(1,nt,1,nt,1,nt);
  R20_re=f3tensor(1,nt,1,nt,1,nt);
  R20_im=f3tensor(1,nt,1,nt,1,nt);
  /* 2D gb + se */
  P1_re=matrix(1,nt,1,nt);
  P1_im=matrix(1,nt,1,nt);
  P2_re=matrix(1,nt,1,nt);
  P2_im=matrix(1,nt,1,nt);
  /* 2D esa */
  P3_re=matrix(1,nt,1,nt);
  P3_im=matrix(1,nt,1,nt);
  P4_re=matrix(1,nt,1,nt);
  P4_im=matrix(1,nt,1,nt);
  /* 1D */
  S_re=vector(1,nt);
  S_im=vector(1,nt);
  
  /* normalized results and final output */
  /* 3D */
  R1f_re  = f3tensor(1,nt,1,nt,1,nt);
  R1f_im  = f3tensor(1,nt,1,nt,1,nt);
  R2f_re  = f3tensor(1,nt,1,nt,1,nt);
  R2f_im  = f3tensor(1,nt,1,nt,1,nt);
  R3f_re  = f3tensor(1,nt,1,nt,1,nt);
  R3f_im  = f3tensor(1,nt,1,nt,1,nt);
  R4f_re  = f3tensor(1,nt,1,nt,1,nt);
  R4f_im  = f3tensor(1,nt,1,nt,1,nt);
  R5f_re  = f3tensor(1,nt,1,nt,1,nt);
  R5f_im  = f3tensor(1,nt,1,nt,1,nt);
  R6f_re  = f3tensor(1,nt,1,nt,1,nt);
  R6f_im  = f3tensor(1,nt,1,nt,1,nt);
  R7f_re  = f3tensor(1,nt,1,nt,1,nt);
  R7f_im  = f3tensor(1,nt,1,nt,1,nt);
  R8f_re  = f3tensor(1,nt,1,nt,1,nt);
  R8f_im  = f3tensor(1,nt,1,nt,1,nt);
  R9f_re  = f3tensor(1,nt,1,nt,1,nt);
  R9f_im  = f3tensor(1,nt,1,nt,1,nt);
  R10f_re = f3tensor(1,nt,1,nt,1,nt);
  R10f_im = f3tensor(1,nt,1,nt,1,nt);
  R11f_re = f3tensor(1,nt,1,nt,1,nt);
  R11f_im = f3tensor(1,nt,1,nt,1,nt);
  R12f_re = f3tensor(1,nt,1,nt,1,nt);
  R12f_im = f3tensor(1,nt,1,nt,1,nt);
  R13f_re = f3tensor(1,nt,1,nt,1,nt);
  R13f_im = f3tensor(1,nt,1,nt,1,nt);
  R14f_re = f3tensor(1,nt,1,nt,1,nt);
  R14f_im = f3tensor(1,nt,1,nt,1,nt);
  R15f_re = f3tensor(1,nt,1,nt,1,nt);
  R15f_im = f3tensor(1,nt,1,nt,1,nt);
  R16f_re = f3tensor(1,nt,1,nt,1,nt);
  R16f_im = f3tensor(1,nt,1,nt,1,nt);
  R17f_re = f3tensor(1,nt,1,nt,1,nt);
  R17f_im = f3tensor(1,nt,1,nt,1,nt);
  R18f_re = f3tensor(1,nt,1,nt,1,nt);
  R18f_im = f3tensor(1,nt,1,nt,1,nt);
  R19f_re = f3tensor(1,nt,1,nt,1,nt);
  R19f_im = f3tensor(1,nt,1,nt,1,nt);
  R20f_re = f3tensor(1,nt,1,nt,1,nt);
  R20f_im = f3tensor(1,nt,1,nt,1,nt);
  R1tot_re  = f3tensor(1,nt,1,nt,1,nt);
  R1tot_im  = f3tensor(1,nt,1,nt,1,nt);
  R2tot_re  = f3tensor(1,nt,1,nt,1,nt);
  R2tot_im  = f3tensor(1,nt,1,nt,1,nt);
  R3tot_re  = f3tensor(1,nt,1,nt,1,nt);
  R3tot_im  = f3tensor(1,nt,1,nt,1,nt);
  R4tot_re  = f3tensor(1,nt,1,nt,1,nt);
  R4tot_im  = f3tensor(1,nt,1,nt,1,nt);


  /* 2D gb + se */
  P1f_re=matrix(1,nt,1,nt);
  P1f_im=matrix(1,nt,1,nt);
  P2f_re=matrix(1,nt,1,nt);
  P2f_im=matrix(1,nt,1,nt);
  /* 2D esa */
  P3f_re=matrix(1,nt,1,nt);
  P3f_im=matrix(1,nt,1,nt);
  P4f_re=matrix(1,nt,1,nt);
  P4f_im=matrix(1,nt,1,nt);
  /* 2D total spectrum */
  P1tot_re=matrix(1,nt,1,nt);
  P1tot_im=matrix(1,nt,1,nt);
  P2tot_re=matrix(1,nt,1,nt);
  P2tot_im=matrix(1,nt,1,nt);
  /* 1D */
  Sf_re=vector(1,nt);
  Sf_im=vector(1,nt);
 
  /* the time integrals to evaluate inside the exponent */
  sdw = vector(1,nt);
  pdw1 = matrix(1,nt,1,nt);
  pdw2 = matrix(1,nt,1,nt);
  pdw3 = matrix(1,nt,1,nt);
  pdw4 = matrix(1,nt,1,nt);
  rdw1  = f3tensor(1,nt,1,nt,1,nt);
  rdw2  = f3tensor(1,nt,1,nt,1,nt);
  rdw3  = f3tensor(1,nt,1,nt,1,nt);
  rdw4  = f3tensor(1,nt,1,nt,1,nt);
  rdw5  = f3tensor(1,nt,1,nt,1,nt);
  rdw6  = f3tensor(1,nt,1,nt,1,nt);
  rdw7  = f3tensor(1,nt,1,nt,1,nt);
  rdw8  = f3tensor(1,nt,1,nt,1,nt);
  rdw9  = f3tensor(1,nt,1,nt,1,nt);
  rdw10 = f3tensor(1,nt,1,nt,1,nt);
  rdw11 = f3tensor(1,nt,1,nt,1,nt);
  rdw12 = f3tensor(1,nt,1,nt,1,nt);
  rdw13 = f3tensor(1,nt,1,nt,1,nt);
  rdw14 = f3tensor(1,nt,1,nt,1,nt);
  rdw15 = f3tensor(1,nt,1,nt,1,nt);
  rdw16 = f3tensor(1,nt,1,nt,1,nt);
  rdw17 = f3tensor(1,nt,1,nt,1,nt);
  rdw18 = f3tensor(1,nt,1,nt,1,nt);
  rdw19 = f3tensor(1,nt,1,nt,1,nt);
  rdw20 = f3tensor(1,nt,1,nt,1,nt);
  
  /*
   *
   * Average dipole and frequencies. The frequencies make it so we can
   * do things in the rotating frame
   *
   */
  if (n_levels>=1) mean_mu01 = mean_mu[0];
  if (n_levels>=2) mean_mu12 = mean_mu[1];
  if (n_levels>=3) mean_mu23 = mean_mu[2];
  if (n_levels>=1) mean_w01 = mean_w[0];
  if (n_levels>=2) mean_w12 = mean_w[1];
  if (n_levels>=3) mean_w23 = mean_w[2];

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
      printf("Time delay %i of %i, t2 = %f t4 = %f\n",ipoptime,n_t2_t4_pairs,t2,t4);
      
      /*initialize arrays*/
      for(it1=1;it1<=nt;it1++)
	{
	  S_re[it1]=0;
	  S_im[it1]=0;
	  for(it3=1;it3<=nt;it3++)
	    {
	      P1_re[it3][it1]=0; /* R4 + R5 in book notation */
	      P1_im[it3][it1]=0;
	      P2_re[it3][it1]=0; /* R1 + R2 */
	      P2_im[it3][it1]=0;
	      P3_re[it3][it1]=0; /* R6 in book notation */
	      P3_im[it3][it1]=0;
	      P4_re[it3][it1]=0; /* R3 */
	      P4_im[it3][it1]=0;
	      for(it5=1;it5<=nt;it5++)
		{
		  //
		  R1_re[it5][it3][it1]=0;
		  R1_im[it5][it3][it1]=0;
		  R2_re[it5][it3][it1]=0;
		  R2_im[it5][it3][it1]=0;
		  R3_re[it5][it3][it1]=0;
		  R3_im[it5][it3][it1]=0;
		  R4_re[it5][it3][it1]=0;
		  R4_im[it5][it3][it1]=0;
		  //
		  R5_re[it5][it3][it1]=0;
		  R5_im[it5][it3][it1]=0;
		  R6_re[it5][it3][it1]=0;
		  R6_im[it5][it3][it1]=0;
		  R7_re[it5][it3][it1]=0;
		  R7_im[it5][it3][it1]=0;
		  R8_re[it5][it3][it1]=0;
		  R8_im[it5][it3][it1]=0;
		  //
		  R9_re[it5][it3][it1]=0;
		  R9_im[it5][it3][it1]=0;
		  R10_re[it5][it3][it1]=0;
		  R10_im[it5][it3][it1]=0;
		  R11_re[it5][it3][it1]=0;
		  R11_im[it5][it3][it1]=0;
		  R12_re[it5][it3][it1]=0;
		  R12_im[it5][it3][it1]=0;
		  //
		  R13_re[it5][it3][it1]=0;
		  R13_im[it5][it3][it1]=0;
		  R14_re[it5][it3][it1]=0;
		  R14_im[it5][it3][it1]=0;
		  R15_re[it5][it3][it1]=0;
		  R15_im[it5][it3][it1]=0;
		  R16_re[it5][it3][it1]=0;
		  R16_im[it5][it3][it1]=0;
		  //
		  R17_re[it5][it3][it1]=0;
		  R17_im[it5][it3][it1]=0;
		  R18_re[it5][it3][it1]=0;
		  R18_im[it5][it3][it1]=0;
		  R19_re[it5][it3][it1]=0;
		  R19_im[it5][it3][it1]=0;
		  R20_re[it5][it3][it1]=0;
		  R20_im[it5][it3][it1]=0;
		}
	    }
	}
      
      if (asprintf(&extension,"_t2_%d_t4_%d.dat",(int) t2,(int) t4) < 0) nrerror("failed to write string");
      printf("Extension: %s\n",extension);
      
      /* nt2 and nt4 are the number of coarse grained points
       * during population times in trajectory */
      nt4=round(t4/ntint);
      nt2=round(t2/ntint);
     
      /* number of points in the coarse-grained trajectory */
      ntraject = nt; /* default for 1D */
      if (order==3) ntraject = 2*nt+nt2; 
      if (order==5) ntraject = 3*nt+nt2+nt4; 

      /* frequency trajectories averaged over ntint for each level */
      dwint1=vector(1,ntraject);
      dwint2=vector(1,ntraject);
      dwint3=vector(1,ntraject);
      muint1=vector(1,ntraject);
      muint2=vector(1,ntraject);
      muint3=vector(1,ntraject);
      
      /* calculate the maximum number of dw trajectories in the total
       * time evolution as <nsamples> */
      if (nsteps<ntraject*ntint)
	  nrerror("ERROR: trajectory is too short for response function calculation!");
      nsamples = floor( (nsteps-ntraject*ntint)/nskip );
      printf("The number of samples for this set of delays is %6.4i\n",nsamples);
      
      /* convert units */
      for (i_level=0;i_level<n_levels;i_level++){
	sum  = 0;
	sum2 = 0;
	for (i=1;i<=nprotons;i++)
	  for (j=1;j<=nsteps;j++){
	    sum+=(double)w_matrices[i_level][i][j];
	    sum2+=(double)w_matrices[i_level][i][j]*w_matrices[i_level][i][j];
	    w_matrices[i_level][i][j] *= wavenumbersToRadiansPerPs;
	  }
	sum/=nprotons*nsteps;
	sum2/=nprotons*nsteps;
	sum2 = sqrt(sum2);

	printf("mean w[%i] %f cm-1 ",i_level,mean_w[i_level]);
	mean_w[i_level] *= wavenumbersToRadiansPerPs;
	printf("%f rad/ps\n",mean_w[i_level]);
	printf("mean dw[%i] %f cm-1 ",i_level,sum);
	sum *= wavenumbersToRadiansPerPs;
	printf("%f rad/ps\n",sum);
	printf("mean dw^2[%i] %f cm-1 ",i_level,sum2);
	sum2 *= wavenumbersToRadiansPerPs;
	printf("%f rad/ps\n",sum2);
	printf("mean mu[%i] %f (units?)\n",i_level,mean_mu[i_level]);
      }
      printf("dt = %f ps\n",dt);

      //loop over each proton <iproton>
      printf("Start processing frequency trajectories calculating response functions.\n");
      for(iproton=1;iproton<=nprotons;iproton++) 
	for(isample=1;isample<=nsamples;isample++)
	  {
	    if((isample+nsamples*(iproton-1))%5000==0) 
	      printf("trajectory #%d\n",(isample+nsamples*(iproton-1)));
        
	    /* point dw at the appropriate row of dw_matrix
	     * offset by an amount (nskip) to make this an independent,
	     * uncorrelated trajectory */

	      i_level = 0; /* 01 transition */
	      dw=&w_matrices[i_level][iproton][1+(isample-1)*nskip];
	      mu=&x_matrices[i_level][iproton][1+(isample-1)*nskip];
	      
	      /* preintegrate the frequency and convert units */
	      for(it=1;it<=ntraject;it++)
		{
		  dwint1[it]=(dw[(it-1)*ntint+1]+dw[it*ntint+1])/2;
		  for(j=2;j<=ntint;j++)
		    dwint1[it]+=dw[(it-1)*ntint+j];
		  
		  /* the code above has an effective length of
		   * ntint+1 (where matlab trapz() has a length of
		   * ntint-1!  but hte denom has to be ntint in the
		   * denom to get agreement with matlab???
		   */
		  dwint1[it]*=dt/ntint;
		}
	      /* preintegrate the dipole */
	      for(it=1;it<=ntraject;it++)
		{
		  muint1[it]=(mu[(it-1)*ntint+1]+mu[it*ntint+1])/2;
		  for(j=2;j<=ntint;j++)
		    muint1[it]+=mu[(it-1)*ntint+j];
		  muint1[it]/=(ntint);
		}
	      
	      /* output the first freq trajectory for debugging purposes */
	      if (DEBUG_LEVEL>=1){
		if (iproton==1 && isample==1){
		  printf("dwint1 = \n");
		  for (i=1;i<=ntraject;i++)
		    printf("%f\n",dwint1[i]);
		}
	      }
	      
	      /* do the same preintegration for the higher states */
	      if (n_levels>=2){
		i_level = 1; /* the 12 transition */
		dw=&w_matrices[i_level][iproton][1+(isample-1)*nskip];
		mu=&x_matrices[i_level][iproton][1+(isample-1)*nskip];
		
		/* preintegrate the frequency and convert units */
		for(it=1;it<=ntraject;it++)
		  {
		    dwint2[it]=(dw[(it-1)*ntint+1]+dw[it*ntint+1])/2;
		    for(j=2;j<=ntint;j++)
		      dwint2[it]+=dw[(it-1)*ntint+j];
		    
		    /* the code above has an effective length of
		     * ntint+1 (where matlab trapz() has a length of
		     * ntint-1!  but hte denom has to be ntint in the
		     * denom to get agreement with matlab???
		     */
		    dwint2[it]*=dt/ntint;
		  }
		
		for(it=1;it<=ntraject;it++)
		  {
		    muint2[it]=(mu[(it-1)*ntint+1]+mu[it*ntint+1])/2;
		    for(j=2;j<=ntint;j++)
		      muint2[it]+=mu[(it-1)*ntint+j];
		    muint2[it]/=(ntint);
		  }
	      } //end if (n_levels>=2) 
	      /* do the same preintegration for the higher states */
	      if (n_levels>=3){
		i_level = 2; /* the 23 transition */
		dw=&w_matrices[i_level][iproton][1+(isample-1)*nskip];
		mu=&x_matrices[i_level][iproton][1+(isample-1)*nskip];
		
		/* preintegrate the frequency and convert units */
		for(it=1;it<=ntraject;it++)
		  {
		    dwint3[it]=(dw[(it-1)*ntint+1]+dw[it*ntint+1])/2;
		    for(j=2;j<=ntint;j++)
		      dwint3[it]+=dw[(it-1)*ntint+j];
		    
		    /* the code above has an effective length of
		     * ntint+1 (where matlab trapz() has a length of
		     * ntint-1!  but hte denom has to be ntint in the
		     * denom to get agreement with matlab???
		     */
		    dwint3[it]*=dt/ntint;
		  }
		
		for(it=1;it<=ntraject;it++)
		  {
		    muint3[it]=(mu[(it-1)*ntint+1]+mu[it*ntint+1])/2;
		    for(j=2;j<=ntint;j++)
		      muint3[it]+=mu[(it-1)*ntint+j];
		    muint3[it]/=(ntint);
		  }
	      } //end if (n_levels>=3) 
		
	    if (flag_noncondon==0){

	      //	      printf("Condon calculation\n");

	      /*
	       * Calculate nonlinear polarisation
	       *
	       */
	      
	      sdw[1]=0;
	      pdw1[1][1]=0;
	      pdw2[1][1]=0;
	      pdw3[1][1]=0;
	      pdw4[1][1]=0;
	      rdw1[1][1][1]  = 0;
	      rdw2[1][1][1]  = 0;
	      rdw3[1][1][1]  = 0;
	      rdw4[1][1][1]  = 0;
	      rdw5[1][1][1]  = 0;
	      rdw6[1][1][1]  = 0;
	      rdw7[1][1][1]  = 0;
	      rdw8[1][1][1]  = 0;
	      rdw9[1][1][1]  = 0;
	      rdw10[1][1][1] = 0;
	      rdw11[1][1][1] = 0;
	      rdw12[1][1][1] = 0;
	      rdw13[1][1][1] = 0;
	      rdw14[1][1][1] = 0;
	      rdw15[1][1][1] = 0;
	      rdw16[1][1][1] = 0;
	      rdw17[1][1][1] = 0;
	      rdw18[1][1][1] = 0;
	      rdw19[1][1][1] = 0;
	      rdw20[1][1][1] = 0;
	      for(it1=1;it1<=nt;it1++)
		{
		  if (it1<nt)
		    {
		      sdw[it1+1]=sdw[it1]-dwint1[it1];
		      pdw1[1][it1+1]=pdw1[1][it1]+dwint1[it1];
		      pdw2[1][it1+1]=pdw2[1][it1]-dwint1[it1];
		      pdw3[1][it1+1]=pdw3[1][it1]+dwint1[it1];
		      pdw4[1][it1+1]=pdw4[1][it1]-dwint1[it1];
		      /* peak 1 */
		      rdw1[1][1][it1+1]  = rdw1[1][1][it1]  - dwint1[it1];
		      rdw2[1][1][it1+1]  = rdw2[1][1][it1]  + dwint1[it1];
		      rdw3[1][1][it1+1]  = rdw3[1][1][it1]  - dwint1[it1];
		      rdw4[1][1][it1+1]  = rdw4[1][1][it1]  + dwint1[it1];
		      /* peak 2 */
		      rdw5[1][1][it1+1]  = rdw5[1][1][it1]  - dwint1[it1];
		      rdw6[1][1][it1+1]  = rdw6[1][1][it1]  + dwint1[it1];
		      rdw7[1][1][it1+1]  = rdw7[1][1][it1]  - dwint1[it1];
		      rdw8[1][1][it1+1]  = rdw8[1][1][it1]  + dwint1[it1];
		      /* peak 3 */
		      rdw9[1][1][it1+1]  = rdw9[1][1][it1]  - dwint1[it1];
		      rdw10[1][1][it1+1] = rdw10[1][1][it1] + dwint1[it1];
		      rdw11[1][1][it1+1] = rdw11[1][1][it1] - dwint1[it1];
		      rdw12[1][1][it1+1] = rdw12[1][1][it1] + dwint1[it1];
		      /* peak 4 */
		      rdw13[1][1][it1+1] = rdw13[1][1][it1] - dwint1[it1];
		      rdw14[1][1][it1+1] = rdw14[1][1][it1] + dwint1[it1];
		      rdw15[1][1][it1+1] = rdw15[1][1][it1] - dwint1[it1];
		      rdw16[1][1][it1+1] = rdw16[1][1][it1] + dwint1[it1];
		      /* peak 5 */
		      rdw17[1][1][it1+1] = rdw17[1][1][it1] - dwint1[it1];
		      rdw18[1][1][it1+1] = rdw18[1][1][it1] + dwint1[it1];
		      rdw19[1][1][it1+1] = rdw19[1][1][it1] - dwint1[it1];
		      rdw20[1][1][it1+1] = rdw20[1][1][it1] + dwint1[it1];
		    }
		  S_re[it1]+=cos(sdw[it1]);
		  S_im[it1]+=sin(sdw[it1]);
		  if (order>=3) 
		    for(it3=1;it3<=nt;it3++)
		      {
			if (it3<nt)
			  {
			    pdw1[it3+1][it1] = pdw1[it3][it1] - dwint1[it3+nt2+it1-1];
			    pdw2[it3+1][it1] = pdw2[it3][it1] - dwint1[it3+nt2+it1-1];
			    /* peak 1 */
			    rdw1[1][it3+1][it1]  = rdw1[1][it3][it1]  - dwint1[it3+nt2+it1-1];
			    rdw2[1][it3+1][it1]  = rdw2[1][it3][it1]  - dwint1[it3+nt2+it1-1];
			    rdw3[1][it3+1][it1]  = rdw3[1][it3][it1]  + dwint1[it3+nt2+it1-1];
			    rdw4[1][it3+1][it1]  = rdw4[1][it3][it1]  + dwint1[it3+nt2+it1-1];
			    if (n_levels>=2){
			      pdw3[it3+1][it1]=pdw3[it3][it1] - dwint2[it3+nt2+it1-1];
			      pdw4[it3+1][it1]=pdw4[it3][it1] - dwint2[it3+nt2+it1-1];
			      /* peak 2 */
			      rdw5[1][it3+1][it1]  = rdw5[1][it3][it1]  - dwint1[it3+nt2+it1-1];
			      rdw6[1][it3+1][it1]  = rdw6[1][it3][it1]  - dwint1[it3+nt2+it1-1];
			      rdw7[1][it3+1][it1]  = rdw7[1][it3][it1]  + dwint1[it3+nt2+it1-1];
			      rdw8[1][it3+1][it1]  = rdw8[1][it3][it1]  + dwint1[it3+nt2+it1-1];
			      if (n_levels>=3){
				/* peak 3 */
				rdw9[1][it3+1][it1]  = rdw9[1][it3][it1]  - dwint2[it3+nt2+it1-1];
				rdw10[1][it3+1][it1] = rdw10[1][it3][it1] - dwint2[it3+nt2+it1-1];
				rdw11[1][it3+1][it1] = rdw11[1][it3][it1] + dwint2[it3+nt2+it1-1];
				rdw12[1][it3+1][it1] = rdw12[1][it3][it1] + dwint2[it3+nt2+it1-1];
			      } /* end n_levels<=3 */
			      /* peak 4 */
			      rdw13[1][it3+1][it1] = rdw13[1][it3][it1] - dwint2[it3+nt2+it1-1];
			      rdw14[1][it3+1][it1] = rdw14[1][it3][it1] - dwint2[it3+nt2+it1-1];
			      rdw15[1][it3+1][it1] = rdw15[1][it3][it1] + dwint2[it3+nt2+it1-1];
			      rdw16[1][it3+1][it1] = rdw16[1][it3][it1] + dwint2[it3+nt2+it1-1];
			      /* peak 5 */
			      rdw17[1][it3+1][it1] = rdw17[1][it3][it1] - dwint2[it3+nt2+it1-1];
			      rdw18[1][it3+1][it1] = rdw18[1][it3][it1] - dwint2[it3+nt2+it1-1];
			      rdw19[1][it3+1][it1] = rdw19[1][it3][it1] + dwint2[it3+nt2+it1-1];
			      rdw20[1][it3+1][it1] = rdw20[1][it3][it1] + dwint2[it3+nt2+it1-1];
			    } /* end n_levels>=2 */
			  } /* end it<nt */
			P1_re[it3][it1]+=cos(pdw1[it3][it1]);
			P1_im[it3][it1]+=sin(pdw1[it3][it1]);
			P2_re[it3][it1]+=cos(pdw2[it3][it1]);
			P2_im[it3][it1]+=sin(pdw2[it3][it1]);
			if (n_levels>=2){
			  P3_re[it3][it1]+=cos(pdw3[it3][it1]);
			  P3_im[it3][it1]+=sin(pdw3[it3][it1]);
			  P4_re[it3][it1]+=cos(pdw4[it3][it1]);
			  P4_im[it3][it1]+=sin(pdw4[it3][it1]);
			}
			if (order>=5)
			  for(it5=1;it5<=nt;it5++)
			    {
			      if (it5<nt)
				{
				  /* peak 1 */
				  rdw1[it5+1][it3][it1]  = rdw1[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				  rdw2[it5+1][it3][it1]  = rdw2[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				  rdw3[it5+1][it3][it1]  = rdw3[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				  rdw4[it5+1][it3][it1]  = rdw4[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				  if (n_levels>=2){
				    /* peak 2 */
				    rdw5[it5+1][it3][it1]  = rdw5[it5][it3][it1]  - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw6[it5+1][it3][it1]  = rdw6[it5][it3][it1]  - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw7[it5+1][it3][it1]  = rdw7[it5][it3][it1]  - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw8[it5+1][it3][it1]  = rdw8[it5][it3][it1]  - dwint2[it5+nt4+it3+nt2+it1-2];
				    if (n_levels>=3){
				      /* peak 3 */
				      rdw9[it5+1][it3][it1]  = rdw9[it5][it3][it1]  - dwint3[it5+nt4+it3+nt2+it1-2];
				      rdw10[it5+1][it3][it1] = rdw10[it5][it3][it1] - dwint3[it5+nt4+it3+nt2+it1-2];
				      rdw11[it5+1][it3][it1] = rdw11[it5][it3][it1] - dwint3[it5+nt4+it3+nt2+it1-2];
				      rdw12[it5+1][it3][it1] = rdw12[it5][it3][it1] - dwint3[it5+nt4+it3+nt2+it1-2];
				    } /* end n_levels<=3 */

				    /* peak 4 */
				    rdw13[it5+1][it3][it1] = rdw13[it5][it3][it1] - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw14[it5+1][it3][it1] = rdw14[it5][it3][it1] - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw15[it5+1][it3][it1] = rdw15[it5][it3][it1] - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw16[it5+1][it3][it1] = rdw16[it5][it3][it1] - dwint2[it5+nt4+it3+nt2+it1-2];
				    /* peak 5 */
				    rdw17[it5+1][it3][it1] = rdw17[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				    rdw18[it5+1][it3][it1] = rdw18[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				    rdw19[it5+1][it3][it1] = rdw19[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				    rdw20[it5+1][it3][it1] = rdw20[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				  } /* end n_levels<=2 */
				} /* end it5<nt */
			      /* peak 1 */
			      R1_re[it5][it3][it1]  += cos(rdw1[it5][it3][it1]);
			      R1_im[it5][it3][it1]  += sin(rdw1[it5][it3][it1]);
			      R2_re[it5][it3][it1]  += cos(rdw2[it5][it3][it1]);
			      R2_im[it5][it3][it1]  += sin(rdw2[it5][it3][it1]);
			      R3_re[it5][it3][it1]  += cos(rdw3[it5][it3][it1]);
			      R3_im[it5][it3][it1]  += sin(rdw3[it5][it3][it1]);
			      R4_re[it5][it3][it1]  += cos(rdw4[it5][it3][it1]);
			      R4_im[it5][it3][it1]  += sin(rdw4[it5][it3][it1]);
			      if (n_levels>=2){
				/* peak 2 */
				R5_re[it5][it3][it1]  += cos(rdw5[it5][it3][it1]);
				R5_im[it5][it3][it1]  += sin(rdw5[it5][it3][it1]);
				R6_re[it5][it3][it1]  += cos(rdw6[it5][it3][it1]);
				R6_im[it5][it3][it1]  += sin(rdw6[it5][it3][it1]);
				R7_re[it5][it3][it1]  += cos(rdw7[it5][it3][it1]);
				R7_im[it5][it3][it1]  += sin(rdw7[it5][it3][it1]);
				R8_re[it5][it3][it1]  += cos(rdw8[it5][it3][it1]);
				R8_im[it5][it3][it1]  += sin(rdw8[it5][it3][it1]);
				if (n_levels>=3){
				  R9_re[it5][it3][it1]  += cos(rdw9[it5][it3][it1]);
				  R9_im[it5][it3][it1]  += sin(rdw9[it5][it3][it1]);
				  R10_re[it5][it3][it1] += cos(rdw10[it5][it3][it1]);
				  R10_im[it5][it3][it1] += sin(rdw10[it5][it3][it1]);
				  R11_re[it5][it3][it1] += cos(rdw11[it5][it3][it1]);
				  R11_im[it5][it3][it1] += sin(rdw11[it5][it3][it1]);
				  R12_re[it5][it3][it1] += cos(rdw12[it5][it3][it1]);
				  R12_im[it5][it3][it1] += sin(rdw12[it5][it3][it1]);
				} /* end n_levels >= 3 */
				/* peak 4 */
				R13_re[it5][it3][it1] += cos(rdw13[it5][it3][it1]);
				R13_im[it5][it3][it1] += sin(rdw13[it5][it3][it1]);
				R14_re[it5][it3][it1] += cos(rdw14[it5][it3][it1]);
				R14_im[it5][it3][it1] += sin(rdw14[it5][it3][it1]);
				R15_re[it5][it3][it1] += cos(rdw15[it5][it3][it1]);
				R15_im[it5][it3][it1] += sin(rdw15[it5][it3][it1]);
				R16_re[it5][it3][it1] += cos(rdw16[it5][it3][it1]);
				R16_im[it5][it3][it1] += sin(rdw16[it5][it3][it1]);
				/* peak 5 */
				R17_re[it5][it3][it1] += cos(rdw17[it5][it3][it1]);
				R17_im[it5][it3][it1] += sin(rdw17[it5][it3][it1]);
				R18_re[it5][it3][it1] += cos(rdw18[it5][it3][it1]);
				R18_im[it5][it3][it1] += sin(rdw18[it5][it3][it1]);
				R19_re[it5][it3][it1] += cos(rdw19[it5][it3][it1]);
				R19_im[it5][it3][it1] += sin(rdw19[it5][it3][it1]);
				R20_re[it5][it3][it1] += cos(rdw20[it5][it3][it1]);
				R20_im[it5][it3][it1] += sin(rdw20[it5][it3][it1]);
			      } /* end n_levels >= 2 */
			    } /* end order >= 5 and it5 */
		      } /* end order >=3 and it3 */
		} /* end for it1 */
	    }/* end if condon */
	    
	    if (flag_noncondon==1){
	      // non condon code goes here

	      //	      printf("Non-condon calculation\n");

	      sdw[1]=0;
	      pdw1[1][1]=0;
	      pdw2[1][1]=0;
	      pdw3[1][1]=0;
	      pdw4[1][1]=0;
	      rdw1[1][1][1]  = 0;
	      rdw2[1][1][1]  = 0;
	      rdw3[1][1][1]  = 0;
	      rdw4[1][1][1]  = 0;
	      rdw5[1][1][1]  = 0;
	      rdw6[1][1][1]  = 0;
	      rdw7[1][1][1]  = 0;
	      rdw8[1][1][1]  = 0;
	      rdw9[1][1][1]  = 0;
	      rdw10[1][1][1] = 0;
	      rdw11[1][1][1] = 0;
	      rdw12[1][1][1] = 0;
	      rdw13[1][1][1] = 0;
	      rdw14[1][1][1] = 0;
	      rdw15[1][1][1] = 0;
	      rdw16[1][1][1] = 0;
	      rdw17[1][1][1] = 0;
	      rdw18[1][1][1] = 0;
	      rdw19[1][1][1] = 0;
	      rdw20[1][1][1] = 0;

		for(it1=1;it1<=nt;it1++){
		    if (it1<nt)
		      {
			sdw[it1+1]=sdw[it1]-dwint1[it1];
			pdw1[1][it1+1]=pdw1[1][it1]+dwint1[it1];
			pdw2[1][it1+1]=pdw2[1][it1]-dwint1[it1];
			pdw3[1][it1+1]=pdw3[1][it1]+dwint1[it1];
			pdw4[1][it1+1]=pdw4[1][it1]-dwint1[it1];
		      /* peak 1 */
		      rdw1[1][1][it1+1]  = rdw1[1][1][it1]  - dwint1[it1];
		      rdw2[1][1][it1+1]  = rdw2[1][1][it1]  + dwint1[it1];
		      rdw3[1][1][it1+1]  = rdw3[1][1][it1]  - dwint1[it1];
		      rdw4[1][1][it1+1]  = rdw4[1][1][it1]  + dwint1[it1];
		      /* peak 2 */
		      rdw5[1][1][it1+1]  = rdw5[1][1][it1]  - dwint1[it1];
		      rdw6[1][1][it1+1]  = rdw6[1][1][it1]  + dwint1[it1];
		      rdw7[1][1][it1+1]  = rdw7[1][1][it1]  - dwint1[it1];
		      rdw8[1][1][it1+1]  = rdw8[1][1][it1]  + dwint1[it1];
		      /* peak 3 */
		      rdw9[1][1][it1+1]  = rdw9[1][1][it1]  - dwint1[it1];
		      rdw10[1][1][it1+1] = rdw10[1][1][it1] + dwint1[it1];
		      rdw11[1][1][it1+1] = rdw11[1][1][it1] - dwint1[it1];
		      rdw12[1][1][it1+1] = rdw12[1][1][it1] + dwint1[it1];
		      /* peak 4 */
		      rdw13[1][1][it1+1] = rdw13[1][1][it1] - dwint1[it1];
		      rdw14[1][1][it1+1] = rdw14[1][1][it1] + dwint1[it1];
		      rdw15[1][1][it1+1] = rdw15[1][1][it1] - dwint1[it1];
		      rdw16[1][1][it1+1] = rdw16[1][1][it1] + dwint1[it1];
		      /* peak 5 */
		      rdw17[1][1][it1+1] = rdw17[1][1][it1] - dwint1[it1];
		      rdw18[1][1][it1+1] = rdw18[1][1][it1] + dwint1[it1];
		      rdw19[1][1][it1+1] = rdw19[1][1][it1] - dwint1[it1];
		      rdw20[1][1][it1+1] = rdw20[1][1][it1] + dwint1[it1];
		      }
		    mu2_mu1 = muint1[it1]*muint1[1];
		    S_re[it1]+=mu2_mu1*cos(sdw[it1]);
		    S_im[it1]+=mu2_mu1*sin(sdw[it1]);
		    if (order>=3) 
		      for(it3=1;it3<=nt;it3++)
			{
			  if (it3<nt)
			    {
			      pdw1[it3+1][it1] = pdw1[it3][it1] - dwint1[it3+nt2+it1-1]; //should this be it3+nt+it1-1???
			      pdw2[it3+1][it1] = pdw2[it3][it1] - dwint1[it3+nt2+it1-1];
			      /* peak 1 */
			      rdw1[1][it3+1][it1]  = rdw1[1][it3][it1]  - dwint1[it3+nt2+it1-1];
			      rdw2[1][it3+1][it1]  = rdw2[1][it3][it1]  - dwint1[it3+nt2+it1-1];
			      rdw3[1][it3+1][it1]  = rdw3[1][it3][it1]  + dwint1[it3+nt2+it1-1];
			      rdw4[1][it3+1][it1]  = rdw4[1][it3][it1]  + dwint1[it3+nt2+it1-1];
			      if (n_levels>=2){
				pdw3[it3+1][it1]=pdw3[it3][it1]-dwint2[it3+nt2+it1-1];
				pdw4[it3+1][it1]=pdw4[it3][it1]-dwint2[it3+nt2+it1-1];
			      /* peak 2 */
			      rdw5[1][it3+1][it1]  = rdw5[1][it3][it1]  - dwint1[it3+nt2+it1-1];
			      rdw6[1][it3+1][it1]  = rdw6[1][it3][it1]  - dwint1[it3+nt2+it1-1];
			      rdw7[1][it3+1][it1]  = rdw7[1][it3][it1]  + dwint1[it3+nt2+it1-1];
			      rdw8[1][it3+1][it1]  = rdw8[1][it3][it1]  + dwint1[it3+nt2+it1-1];
			      if (n_levels>=3){
				/* peak 3 */
				rdw9[1][it3+1][it1]  = rdw9[1][it3][it1]  - dwint2[it3+nt2+it1-1];
				rdw10[1][it3+1][it1] = rdw10[1][it3][it1] - dwint2[it3+nt2+it1-1];
				rdw11[1][it3+1][it1] = rdw11[1][it3][it1] + dwint2[it3+nt2+it1-1];
				rdw12[1][it3+1][it1] = rdw12[1][it3][it1] + dwint2[it3+nt2+it1-1];
			      } /* end n_levels>=3 */
			      /* peak 4 */
			      rdw13[1][it3+1][it1] = rdw13[1][it3][it1] - dwint2[it3+nt2+it1-1];
			      rdw14[1][it3+1][it1] = rdw14[1][it3][it1] - dwint2[it3+nt2+it1-1];
			      rdw15[1][it3+1][it1] = rdw15[1][it3][it1] + dwint2[it3+nt2+it1-1];
			      rdw16[1][it3+1][it1] = rdw16[1][it3][it1] + dwint2[it3+nt2+it1-1];
			      /* peak 5 */
			      rdw17[1][it3+1][it1] = rdw17[1][it3][it1] - dwint2[it3+nt2+it1-1];
			      rdw18[1][it3+1][it1] = rdw18[1][it3][it1] - dwint2[it3+nt2+it1-1];
			      rdw19[1][it3+1][it1] = rdw19[1][it3][it1] + dwint2[it3+nt2+it1-1];
			      rdw20[1][it3+1][it1] = rdw20[1][it3][it1] + dwint2[it3+nt2+it1-1];
			      } /* end n_levels>=2 */
			    }
			  
			  mu4_mu3_mu2_mu1 = muint1[it1+nt2+it3-1]*muint1[it1+nt2]*mu2_mu1;
			  
			  P1_re[it3][it1]+=mu4_mu3_mu2_mu1*cos(pdw1[it3][it1]);
			  P1_im[it3][it1]+=mu4_mu3_mu2_mu1*sin(pdw1[it3][it1]);
			  P2_re[it3][it1]+=mu4_mu3_mu2_mu1*cos(pdw2[it3][it1]);
			  P2_im[it3][it1]+=mu4_mu3_mu2_mu1*sin(pdw2[it3][it1]);
			  if (n_levels>=2){
			    mu4_mu3_mu2_mu1 = muint2[it1+nt2+it3-1]*muint2[it1+nt2]*mu2_mu1;
			    P3_re[it3][it1]+=mu4_mu3_mu2_mu1*cos(pdw3[it3][it1]);
			    P3_im[it3][it1]+=mu4_mu3_mu2_mu1*sin(pdw3[it3][it1]);
			    P4_re[it3][it1]+=mu4_mu3_mu2_mu1*cos(pdw4[it3][it1]);
			    P4_im[it3][it1]+=mu4_mu3_mu2_mu1*sin(pdw4[it3][it1]);
			  }
			  if (order>=5)
			  for(it5=1;it5<=nt;it5++)
			    {
			      if (it5<nt)
				{
				  /* peak 1 */
				  rdw1[it5+1][it3][it1]  = rdw1[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				  rdw2[it5+1][it3][it1]  = rdw2[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				  rdw3[it5+1][it3][it1]  = rdw3[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				  rdw4[it5+1][it3][it1]  = rdw4[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				  if (n_levels>=2){
				    /* peak 2 */
				    rdw5[it5+1][it3][it1]  = rdw5[it5][it3][it1]  - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw6[it5+1][it3][it1]  = rdw6[it5][it3][it1]  - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw7[it5+1][it3][it1]  = rdw7[it5][it3][it1]  - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw8[it5+1][it3][it1]  = rdw8[it5][it3][it1]  - dwint2[it5+nt4+it3+nt2+it1-2];
				    if (n_levels>=3){
				      /* peak 3 */
				      rdw9[it5+1][it3][it1]  = rdw9[it5][it3][it1]  - dwint3[it5+nt4+it3+nt2+it1-2];
				      rdw10[it5+1][it3][it1] = rdw10[it5][it3][it1] - dwint3[it5+nt4+it3+nt2+it1-2];
				      rdw11[it5+1][it3][it1] = rdw11[it5][it3][it1] - dwint3[it5+nt4+it3+nt2+it1-2];
				      rdw12[it5+1][it3][it1] = rdw12[it5][it3][it1] - dwint3[it5+nt4+it3+nt2+it1-2];
				    } /* end n_levels<=3 */

				    /* peak 4 */
				    rdw13[it5+1][it3][it1] = rdw13[it5][it3][it1] - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw14[it5+1][it3][it1] = rdw14[it5][it3][it1] - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw15[it5+1][it3][it1] = rdw15[it5][it3][it1] - dwint2[it5+nt4+it3+nt2+it1-2];
				    rdw16[it5+1][it3][it1] = rdw16[it5][it3][it1] - dwint2[it5+nt4+it3+nt2+it1-2];
				    /* peak 5 */
				    rdw17[it5+1][it3][it1] = rdw17[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				    rdw18[it5+1][it3][it1] = rdw18[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				    rdw19[it5+1][it3][it1] = rdw19[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				    rdw20[it5+1][it3][it1] = rdw20[it5][it3][it1] - dwint1[it5+nt4+it3+nt2+it1-2];
				  } /* end n_levels<=2 */
				} /* end it5<nt */

			      mu6_mu5_mu4_mu3_mu2_mu1 = muint1[it1+nt2+it3+nt4+it5-2]*muint1[it1+nt2+it3+nt4-1]*muint1[it1+nt2+it3-1]*muint1[it1+nt2]*mu2_mu1;
			      //mu6_mu5_mu4_mu3_mu2_mu1 = mu_01_2 * mu_01_2 * mu_01_2;
			      /* peak 1 */
			      R1_re[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw1[it5][it3][it1]);
			      R1_im[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw1[it5][it3][it1]);
			      R2_re[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw2[it5][it3][it1]);
			      R2_im[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw2[it5][it3][it1]);
			      R3_re[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw3[it5][it3][it1]);
			      R3_im[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw3[it5][it3][it1]);
			      R4_re[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw4[it5][it3][it1]);
			      R4_im[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw4[it5][it3][it1]);
			      if (n_levels>=2){
				/* peak 2 */
				mu6_mu5_mu4_mu3_mu2_mu1 = muint2[it1+nt2+it3+nt4+it5-2]*muint2[it1+nt2+it3+nt4-1]*muint1[it1+nt2+it3-1]*muint1[it1+nt2]*mu2_mu1;
				//mu6_mu5_mu4_mu3_mu2_mu1 = mu_12_2 * mu_01_2 * mu_01_2;
				R5_re[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw5[it5][it3][it1]);
				R5_im[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw5[it5][it3][it1]);
				R6_re[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw6[it5][it3][it1]);
				R6_im[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw6[it5][it3][it1]);
				R7_re[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw7[it5][it3][it1]);
				R7_im[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw7[it5][it3][it1]);
				R8_re[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw8[it5][it3][it1]);
				R8_im[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw8[it5][it3][it1]);
				if (n_levels>=3){
				  mu6_mu5_mu4_mu3_mu2_mu1 = muint3[it1+nt2+it3+nt4+it5-2]*muint3[it1+nt2+it3+nt4-1]*muint2[it1+nt2+it3-1]*muint2[it1+nt2]*mu2_mu1;
				  //mu6_mu5_mu4_mu3_mu2_mu1 = mu_23_2 * mu_12_2 * mu_01_2;
				  R9_re[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw9[it5][it3][it1]);
				  R9_im[it5][it3][it1]  += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw9[it5][it3][it1]);
				  R10_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw10[it5][it3][it1]);
				  R10_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw10[it5][it3][it1]);
				  R11_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw11[it5][it3][it1]);
				  R11_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw11[it5][it3][it1]);
				  R12_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw12[it5][it3][it1]);
				  R12_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw12[it5][it3][it1]);
				} /* end n_levels >= 3 */
				/* peak 4 */
				mu6_mu5_mu4_mu3_mu2_mu1 = muint2[it1+nt2+it3+nt4+it5-2]*muint2[it1+nt2+it3+nt4-1]*muint2[it1+nt2+it3-1]*muint2[it1+nt2]*mu2_mu1;
				//mu6_mu5_mu4_mu3_mu2_mu1 = mu_12_2 * mu_12_2 * mu_01_2;
				R13_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw13[it5][it3][it1]);
				R13_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw13[it5][it3][it1]);
				R14_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw14[it5][it3][it1]);
				R14_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw14[it5][it3][it1]);
				R15_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw15[it5][it3][it1]);
				R15_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw15[it5][it3][it1]);
				R16_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw16[it5][it3][it1]);
				R16_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw16[it5][it3][it1]);
				/* peak 5 */
				mu6_mu5_mu4_mu3_mu2_mu1 = muint1[it1+nt2+it3+nt4+it5-2]*muint1[it1+nt2+it3+nt4-1]*muint2[it1+nt2+it3-1]*muint2[it1+nt2]*mu2_mu1;
				//mu6_mu5_mu4_mu3_mu2_mu1 = mu_01_2 * mu_12_2 * mu_01_2;
				R17_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw17[it5][it3][it1]);
				R17_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw17[it5][it3][it1]);
				R18_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw18[it5][it3][it1]);
				R18_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw18[it5][it3][it1]);
				R19_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw19[it5][it3][it1]);
				R19_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw19[it5][it3][it1]);
				R20_re[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * cos(rdw20[it5][it3][it1]);
				R20_im[it5][it3][it1] += mu6_mu5_mu4_mu3_mu2_mu1 * sin(rdw20[it5][it3][it1]);
			      } /* end n_levels >= 2 */
			    } /* end order >= 5 and it5 */
			}//end order >=3 and it3
		}//end for it1
	    }//end if noncondon
	    
	if((isample+nsamples*(iproton-1))%20000==0)
            {
              /*write intermediate results*/
              normalizeResults(nt,dt,(isample+nsamples*(iproton-1)),
			       n_levels,flag_noncondon,mean_w,mean_mu,
                               S_re,S_im,
                               P1_re,P1_im,P2_re,P2_im,
			       P3_re,P3_im,P4_re,P4_im,
                               R1_re,R1_im,R2_re,R2_im,
			       R3_re,R3_im,R4_re,R4_im,
                               R5_re,R5_im,R6_re,R6_im,
			       R7_re,R7_im,R8_re,R8_im,
                               R9_re,R9_im,R10_re,R10_im,
			       R11_re,R11_im,R12_re,R12_im,
                               R13_re,R13_im,R14_re,R14_im,
			       R15_re,R15_im,R16_re,R16_im,
                               R17_re,R17_im,R18_re,R18_im,
			       R19_re,R19_im,R20_re,R20_im,
                               Sf_re,Sf_im,
                               P1f_re,P1f_im,P2f_re,P2f_im,
			       P3f_re,P3f_im,P4f_re,P4f_im,
			       P1tot_re,P1tot_im,P2tot_re,P2tot_im,
                               R1f_re,R1f_im,R2f_re,R2f_im,
			       R3f_re,R3f_im,R4f_re,R4f_im,
                               R5f_re,R5f_im,R6f_re,R6f_im,
			       R7f_re,R7f_im,R8f_re,R8f_im,
                               R9f_re,R9f_im,R10f_re,R10f_im,
			       R11f_re,R11f_im,R12f_re,R12f_im,
                               R13f_re,R13f_im,R14f_re,R14f_im,
			       R15f_re,R15f_im,R16f_re,R16f_im,
                               R17f_re,R17f_im,R18f_re,R18f_im,
			       R19f_re,R19f_im,R20f_re,R20f_im,
                               R1tot_re,R1tot_im,R2tot_re,R2tot_im,R3tot_re,R3tot_im,R4tot_re,R4tot_im);
	      writeResultsTime(base_name,nt,extension,n_levels,
			       Sf_re,Sf_im,
			       P1f_re,P1f_im,P2f_re,P2f_im,
			       P3f_re,P3f_im,P4f_re,P4f_im,
			       P1tot_re,P1tot_im,P2tot_re,P2tot_im,
                               R1f_re,R1f_im,R2f_re,R2f_im,
			       R3f_re,R3f_im,R4f_re,R4f_im,
                               R5f_re,R5f_im,R6f_re,R6f_im,
			       R7f_re,R7f_im,R8f_re,R8f_im,
                               R9f_re,R9f_im,R10f_re,R10f_im,
			       R11f_re,R11f_im,R12f_re,R12f_im,
                               R13f_re,R13f_im,R14f_re,R14f_im,
			       R15f_re,R15f_im,R16f_re,R16f_im,
                               R17f_re,R17f_im,R18f_re,R18f_im,
			       R19f_re,R19f_im,R20f_re,R20f_im,
                               R1tot_re,R1tot_im,R2tot_re,R2tot_im,R3tot_re,R3tot_im,R4tot_re,R4tot_im); 

            } //end if (isample+nsamples*(iproton-1)%2000==0 

	    
	  } //end isample loop
          
      //final results
      normalizeResults(nt,dt,nsamples*nprotons,
		       n_levels,flag_noncondon,mean_w,mean_mu,
		       S_re,S_im,
		       P1_re,P1_im,P2_re,P2_im,
		       P3_re,P3_im,P4_re,P4_im,
		       R1_re,R1_im,R2_re,R2_im,
		       R3_re,R3_im,R4_re,R4_im,
		       R5_re,R5_im,R6_re,R6_im,
		       R7_re,R7_im,R8_re,R8_im,
		       R9_re,R9_im,R10_re,R10_im,
		       R11_re,R11_im,R12_re,R12_im,
		       R13_re,R13_im,R14_re,R14_im,
		       R15_re,R15_im,R16_re,R16_im,
		       R17_re,R17_im,R18_re,R18_im,
		       R19_re,R19_im,R20_re,R20_im,
		       Sf_re,Sf_im,
		       P1f_re,P1f_im,P2f_re,P2f_im,
		       P3f_re,P3f_im,P4f_re,P4f_im,
		       P1tot_re,P1tot_im,P2tot_re,P2tot_im,
		       R1f_re,R1f_im,R2f_re,R2f_im,
		       R3f_re,R3f_im,R4f_re,R4f_im,
		       R5f_re,R5f_im,R6f_re,R6f_im,
		       R7f_re,R7f_im,R8f_re,R8f_im,
		       R9f_re,R9f_im,R10f_re,R10f_im,
		       R11f_re,R11f_im,R12f_re,R12f_im,
		       R13f_re,R13f_im,R14f_re,R14f_im,
		       R15f_re,R15f_im,R16f_re,R16f_im,
		       R17f_re,R17f_im,R18f_re,R18f_im,
		       R19f_re,R19f_im,R20f_re,R20f_im,
		       R1tot_re,R1tot_im,R2tot_re,R2tot_im,R3tot_re,R3tot_im,R4tot_re,R4tot_im);
      writeResultsTime(base_name,nt,extension,n_levels,
		       Sf_re,Sf_im,
		       P1f_re,P1f_im,P2f_re,P2f_im,
		       P3f_re,P3f_im,P4f_re,P4f_im,
		       P1tot_re,P1tot_im,P2tot_re,P2tot_im,
		       R1f_re,R1f_im,R2f_re,R2f_im,
		       R3f_re,R3f_im,R4f_re,R4f_im,
		       R5f_re,R5f_im,R6f_re,R6f_im,
		       R7f_re,R7f_im,R8f_re,R8f_im,
		       R9f_re,R9f_im,R10f_re,R10f_im,
		       R11f_re,R11f_im,R12f_re,R12f_im,
		       R13f_re,R13f_im,R14f_re,R14f_im,
		       R15f_re,R15f_im,R16f_re,R16f_im,
		       R17f_re,R17f_im,R18f_re,R18f_im,
		       R19f_re,R19f_im,R20f_re,R20f_im,
		       R1tot_re,R1tot_im,R2tot_re,R2tot_im,R3tot_re,R3tot_im,R4tot_re,R4tot_im); 
      
          //free dwint1 for the next t2,t4
          free_vector(dwint1,1,ntraject);
          free_vector(dwint2,1,ntraject);
          free_vector(dwint3,1,ntraject);
          free_vector(muint1,1,ntraject);
          free_vector(muint2,1,ntraject);
          free_vector(muint3,1,ntraject);
	  free(extension);

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
    free_f3tensor(R5_re,1,nt,1,nt,1,nt);
    free_f3tensor(R5_im,1,nt,1,nt,1,nt);
    free_f3tensor(R6_re,1,nt,1,nt,1,nt);
    free_f3tensor(R6_im,1,nt,1,nt,1,nt);
    free_f3tensor(R7_re,1,nt,1,nt,1,nt);
    free_f3tensor(R7_im,1,nt,1,nt,1,nt);
    free_f3tensor(R8_re,1,nt,1,nt,1,nt);
    free_f3tensor(R8_im,1,nt,1,nt,1,nt);
    free_f3tensor(R9_re,1,nt,1,nt,1,nt);
    free_f3tensor(R9_im,1,nt,1,nt,1,nt);
    free_f3tensor(R10_re,1,nt,1,nt,1,nt);
    free_f3tensor(R10_im,1,nt,1,nt,1,nt);
    free_f3tensor(R11_re,1,nt,1,nt,1,nt);
    free_f3tensor(R11_im,1,nt,1,nt,1,nt);
    free_f3tensor(R12_re,1,nt,1,nt,1,nt);
    free_f3tensor(R12_im,1,nt,1,nt,1,nt);
    free_f3tensor(R13_re,1,nt,1,nt,1,nt);
    free_f3tensor(R13_im,1,nt,1,nt,1,nt);
    free_f3tensor(R14_re,1,nt,1,nt,1,nt);
    free_f3tensor(R14_im,1,nt,1,nt,1,nt);
    free_f3tensor(R15_re,1,nt,1,nt,1,nt);
    free_f3tensor(R15_im,1,nt,1,nt,1,nt);
    free_f3tensor(R16_re,1,nt,1,nt,1,nt);
    free_f3tensor(R16_im,1,nt,1,nt,1,nt);
    free_f3tensor(R17_re,1,nt,1,nt,1,nt);
    free_f3tensor(R17_im,1,nt,1,nt,1,nt);
    free_f3tensor(R18_re,1,nt,1,nt,1,nt);
    free_f3tensor(R18_im,1,nt,1,nt,1,nt);
    free_f3tensor(R19_re,1,nt,1,nt,1,nt);
    free_f3tensor(R19_im,1,nt,1,nt,1,nt);
    free_f3tensor(R20_re,1,nt,1,nt,1,nt);
    free_f3tensor(R20_im,1,nt,1,nt,1,nt);
    printf("...response functions P1\n");
    free_matrix(P1_re,1,nt,1,nt);
    free_matrix(P1_im,1,nt,1,nt);
    printf("...response functions P2\n");
    free_matrix(P2_re,1,nt,1,nt);
    free_matrix(P2_im,1,nt,1,nt);
    printf("...response functions P3\n");
    free_matrix(P3_re,1,nt,1,nt);
    free_matrix(P3_im,1,nt,1,nt);
    printf("...response functions P4\n");
    free_matrix(P4_re,1,nt,1,nt);
    free_matrix(P4_im,1,nt,1,nt);
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
    free_f3tensor(R5f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R5f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R6f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R6f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R7f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R7f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R8f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R8f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R9f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R9f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R10f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R10f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R11f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R11f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R12f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R12f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R13f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R13f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R14f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R14f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R15f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R15f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R16f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R16f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R17f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R17f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R18f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R18f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R19f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R19f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R20f_re,1,nt,1,nt,1,nt);
    free_f3tensor(R20f_im,1,nt,1,nt,1,nt);
    free_f3tensor(R1tot_re,1,nt,1,nt,1,nt);
    free_f3tensor(R1tot_im,1,nt,1,nt,1,nt);
    free_f3tensor(R2tot_re,1,nt,1,nt,1,nt);
    free_f3tensor(R2tot_im,1,nt,1,nt,1,nt);
    free_f3tensor(R3tot_re,1,nt,1,nt,1,nt);
    free_f3tensor(R3tot_im,1,nt,1,nt,1,nt);
    free_f3tensor(R4tot_re,1,nt,1,nt,1,nt);
    free_f3tensor(R4tot_im,1,nt,1,nt,1,nt);
    free_matrix(P1f_re,1,nt,1,nt);
    free_matrix(P1f_im,1,nt,1,nt);
    free_matrix(P2f_re,1,nt,1,nt);
    free_matrix(P2f_im,1,nt,1,nt);
    free_matrix(P3f_re,1,nt,1,nt);
    free_matrix(P3f_im,1,nt,1,nt);
    free_matrix(P4f_re,1,nt,1,nt);
    free_matrix(P4f_im,1,nt,1,nt);
    free_matrix(P1tot_re,1,nt,1,nt);
    free_matrix(P1tot_im,1,nt,1,nt);
    free_vector(Sf_re,1,nt);
    free_vector(Sf_im,1,nt);
    
    printf("...response intermediates\n");
    free_vector(sdw,1,nt);
    free_matrix(    pdw1,1,nt,1,nt);
    free_matrix(    pdw2,1,nt,1,nt);
    free_f3tensor(    rdw1, 1,nt,1,nt,1,nt);
    free_f3tensor(    rdw2, 1,nt,1,nt,1,nt);
    free_f3tensor(    rdw3, 1,nt,1,nt,1,nt);
    free_f3tensor(    rdw4, 1,nt,1,nt,1,nt);
    free_f3tensor(    rdw5, 1,nt,1,nt,1,nt);
    free_f3tensor(    rdw6, 1,nt,1,nt,1,nt);
    free_f3tensor(    rdw7, 1,nt,1,nt,1,nt);
    free_f3tensor(    rdw8, 1,nt,1,nt,1,nt);
    free_f3tensor(    rdw9, 1,nt,1,nt,1,nt);
    free_f3tensor(    rdw10,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw11,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw12,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw13,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw14,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw15,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw16,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw17,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw18,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw19,1,nt,1,nt,1,nt);
    free_f3tensor(    rdw20,1,nt,1,nt,1,nt);
    
    //    printf("...dw_matrix\n");
    //    free_matrix(dw_matrix,1,nprotons,1,nsteps);
    

    }


int main(int argc, char *argv[]) {
  int i;
  time_t my_time=time(0); // time process started
  time_t end_time; // time process ended
  
  char *base_name,*parameter_file_name,*time_file_name,*string;
  float **t2_t4_pairs;
  int n_t2_t4_pairs;
  float ***w_matrices;
  float ***x_matrices;
  float *mean_x,*mean_w;
  float elapsed;

  int opt = 0, longIndex=0; //for getopt_long
  int nprotons;
  int nsteps;
  int nprotons_in_file;
  int nsteps_in_file;
  int proton_offset = 0; //TEST TEST
  int step_offset = 0; //TEST TEST

  printf("This is %s called as\n",argv[0]);
  for (i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  printf("On: %s\n",ctime(&my_time));
  printf("\n"); 

  /* 
   * initialize global arguments from the parameter file
   *
   */
  //initialize parameter file name
  if (asprintf(&parameter_file_name,"") < 0) nrerror("Failed to write string.");

  //find the parameter file from the input options
  /*
  while((opt = getopt_long( argc, argv, optString,longOpts, &longIndex))!=-1)
    switch(opt){
    case 'p': //parameter file
      free(parameter_file_name);
      if (asprintf(&parameter_file_name,"%s",optarg) < 0) nrerror("failed to write string");
      break;
    }
  */
  int count = 0;
  while(count<argc){
    count++;
    if (fnmatch(argv[count],"-p",FNM_CASEFOLD) == 0 || fnmatch(argv[count],"--param*",FNM_CASEFOLD) == 0){
      free(parameter_file_name);
      if (asprintf(&parameter_file_name,"%s",argv[count+1]) < 0) nrerror("Failed to write string");
      break;
    }
  }

  // error if parameter_file_name is still empty
  if (fnmatch(parameter_file_name,"",FNM_CASEFOLD) == 0) nrerror("please supply a parameter file name with the -p flag!");
  initialize_globalArgs(parameter_file_name);


  /* initialize base name */
  if (asprintf(&base_name,"",FNM_CASEFOLD) < 0) nrerror("Failed to write string.");

  //initialize time_file_name to empty
  if (asprintf(&time_file_name,"",FNM_CASEFOLD) < 0) nrerror("Failed to write string.");

  // initialize number of protons and steps
  nprotons = globalArgs.nprotons_in_file;
  nsteps = globalArgs.nsteps_in_file;

  /* read the command line parameters and set values as appropriate */
  while ((opt = getopt_long(argc, argv, optString, longOpts, &longIndex))!=-1){
    switch(opt){
    case 'x': /* set nprotons_in_file */
      if(sscanf(optarg,"%i",&nsteps_in_file)!=1) nrerror("failed reading commandline option -e");
      break;
    case 'f': /* set nprotons_in_file */
      if(sscanf(optarg,"%i",&nprotons_in_file)!=1) nrerror("failed reading commandline option -f");
      break;
    case 'g': /* set nt number of grid points */
      if(sscanf(optarg,"%i",&globalArgs.nt)!=1) nrerror("failed reading commandline option");
      break;
    case 'i': /* set ntint number of (pre)integration steps */
      if(sscanf(optarg,"%i",&globalArgs.ntint)!=1) nrerror("failed reading commandline option");
      break;
    case 'n': /* set nsteps number of steps to use in the md simulation */
      if(sscanf(optarg,"%i",&nsteps)!=1) nrerror("failed reading commandline option");
      break;
    case 'o': /* output base name */
      free(base_name);
      if (asprintf(&base_name,optarg,FNM_CASEFOLD) < 0) nrerror("Failed to write string.");
      break;
    case 'p': /* parameter file name */
      free(parameter_file_name);
      if (asprintf(&parameter_file_name,"%s",optarg) < 0) nrerror("failed to write string");
      break;
    case 'q': /* set nproton number of protons in trajectory */
      if(sscanf(optarg,"%i",&nprotons)!=1) nrerror("failed reading commandline option");
      break;
    case 'r': /* set proton offset in trajectory */
      if(sscanf(optarg,"%i",&proton_offset)!=1) nrerror("failed reading commandline option");
      break;
    case 's': /* set nskip number of steps to skip to make independent samples */
      if(sscanf(optarg,"%i",&globalArgs.nskip)!=1) nrerror("failed reading commandline option");
      break;
    case 't': /* set dt size of base time-step from md simulation */
      if (asprintf(&time_file_name,"%s",optarg) < 0 ) nrerror("failed to write string");
      if (DEBUG_LEVEL>=1) printf("time file %s\n",time_file_name);
      break;
    case 'v': /* increase verbosity */
      DEBUG_LEVEL+=1;
      break;
    case 'z': /* set flag for the FFT */
      if(sscanf(optarg,"%i",&globalArgs.flag_fft)!=1) nrerror("failed reading commandline option -z, make sure to specify a value such as -z 0 (false, no fft) or -z 1 (true, fft) ");
      break;
    case 'h': /* fall-through is intentional */
    case '?':
      display_usage(argc,argv);
      break;
    case 0: /* all options that don't have a short one character version */
      if (fnmatch("compress_output",longOpts[longIndex].name, FNM_CASEFOLD ) == 0 ){
	globalArgs.flag_compressoutput = 1;
      }
      if (fnmatch("dt",longOpts[longIndex].name, FNM_CASEFOLD ) == 0 ){
	sscanf(optarg,"%f",&globalArgs.dt);
      }      
      break;
    default:
      /* should never get here */
      fprintf(stderr,"WARNING!!! Commandline option %s not recognized!!! \n",longOpts[longIndex].name);
      exit(EXIT_FAILURE);
      break;
    }
  }

  // error if base_name is empty
  if (fnmatch(base_name,"",FNM_CASEFOLD) == 0) nrerror("please supply a base name with the -o flag!");

  //make sure options are okay
  display_options();

  // read time file 
  read_time_file(time_file_name,&t2_t4_pairs,&n_t2_t4_pairs);

  // read steps and protons if necessary
  if (globalArgs.nsteps_in_file == -1 || globalArgs.nprotons_in_file == -1)
    read_nsteps_and_nprotons_from_freq(parameter_file_name);

  // initialize w and x matrices
  w_matrices = malloc(globalArgs.n_levels*sizeof(float**));
  x_matrices = malloc(globalArgs.n_levels*sizeof(float**));
  mean_w = vector(0,globalArgs.n_levels-1);
  mean_x = vector(0,globalArgs.n_levels-1);

  /* read frequency file */
  read_w_files(nprotons,proton_offset,nsteps,step_offset,w_matrices,mean_w);
  read_x_files(nprotons,proton_offset,nsteps,step_offset,x_matrices,mean_x);

  /* output the mean freq shift to the param file */
  for (i=0;i<globalArgs.n_levels;i++){
    if (asprintf(&string,"mean_w_%i",i) < 0) nrerror("failed to write string");
    gaWriteFloat(parameter_file_name,string,mean_w[i]);
    free(string);
  }

  // main calculation
  freqTrajToR5(base_name,t2_t4_pairs,n_t2_t4_pairs,nsteps,nprotons,w_matrices,x_matrices,mean_w,mean_x);

  // clean up
  free_matrix(t2_t4_pairs,1,n_t2_t4_pairs,1,n_t2_t4_pairs);
  for (i=0;i<globalArgs.n_levels;i++)
    free_matrix(w_matrices[i],1,nprotons,1,nsteps);
  //if (globalArgs.flag_noncondon==1)
  for (i=0;i<globalArgs.n_levels;i++)
    free_matrix(x_matrices[i],1,nprotons,1,nsteps);
  free(w_matrices);
  free(x_matrices);
  free_vector(mean_w,0,globalArgs.n_levels);
  free_vector(mean_x,0,globalArgs.n_levels);

  // we're done
  end_time = time(0);
  printf("done: %s\n",ctime(&end_time));
  elapsed = difftime(end_time , my_time);
  printf("elapsed: %f s\n",elapsed);
  printf("\n"); 

  exit(EXIT_SUCCESS);
    
}//end main()
