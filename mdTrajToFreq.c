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
#include "globalArgs.h"
#include "mymath.h"

/*
 *  7.2.2007 S. Garrett-Roe to read an MD trajectory (one file of coordinates and 
 *  one file of forces) and calculate a frequnecy trajectory.
 *
 *  26 May 2011 start to rewrite with better input of options, the
 *  possibility to use other mappings of coords and forces to
 *  frequency, the possibility to output dipole and anharmonicity
 */

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


void mdTrajToFreq(const char* parameter_file_name, const char* coord_file_name, const char* force_file_name, const char* base_name)
{
  FILE *coord_fid,*force_fid;
  FILE *fid; 
  FILE **w_fid_array;// = malloc(3*sizeof(FILE*));//TEST TEST
  FILE **x_fid_array;// = malloc(3*sizeof(FILE*));//TEST TEST!!!
  char *string,*fname,*pname,*param_name;
  //,coord_file_name[100],force_file_name[100],base_name[100],
  char *command_template;
  
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

  const int flag_massweightedforces = globalArgs.flag_massweightedforces;
  const int flag_noncondon = globalArgs.flag_noncondon;
  const int flag_twolevelsystem = globalArgs.flag_twolevelsystem;
  float **a = globalArgs.a; //why can't these be const?
  float **b = globalArgs.b;
  const float *mu_mug = globalArgs.mu_mug;
  const int n_levels = (globalArgs.flag_twolevelsystem == 0 ? (globalArgs.order+1)/2 : 1);
  const int fit_order = globalArgs.fit_order;
  const int flag_compressoutput = globalArgs.flag_compressoutput;
  const int flag_compressedinput = globalArgs.flag_compressedinput;
  const float q_H = globalArgs.q_H;

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
  if (flag_compressedinput == 0){
    //if input is ascii
    if (asprintf(&command_template,"wc -l %%s") < 0){
      fprintf(stderr,"failed to write string");
      exit(EXIT_FAILURE);
    }
  } else {
    //if compressed input
    if (asprintf(&command_template,"gzip -cd %%s | wc -l") < 0){
      fprintf(stderr,"failed to write string");
      exit(EXIT_FAILURE);
    }
  }
  //  if (asprintf(&string,"wc -l %s",coord_file_name)<0){
  if (asprintf(&string,command_template,coord_file_name)<0){
    fprintf(stderr,"failed to write string");
    exit(EXIT_FAILURE);
  }
  free(command_template);
  fid = popen(string,"r");
  if(fid==NULL){
    fprintf(stderr,"Opening '%s' pipe failed.",string);
    exit(EXIT_FAILURE);
  }
  fscanf(fid,"%ld",&nsteps);
  pclose(fid);
  free(string);
  printf("Found %ld steps in file %s\n", nsteps, coord_file_name);
  
  if (flag_compressedinput == 0){
    //if input is ascii
    if (asprintf(&command_template,"head -n1 %%s | wc -w") < 0){
      fprintf(stderr,"failed to write string");
      exit(EXIT_FAILURE);
    }
  } else {
    //if compressed input
    if (asprintf(&command_template,"gzip -cd %%s | head -n1 | wc -w") < 0){
      fprintf(stderr,"failed to write string");
      exit(EXIT_FAILURE);
    }
  }
  //  if (asprintf(&string,"head -n 1 %s | wc -w",coord_file_name) < 0){
  if (asprintf(&string,command_template,coord_file_name)<0){
    fprintf(stderr,"failed to write string");
    exit(EXIT_FAILURE);
  }
  free(command_template);
  fid = popen(string,"r");
  if(fid==NULL){
    fprintf(stderr,"Opening '%s' pipe failed.",string);
    exit(EXIT_FAILURE);
  }
  fscanf(fid,"%ld",&natoms);
  pclose(fid);
  free(string);
  
  natoms = (natoms-1)/3; //substract 1 for the time stamp 
  //and div by 3 for x y z coords
  nmols = natoms/3; //3 atoms per H2O molecule
  nprotons = nmols*2;//2 protons
  printf("Found %ld atoms, %ld molecules in file %s\n", natoms, nmols, coord_file_name);

  // add these values to the parameter file
  globalArgs.nprotons_in_file = nprotons;
  gaWriteInt(parameter_file_name,"nprotons_in_file",nmols); //only output one freq per molecule
  globalArgs.nsteps_in_file = nprotons;
  gaWriteInt(parameter_file_name,"nsteps_in_file",nsteps);

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
   */
  
  
  // open coordinate and force files
  printf("Open trajectory files\n");  
  if (flag_compressedinput==0){
    // plain text
    coord_fid = fopen(coord_file_name,"rt");
    if(coord_fid==NULL) nrerror("opening coordinate file failed.");
    force_fid = fopen(force_file_name,"rt");
    if(force_fid==NULL) nrerror("opening force file failed.");
  } else {
    // compressed
    if (asprintf(&string,"gzip -cd %s",coord_file_name) < 0)
      {
	fprintf(stderr,"failed to write string");
	exit(EXIT_FAILURE);
      }
    coord_fid = popen(string,"r");
    free(string);
    if(coord_fid==NULL) nrerror("opening coordinate file gzip pipe failed.");

    if (asprintf(&string,"gzip -cd %s",force_file_name) < 0)
      {
	fprintf(stderr,"failed to write string");
	exit(EXIT_FAILURE);
      }
    force_fid = popen(string,"r");
    free(string);
    if(force_fid==NULL) nrerror("opening force file gzip pipe failed.");    
  }

  //set these up as /dev/NULL if not desired and then just calculate everything
  w_fid_array = malloc(n_levels*sizeof(FILE*));
  x_fid_array = malloc(n_levels*sizeof(FILE*));

  //open output files
  if (flag_compressoutput==0){
    //open text files
    if (DEBUG_LEVEL>=0) printf("open files for text output\n");
    for (i_level=0;i_level<n_levels;i_level++){
      //freq files

      //build the name of the output file
      if (asprintf(&fname,"%s_dw%i.dat",base_name,(int)i_level) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}

      //build the name of the variable in the parameter file param_name
      if (asprintf(&param_name,"w_file_%i",i_level) < 0) {
	fprintf(stderr,"failed to write string");
	exit(EXIT_FAILURE);
      }
      
      //open the file fname
      w_fid_array[i_level]=fopen(fname,"wt"); 
      if (w_fid_array[i_level]==NULL){
	fprintf(stderr,"Error opening file %s for output\n",fname);
	exit(EXIT_FAILURE);
      }

      //save the file name to the parameter file 
      gaWriteString(parameter_file_name,param_name,fname);
      free(param_name);
      free(fname);

      //mu files
      //build the name of the variable in the parameter file
      if (asprintf(&param_name,"mu_file_%i",i_level) < 0) {
	fprintf(stderr,"failed to write string");
	exit(EXIT_FAILURE);
      }
      if (flag_noncondon==1){
	//build the file name
	if (asprintf(&fname,"%s_mu%i.dat",base_name,(int)i_level) < 0)
	  {
	    fprintf(stderr,"failed to write string");
	    exit(EXIT_FAILURE);
	  }
    
	//save the file name to the parameter file 
	gaWriteString(parameter_file_name,param_name,fname);
      } else {
	//use dev null if no file
	if (asprintf(&fname,"/dev/null") < 0)
	  {
	    fprintf(stderr,"failed to write string");
	    exit(EXIT_FAILURE);
	  }

	//clear the name if it is in the parameter file
	gaRemoveParameter(parameter_file_name,param_name);
      }      
      //open the file (which can be a real file or dev null)
      x_fid_array[i_level]=fopen(fname,"wt"); 
      if (x_fid_array[i_level]==NULL){
	fprintf(stderr,"Error opening file %s for output\n",fname);
	exit(EXIT_FAILURE);
      }
      free(param_name);
      free(fname);
    }
  } else {
    //open a pipe to gz to compress the output
    if (DEBUG_LEVEL>=0) printf("open pipe to gzip files for compressed output\n");
    for (i_level=0;i_level<n_levels;i_level++){
      //freq files
      
      //build the name of the variable in the parameter file param_name
      if (asprintf(&param_name,"w_file_%i",i_level) < 0) {
	fprintf(stderr,"failed to write string");
	exit(EXIT_FAILURE);
      }
      
      //build the name of the output file
      if (asprintf(&fname,"%s_dw%i.dat.gz",base_name,(int)i_level) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}

      //build the pipe
      if (asprintf(&pname,"gzip > %s",fname,(int)i_level) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
   
      if (DEBUG_LEVEL>=2) printf("%s\n",pname);
      w_fid_array[i_level] = popen(pname,"w"); //TEST TEST
      if (w_fid_array[i_level]==NULL){
	fprintf(stderr,"Error opening pipe %s for output\n",pname);
	exit(EXIT_FAILURE);
      }

      //save the file name to the parameter file 
      gaWriteString(parameter_file_name,param_name,fname);
      free(fname);
      free(pname);
      free(param_name);

      //mu files

      //build the name of the variable in the parameter file
      if (asprintf(&param_name,"mu_file_%i",i_level) < 0) {
	fprintf(stderr,"failed to write string");
	exit(EXIT_FAILURE);
      }

      if (flag_noncondon==1){
	//build file name
	if (asprintf(&fname,"%s_mu%i.dat.gz",base_name,(int)i_level) < 0)
	  {
	    fprintf(stderr,"failed to write string");
	    exit(EXIT_FAILURE);
	  }

	//save the file name to the parameter file 
	gaWriteString(parameter_file_name,param_name,fname);

      } else {

	// use dev null if not outputting to a real file
	if (asprintf(&fname,"/dev/null") < 0)
	  {
	    fprintf(stderr,"failed to write string");
	    exit(EXIT_FAILURE);
	  }

	//clear the name if it is in the parameter file
	gaRemoveParameter(parameter_file_name,param_name);

      }        

      //build pipe name
      if (asprintf(&pname,"gzip > %s",fname,(int)i_level) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}

      x_fid_array[i_level]=popen(pname,"w");
      if (x_fid_array[i_level]==NULL){
	fprintf(stderr,"Error opening pipe %s for output\n",pname);
	exit(EXIT_FAILURE);
      }
      free(fname);
      free(pname);
      free(param_name);

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
  if (DEBUG_LEVEL>=1) printf("start reading coord and force files\n");
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


int main( int argc, char *argv[] ) {
  //This is the main function, the glue that sticks the pieces
  //together

  //declare main() variables
  int opt = 0, longIndex=0; //for getopt_long
  char *parameter_file_name,*coord_file_name,*force_file_name,*base_name;

  /* initialize global arguments */
  initialize_globalArgs();

  /* initialize base name */
  if (asprintf(&base_name,"",FNM_CASEFOLD) < 0) nrerror("Failed to write string.");

  printf("process options\n");
  /* process command line arguments */
  while((opt = getopt_long( argc, argv, optString,longOpts, &longIndex))!=-1)
    //if (DEBUG_LEVEL>2) printf("%c",opt);
    switch(opt){
    case 'c': //coord file
      if (asprintf(&coord_file_name,"%s",optarg) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
      if (fnmatch("*.gz",coord_file_name,FNM_CASEFOLD)==0){
	globalArgs.flag_compressedinput = 1;
      }
      break;
    case 'f': //forces file
      if (asprintf(&force_file_name,"%s",optarg) < 0)
	{
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
      if (fnmatch("*.gz",force_file_name,FNM_CASEFOLD)==0){
	globalArgs.flag_compressedinput = 1;
      }
      break;
    case 'o': //output base name
      free(base_name);
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

  // error if base_name is empty
  if (fnmatch(base_name,"",FNM_CASEFOLD) == 0) nrerror("please supply a base name with the -o flag!");

  //input required data
  //read_input_data(parameter_file_name);

  //do the main calculation
  mdTrajToFreq(parameter_file_name,coord_file_name,force_file_name,base_name);

  //output the results

  //cleanup

  // we're done!
  printf("done\n");
  exit(EXIT_SUCCESS);

}//end main()

