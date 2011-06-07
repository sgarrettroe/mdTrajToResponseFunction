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
  printf("flag_compressedinput \t%d\n",globalArgs.flag_compressedinput); 

  printf("\n\n");
}

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
    //if (sscanf(line,"%as = %f",&name,&val) < 2){ //only works in gnu's libcnot freebsd!!!
    if (sscanf(line,"%s = %f",name,&val) < 2){
      if (DEBUG_LEVEL>=1) printf("skip input parameter line: \n%s",line);
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
    if (fnmatch("flag_compressedinput",name,FNM_CASEFOLD)==0){
      globalArgs.flag_compressedinput= (int) val;
    }
    
  } //end while(feof(fid)==0)

  fclose(fid);

}
