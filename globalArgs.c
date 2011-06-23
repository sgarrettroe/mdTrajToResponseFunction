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
#include <sys/param.h>
#include "globalArgs.h"
#include <limits.h>
#include "mymath.h"

int DEBUG_LEVEL = 0;

void display_options( void ){
  // show the user what the current global arguments (options) are
  int i,j;

  printf("\n");
  printf("Global options (globalArgs)\n");
  printf("\n");

  printf("test =           \t%12d\n",globalArgs.test); 
  printf("dt =             \t%12f ps\n",globalArgs.dt); 
  printf("nt =             \t%12i\n",globalArgs.nt); 
  printf("nsteps_in_file = \t%12i\n",globalArgs.nsteps_in_file); 
  printf("nprotons_in_file =\t%12i\n",globalArgs.nprotons_in_file); 
  printf("nmols_in_file =\t%12i\n",globalArgs.nmols_in_file); 
  printf("ntint =          \t%12i\n",globalArgs.ntint); 
  printf("nskip =          \t%12i\n",globalArgs.nskip); 
  printf("order =          \t%12i\n",globalArgs.order); 
  printf("n_levels =       \t%12i\n",globalArgs.n_levels); 
  printf("fit_order =      \t%12i\n",globalArgs.fit_order); 
  printf("q_H =            \t%12i\n",globalArgs.q_H); 
  for (i=0;i<=(globalArgs.order-1)/2;i++){
    for (j=0;j<=globalArgs.fit_order;j++){
      printf("a_%i_%i =              \t%12f\n",i,j,globalArgs.a[i][j]);
    }
  }
  for (i=0;i<=(globalArgs.order-1)/2;i++){
    for (j=0;j<=globalArgs.fit_order;j++){
      printf("b_%i_%i =              \t%12f\n",i,j,globalArgs.b[i][j]);
    }
  }
  for (j=0;j<=globalArgs.fit_order;j++){
    printf("mu_mug_%i =      \t%12f\n",j,globalArgs.mu_mug[j]); 
  }
  for (i=0;i<globalArgs.n_levels;i++){
    printf("w_file_%i =              \t%s \n",i,globalArgs.w_file_names[i]); 
  }
  for (i=0;i<globalArgs.n_levels;i++){
    printf("mu_file_%i =             \t%s \n",i,globalArgs.mu_file_names[i]); 
  }
    
  printf("flag_massweightedforces \t%d\n",globalArgs.flag_massweightedforces); 
  printf("flag_noncondon          \t%d\n",globalArgs.flag_noncondon); 
  printf("flag_twolevelsystem     \t%d\n",globalArgs.flag_twolevelsystem); 
  printf("flag_compressoutput     \t%d\n",globalArgs.flag_compressoutput); 
  printf("flag_compressedinput    \t%d\n",globalArgs.flag_compressedinput); 

  printf("\n\n");
}

void read_input_parameters(const char *parameter_file_name){
  //read the input paramter file
  FILE *fid;
  char *line;
  char name[100];
  //char *name;
  float val;
  char sval[100];
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



#ifdef BSD /* defined in sys/param.h */
#warning "BSD fgetln()"
    char *unparsed_line;
    unparsed_line = fgetln(fid,&len);
    if (DEBUG_LEVEL>=3) printf("len  %i\n",(int)len);
    line = calloc(len,sizeof(char));
    strncpy(line,unparsed_line,(int)len);
#else
#warning "GLIBC getline()"
      int bytes_read;
     //getline is not in FreeBSD's libc yet ... 
      len = 100;
    line = (char *) malloc(len * sizeof(char));
    bytes_read = getline(&line, &len, fid);
    /* //dont ask me why this isn't working correctly and always gives -1 even though it reads the line...???
    if (bytes_read == -1){
      fprintf(stderr,"Unable to read line of parameter file %s\n",parameter_file_name);
      exit(EXIT_FAILURE);
    }
    */
#endif /* BSD */
    if (DEBUG_LEVEL>=2)
      printf("line %i: %s",count,line);
    
    if (len == 0){
      continue;
    } //end if len>0

    if(strncmp(line,"#",1)==0) {
      // if the first character is a comment character do nothing on this loop
      if (DEBUG_LEVEL>=1) printf("comment: %s",line);
      continue;
    }

    // if it is not a comment character try to process it
    //if (sscanf(line,"%as = %f",&name,&val) < 2){ //only works in gnu's libc not freebsd!!!
    if (sscanf(line,"%s = %s",name,sval) < 2){
      if (DEBUG_LEVEL>=1) printf("skip input parameter line: \n%s",line);
      continue;
    }
    sscanf(sval,"%f",&val);

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
      globalArgs.n_levels = (val+1)/2;
      globalArgs.w_file_names = malloc(globalArgs.n_levels*sizeof(char*));
      globalArgs.mu_file_names = malloc(globalArgs.n_levels*sizeof(char*));
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
      if (globalArgs.flag_twolevelsystem == 1){
	globalArgs.n_levels = 1;
	//realloc w_file_names and mu_file_names?
	realloc(globalArgs.w_file_names,  globalArgs.n_levels*sizeof(char*));
	realloc(globalArgs.mu_file_names, globalArgs.n_levels*sizeof(char*));
      }
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
    if (fnmatch("nt",name,FNM_CASEFOLD)==0){
      globalArgs.nt = (int) val;
    }
    if (fnmatch("ntint",name,FNM_CASEFOLD)==0){
      globalArgs.ntint = (int) val;
    }
    if (fnmatch("nskip",name,FNM_CASEFOLD)==0){
      globalArgs.nskip = (int) val;
    }
    if (fnmatch("w_file_*",name,FNM_CASEFOLD)==0){
      sscanf(name,"w_file_%d",&i);
      if ( DEBUG_LEVEL >= 1 ) printf("assigning w_file_%i = %s\n",i,sval);
      if ( i < globalArgs.n_levels ) 
	if (asprintf(&globalArgs.w_file_names[i],"%s",sval) < 0){
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
    }
    if (fnmatch("mu_file_*",name,FNM_CASEFOLD)==0){
      sscanf(name,"mu_file_%d",&i);
      if ( DEBUG_LEVEL >= 1 ) printf("assigning mu_file_%i = %s\n",i,sval);
      if ( i < globalArgs.n_levels ) 
	if (asprintf(&globalArgs.mu_file_names[i],"%s",sval) < 0){
	  fprintf(stderr,"failed to write string");
	  exit(EXIT_FAILURE);
	}
    }
    if (fnmatch("nmols_in_file",name,FNM_CASEFOLD)==0){
      globalArgs.nmols_in_file = (int) val;
    }
    if (fnmatch("nprotons_in_file",name,FNM_CASEFOLD)==0){
      globalArgs.nprotons_in_file = (int) val;
    }
    if (fnmatch("nsteps_in_file",name,FNM_CASEFOLD)==0){
      globalArgs.nsteps_in_file = (int) val;
    }

    
  } //end while(feof(fid)==0)

  fclose(fid);

}

void initialize_globalArgs(const char *parameter_file_name ){ 
  // generate the default values of the global arguments
  globalArgs.test = 1;
  globalArgs.flag_compressoutput = 0;
  globalArgs.flag_compressedinput = 0;
  globalArgs.nsteps_in_file = -1;
  globalArgs.nprotons_in_file = -1;
  globalArgs.nmols_in_file = -1;
  globalArgs.nskip = 15;
  globalArgs.flag_fft = 1;

  read_input_parameters(parameter_file_name);
}


/*
 * general functions
 */
int gaRemoveParameter(const char *parameter_file_name,const char *name){
  char *command;
  int ret;

  //  if (asprintf(&command,"perl -i.bak -pe \'s/^%s\\s*=.*//g\' %s",name,parameter_file_name) < 0){
  if (asprintf(&command,"perl -i.bak -ne \'print unless /^%s[\\s=]/\' %s",name,parameter_file_name) < 0){
    fprintf(stderr,"failed to write string");
    exit(EXIT_FAILURE);
  }
  if (DEBUG_LEVEL>=1)  printf("%s\n",command);
  ret = system(command);
      if(ret!=0)
	fprintf(stderr,"Error %d executing (system()) %s",ret,command);

  return(0);
}

/*
 * string functions
 */
int gaRewriteString(const char *parameter_file_name,const char *name,const char *val){
  /* This is still in development so test it carefully */
  char *command;
  int ret;

  if (asprintf(&command,"perl -i .bak -pe \'s/^%s\\s*=\\s*(\\w*.\\w*)\\s*#(.*) / %s = %s #?$2/g\' %s",name,name,val,parameter_file_name) < 0){
    fprintf(stderr,"failed to write string");
    exit(EXIT_FAILURE);
  }
  if (DEBUG_LEVEL>=1) printf("%s\n",command);
  ret = system(command);
      if(ret!=0)
	fprintf(stderr,"Error %d executing (system()) %s",ret,command);

  
  return(0);
}
int gaAppendString(const char *parameter_file_name,const char *name, const char *val){
  char *command;
   int ret;

  ret = 0;
  if (asprintf(&command,"echo %s = %s >> %s",name,val,parameter_file_name) < 0){
    fprintf(stderr,"failed to write string");
    exit(EXIT_FAILURE);
  }
  if (DEBUG_LEVEL>=1) printf("%s\n",command);
  ret = system(command);
  if(ret!=0)
    {
      fprintf(stderr,"Error %d executing (system()) %s",ret,command);
      exit(EXIT_FAILURE);
    }

  return(ret);
}
int gaWriteString(const char *parameter_file_name,const char *name,const char *val){
  int ret=0;
  
  ret = gaRemoveParameter(parameter_file_name,name);
  ret = gaAppendString(parameter_file_name,name,val);
  return(ret);
}

/*
 * int functions
 */
int gaAppendInt(const char *parameter_file_name,const char *name, const int val){
  char *command;
   int ret;

  ret = 0;
  if (asprintf(&command,"echo %s = %i >> %s",name,val,parameter_file_name) < 0){
    fprintf(stderr,"failed to write string");
    exit(EXIT_FAILURE);
  }
  if (DEBUG_LEVEL>=1) printf("%s\n",command);
  ret = system(command);
  if(ret!=0)
    {
      fprintf(stderr,"Error %d executing (system()) %s",ret,command);
      exit(EXIT_FAILURE);
    }

  return(ret);
}
int gaWriteInt(const char *parameter_file_name,const char *name,int val){
  int ret = 0;
  ret = gaRemoveParameter(parameter_file_name,name);
  ret = gaAppendInt(parameter_file_name,name,val);
  return(ret);
}

/*
 * float functions
 */
int gaAppendFloat(const char *parameter_file_name,const char *name, const float val){
  char *command;
   int ret;

  ret = 0;
  if (asprintf(&command,"echo %s = %f >> %s",name,val,parameter_file_name) < 0){
    fprintf(stderr,"failed to write string");
    exit(EXIT_FAILURE);
  }
  if (DEBUG_LEVEL>=1) printf("%s\n",command);
  ret = system(command);
  if(ret!=0)
    {
      fprintf(stderr,"Error %d executing (system()) %s",ret,command);
      exit(EXIT_FAILURE);
    }

  return(ret);
}
int gaWriteFloat(const char *parameter_file_name,const char * name,float val){
  int ret = 0;
  ret = gaRemoveParameter(parameter_file_name,name);
  ret = gaAppendFloat(parameter_file_name,name,val);
  return(ret);
}

