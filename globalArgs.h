
extern int DEBUG_LEVEL;

struct globalArgs_t {
  //the parameters that need to be globally visible can be stord in
  //this data structure
  int test; 
  float dt;
  int nt;
  int nmols_in_file;
  int nprotons_in_file;
  int nsteps_in_file;
  int ntint;
  int nskip;
  int order;
  int fit_order; //order of expansion
  int n_levels; //max number of excited states
  float **a; //zero based matrix size (0,(order-1)/2,0,fit_order)
  float **b; //zero based matrix size (0,(order-1)/2,0,fit_order)
  float *mu_mug; 
  float q_H; //charge on the proton (to convert force to field)
  int flag_massweightedforces;
  int flag_noncondon;
  int flag_twolevelsystem;
  int flag_compressoutput;
  int flag_compressedinput;
  int flag_fft;
  char **w_file_names;
  char **mu_file_names;
} globalArgs;

void display_options( void );

//void read_input_parameters( char *parameter_file_name );

void initialize_globalArgs( const char *parameter_file_name );

//int gaRewriteString(const char *parameter_file_name,const char *name,const char *val);
//int gaAppendString(const char *parameter_file_name,const char *name,const char *val);

int gaRemoveParameter(const char *parameter_file_name,const char *name);
int gaWriteString(const char *parameter_file_name,const char *name,const char *val);
int gaWriteInt(const char *parameter_file_name,const char *name,const int val);
int gaWriteFloat(const char *parameter_file_name,const char *name,const float val);
