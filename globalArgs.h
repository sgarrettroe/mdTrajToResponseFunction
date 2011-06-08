
static int DEBUG_LEVEL = 0;

struct globalArgs_t {
  //the parameters that need to be globally visible can be stord in
  //this data structure
  int test; 
  float dt;
  int nt;
  int ntint;
  int nskip;
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
  int flag_compressedinput;
} globalArgs;

void display_options( void );

void read_input_parameters(char *parameter_file_name);

int gaRewriteString(const char *parameter_file_name,const char *name,const char *val);
int gaRemoveString(const char *parameter_file_name,const char *name);
int gaAppendString(const char *parameter_file_name,const char *name,const char *val);
int gaWriteString(const char *parameter_file_name,const char *name,const char *val);
