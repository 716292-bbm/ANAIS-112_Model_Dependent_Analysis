void write_difrate (char* coefficients, char* path, char* halo_path, double rhochi, void* input_difcros, char* model);
void write_difrate_total(char* coefficients, char* path, char* halo_path, double rhochi, void* input_difcros, char* model);
void write_binned_data( char* coefficients, char* exp_path, char* path, char * halo_path, double rhochi, void * input_difcros, char* model);
void write_binned_data_each( char* coefficients, char* exp_path, char* path, char * halo_path, double rhochi, void * input_difcros, char* model);
const char* write_header_NP(char* header);
void write_binned_data_NP( char* coefficients, char* exp_path, char* path, char * halo_path, double rhochi, void * input_difcros, char* model);