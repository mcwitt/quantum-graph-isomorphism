typedef struct
{
    double smin, smax;
    double emin, emax;
    int ns, nh;

    int itermax;
    double eps;

    int num_files;
    char **files;
} params_t;

void params_defaults(params_t *p);
void params_from_cmd(params_t *p, int argc, char *argv[]);