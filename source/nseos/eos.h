#ifndef H_EOS
#define H_EOS

void calc_equation_of_state(struct parameters satdata, double p, char *outfile[]);

struct transition_qtt
{
    double nt;
    double pt;
};
void eval_transition_qtt(struct parameters satdata, double p,
        struct transition_qtt *tqtt);

#endif // H_EOS
