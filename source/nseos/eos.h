#ifndef H_EOS
#define H_EOS

struct transition_qtt
{
    double nt;
    double pt;
};

int calc_equation_of_state(struct parameters satdata, double p, 
        struct transition_qtt *tqtt, char *outfile[]);

void eval_transition_qtt(struct parameters satdata, double p,
        struct transition_qtt *tqtt);

#endif // H_EOS
