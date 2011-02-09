#include "tempo2pred.h"

long double get_period(const struct arch_params ap, const T2Predictor *t2p,
        const int subint)
{
	double dmjd = (double)ap.imjd + ap.offs + (double)ap.smjd/86400.0 +
		ap.lst_subs[subint]/86400.0;

	long double period = 1.0/(T2Predictor_GetFrequency(t2p, dmjd, ap.freq0));
	return period;
}

long double get_phase(const struct arch_params ap, const T2Predictor *t2p,
        const int subint, const float offs_sub)
{
	double dmjd = (double)ap.imjd + ap.offs + (double)ap.smjd/86400.0 +
		ap.lst_subs[subint]/86400.0 + offs_sub/86400.0;

    long double phase = T2Predictor_GetPhase(t2p, dmjd, ap.freq0);

    return phase - floor(phase);
}

T2Predictor create_T2Predicator(char *pred_file)
{
    T2Predictor t2p;
    T2Predictor_Read(&t2p, pred_file);

    return t2p;
}

float compute_offs_sub(const long double reference_phase,
        const long double subint_phase, const long double subint_period,
        const float search_offs_sub)
{
    return search_offs_sub - (((float)subint_phase - (float)reference_phase) * (float)subint_period);
}
