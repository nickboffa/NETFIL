#ifndef STRATEGY_H
#define STRATEGY_H
#include "network.h"
#include "mda.h"

enum class AreaType { AS_ALL, VILLAGE, SUBVILLAGE };

struct StrategyInputs {
    int tas_every_k_years = 2;   // k
    AreaType area_type = AreaType::AS_ALL;
    double start_threshold_p = 1.0;  // p%
    int tas_sample_size = 1500;      // adults 15+, antigen TAS
    int max_years = 21;              // 2010..2030 inclusive
};

struct StrategyOutputs {
    int years_to_elimination = -1;          // first ABS year that met rule (record-only)
    long long total_infected_people = 0;    // sum over years of I (year-end snapshot)
    long long total_years_with_infection = 0; // PERSON-YEARS infected (sum of I each year)
    long long total_mda_doses = 0;          // total doses given (sum over years)
    int mda_years = 0;                      // years where any MDA actually occurred
    int tas_surveys = 0;                    // TAS cadence years conducted post-threshold
    int mf_zero_at_end = 0;                 // 1 if final year’s regional MF == 0, else 0
};

void run_elimination_strategy(region& rgn,
                              mda_strat base_mda,
                              const StrategyInputs& in,
                              StrategyOutputs& out);

#endif
