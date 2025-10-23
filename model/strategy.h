#ifndef STRATEGY_H
#define STRATEGY_H

#include "network.h"
#include "mda.h"

enum class AreaType { AS_ALL, VILLAGE, SUBVILLAGE };

struct StrategyInputs {
    int tas_every_Y_years = 2;        // TAS cadence Y (years)
    AreaType area_type = AreaType::AS_ALL;
    double start_threshold_p = 1.0;   // p_MDA threshold in percent (e.g., 1.0)
    int tas_sample_size = 1500;       // adults 15+, antigen TAS (still used for TAS rounds)
    int max_years = 31;               // safety; controller caps to END_YEAR
};

struct StrategyOutputs {
    // NEW: end-of-2010 prevalences (after sim(year=0))
    double mf_prev_2010_end = -1.0;              // % MF, all ages, region-wide
    double ant_prev_2010_end_allages = -1.0;     // % antigen, ALL ages, exact (not sampled)

    int years_to_elimination = -1;               // first ABS year that met rule (record-only)
    long long total_infected_people = 0;         // UNIQUE persons ever MF-positive during run
    long long total_years_with_infection = 0;    // PERSON-YEARS MF-positive (sum of I each year)
    long long total_mda_doses = 0;               // total doses given (sum over years)
    int mda_years = 0;                           // years where any MDA actually occurred
    int tas_surveys = 0;                         // TAS cadence years conducted post-threshold
    int mf_zero_at_end = 0;                      // 1 if final year’s regional MF == 0, else 0
    int final_year = -1;                         // ABS year when MF first hits 0, else END_YEAR
};

void run_elimination_strategy(region& rgn,
                              mda_strat base_mda,
                              const StrategyInputs& in,
                              StrategyOutputs& out);

#endif // STRATEGY_H
