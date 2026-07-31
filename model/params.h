#ifndef headers_h
#define headers_h

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <random>
#include <ctime>
#include <set>
#include <map>

using namespace std;

constexpr double IMMATURE_PERIOD_MEAN     = 30*9;   //mean immature period
constexpr double IMMATURE_PERIOD_MEAN_STD = 30*0.1; //STD dev of immature period

constexpr double MATURE_PERIOD_MEAN       = 364*5;   //mean mature period
constexpr double MATURE_PERIOD_MEAN_STD   = 364*0.1; //STD dev of mature period

constexpr double PROPORTION_MALE_WORM  = 0.5; //proportion of worms that are male
constexpr double PROPORTION_MALE_AGENT = 0.5; //proportion of agents that are male

// Potential improvement: infer number of age groups from pop_age_dists.csv?
constexpr int N_AGE_GROUPS    = 16; //number of 5-year age brackets (for seeding pop)
constexpr int WIDTH_AGE_GROUPS = 5; // 0-4, 5-9, ... 75-79

// ABC_FITTING must remain a #define — it is used in a preprocessor #if directive
#define ABC_FITTING true

#if ABC_FITTING
constexpr int SIM_YEARS = 7;
#else
constexpr int SIM_YEARS = 21;
#endif


constexpr int START_YEAR = 2010; // Model starting year

constexpr double COMMUTING_PROP      = 0.5;          //proportion of group that commute daily (over 5 years old)
constexpr int    RECALC_YEARS        = 100;           //how often we want to recalc commuters
constexpr char   DISTANCE_TYPE       = 'r';           // r for road distance, e for euclidean
constexpr double DAILY_PROB_LOSE_ANT = 0.992327946;  //set so the half-life is 90 days i.e. pow(0.5,1/90)



double random_real();
double normal(double mean, double stddev);
int poisson(double rate);
double bite_gamma(double shape, double scale);
double init_beta(double a, double b); 
void partial_shuffle(vector<double>& vec, int start, int end);

#define DATADIR "../data/"
#define OUTDIR "../output/"
#define CONFIG "../$config/"
#define CONFIG_POP "../$config/pop/"
#define TRAN_PARAM "tran_params.csv"

#define GROUP_DATA                  "groups.csv"

#define EXPOSURE_AGE                "exposure_age.csv"

#define BIRTH_FILE                  "birth_rates.csv"
#define MORTALITY_FILE              "mortality_rates.csv"
#define AGE_BRACKETS                "pop_age_dist.csv"

#define CROW_DISTANCE               "euc_dist.csv"
#define CAR_DISTANCE                "road_dist.csv"

#define MDA_PARAMS                  "mda_params.csv"
#define INIT_PARAMS                 "init_params.csv"
#define CLUSTERING_PARAMS           "clustering_params.csv"
#define COUNTRY_PARAMS              "country_params.csv"
#endif /* headers_h */