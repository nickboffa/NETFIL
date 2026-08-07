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

// Potential improvement: infer number of age groups from pop_age_dists.csv?
constexpr int N_AGE_GROUPS    = 16; //number of 5-year age brackets (for seeding pop)
constexpr int WIDTH_AGE_GROUPS = 5; // 0-4, 5-9, ... 75-79

// ABC_FITTING must remain a #define — it is used in a preprocessor #if directive
#define ABC_FITTING true

constexpr int    RECALC_YEARS        = 100;           //how often we want to recalc commuters
constexpr char   DISTANCE_TYPE       = 'r';           // r for road distance, e for euclidean. If changed, update setup.sh to copy the matching distance files.



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

#define CROW_DISTANCE               "euc_dist.csv"
#define CAR_DISTANCE                "road_dist.csv"
#define CROW_DISTANCE_BIN           "euc_dist.bin"
#define CAR_DISTANCE_BIN            "road_dist.bin"

#define COUNTRY_CONFIG              "country.json"

#define CLUSTERING_PARAMS           "clustering_params.csv"
#endif /* headers_h */