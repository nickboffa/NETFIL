#ifndef network_hpp
#define network_hpp
#include <array>
#include <cmath>
#include <string>

#include "mda.h"
#include "agent.h"

using namespace std;

// CLI/tran_params.csv express AntLoss as an antigen half-life in days (more
// interpretable/fittable than the raw daily retention probability) — this
// converts once at the point of assignment; daily_prob_lose_ant itself stays
// the raw per-day probability everywhere it's actually used (statistics.cpp).
inline double daily_retention_from_halflife(double half_life_days) {
    return pow(0.5, 1.0 / half_life_days);
}

class Group;                           // Groups of people akin to villages
class Region;                          // Region which is comprised of the groups!

class Group{
public:

    int gid;                           // Group ID
    
    double day_strength;            // Strength of infection during the day
    double night_strength;           // Strength of infection during the day
    
    double day_bites;
    double night_bites;

    double lat, lon;                    // Latitude & longitude
    Region *rgn;                       // Region
    double sum_mf;                      // NEED TO DEFINE

    map<int, Agent*> group_pop;       // Group population (out of work hours)

    //commuting data
    struct c_node{ //used to store distances to all other groups from current group
        int gid; //other group idea
        double dis; //the distance!
        c_node(int gid, double dis): gid(gid), dis(dis) {}
    };

    double total_commute = 0; //total commuters

    vector<c_node*> commuting_dist; //storing the distances
    map<int, double> commuting_pop; //the prop of commuters from each location
    map<int,double> commuting_cumsum; //cumsum of commuters from each location
    map<int, Agent*> day_population;    //the commuters to current group! and agents from group that did not commute!

    void add_member(Agent *agt);
    void rmv_member(Agent *agt);
    
    void bld_group_pop();  //build initial population
  

    Group(int gid, Region *rgn, double lat, double lon);
    ~Group();
};

class Region{
public:
    int rid;                           //region ID
    string rname;                      //region name
    int rpop;                          //region population
    int next_aid;                      //agent ID tracker for births
    bool init;                         // Has the population been built before?    
    
    double theta1;                      //transmission parameters for the different mf maturation scalings!
    double theta2;
    double theta3;

    double worktonot = 0;               //where do a majortiy of bites occur?

    double mf_to_ant_2014; //used to save fitting data

    double agg_param;
    double agg_scale;

    double sigma_group;
    double beta_0;

    // Computed by seed_lf() and used in the acceptance loop:
    double init_ant_prev_obs;    // observed ANT prevalence  (proportion 0–1, = ant_pos/rpop)
    double init_mf_prev_obs;     // observed MF prevalence   (proportion 0–1, = inf_indiv/rpop)
    double init_ratio_obs;       // observed MF:ANT ratio    (= inf_indiv/ant_pos)

    int age_dist_lower[N_AGE_GROUPS];
    int age_dist_upper[N_AGE_GROUPS];
    //used to keep track of total population for easy analysis
    map<int, Agent*> pre_indiv;        //collection of immautre worms individuals
    map<int, Agent*> inf_indiv;        //collection ofinfectious individuals
    map<int, Agent*> uninf_indiv;      //collection of peple with adult worms but are uninfectious individuals (single gender or sterile)
    map<int, Agent*> no_worms_indiv;   //collection of people with no worms!
    vector<Agent*> pvec[N_AGE_GROUPS]; //storing all people of certain age group
    vector<double> cum_sum_prob_worm {};
    //now all the information about the groups
    int next_gid, group_blocks;
    map<int, Group*> groups;            //storing all groups in region
    map<string, int> group_names;       //each group assigned name to index
    map<int, string> group_numbers;     //each group assigned number to index
    map<int, double*> group_coords;     //coords of each group

    // For fitting
    double mf_2014 = 0;
    double ant_2014 = 0;

    double mf_2016 = 0;
    double ant_2016 = 0;

    //distances 
    double *euclid_dst;                 //euclidean (L2) distance between groups
    double *road_dst;                   //road (L1) distance between groups

    map<int, int> group_pops;           //pop in each group
 
    // Simulation years (loaded from country.json via load_config())
    int init_year;
    int end_year;
    int sim_years;
    int end_day = 363;  // optional "end_day" in country.json: day-of-year (0-363) to
                         // stop at within the final simulated year, skipping the rest
                         // of that year. Defaults to 363 (simulate the full final year).

    // Seeding
    double init_ant_prev;      // target ANT prevalence for initial seeding (proportion 0–1)
    double init_ant_prev_min;  // acceptance lower bound for seeded ANT prevalence (proportion 0–1)
    double init_ant_prev_max;  // acceptance upper bound for seeded ANT prevalence (proportion 0–1)
    // MF seeding mode set by which JSON keys are present in "seeding":
    //   init_mf_prev_*      → use_mf_prev_mode = true  (bounds on MF prevalence = inf/rpop)
    //   init_mf_ant_ratio_* → use_mf_prev_mode = false (bounds on MF:ANT ratio  = inf/ant_pos)
    bool   use_mf_prev_mode;
    double init_mf_prev            = 0;
    double init_mf_prev_min        = 0;
    double init_mf_prev_max        = 0;
    double init_mf_ant_ratio_min   = 0;
    double init_mf_ant_ratio_max   = 0;

    vector<MDARound> mda_rounds;

    // From default_params.json (any of these can be overridden in country.json)
    double init_beta_b;           // worm lifespan shape parameter
    double init_poisson;          // mean immature worm count at seeding
    double immature_to_antigen;   // ratio of prepatent-only to antigen-positive at seeding
    double immature_and_ant;      // proportion of antigen-positive with co-occurring immature worms
    array<double, 16> mortality_rates;  // daily mortality rate per 5-year age bracket
    array<double, 16> birth_rates;      // daily birth rate per 5-year age bracket
    array<double, 16> exposure_kernel;  // mosquito exposure kernel by age bracket
    array<double, 17> age_dist;         // cumulative age distribution CDF (17 boundary points)

    double proportion_male_worm;        // proportion of new worms that are male
    double proportion_male_agent;       // proportion of agents that are male
    double daily_prob_lose_ant;         // daily probability of RETAINING antigen positivity — set via daily_retention_from_halflife(), not directly from CLI/file (which express this as a half-life in days)
    double immature_period_mean;        // mean immature worm period (days)
    double immature_period_mean_std;    // std dev of immature worm period (days)
    double mature_period_mean;          // mean mature worm lifespan (days)
    double mature_period_mean_std;      // std dev of mature worm lifespan (days)
    double commuting_prop;              // proportion of >5yo agents that commute daily
    double p_mda = 1.0;                 // effective-coverage discount: actual coverage = country.json's round coverage * p_mda
    // Re-splits DA/IDA's combined kill+sterilise effect (mda_total_effect
    // below) between the two: kill_prob = (1-ster_to_kill) * total,
    // full_ster_prob = ster_to_kill * total. Default reproduces the
    // literature-placeholder split (0.246 sterilise / 0.649 total).
    double ster_to_kill = 0.246 / 0.649;
    // Combined probability that MDA either kills or sterilises a DA/IDA
    // worm at all — the fraction NOT covered by this remains fully fertile
    // post-treatment. Default is the Hardy et al. literature placeholder
    // (see mda.h's DA_IDA_TOTAL_EFFECT); distinct from ster_to_kill, which
    // only re-splits this total rather than changing how many worms it covers.
    double mda_total_effect = 0.649;

    int start_year;
    int n_sims_run = 0;             // set in main() after all sims complete
    int total_seed_attempts = 0;    // accumulated across all runs in sim()
    vector<double> achieved_coverage;  // indexed by model year (0-based)
    vector<int>    number_treated;

    Region(int rid, string rname);

    void load_config();
    MDAEvent make_mda_event(int calendar_year);

    //Functions that run on region
    void sim(int year, MDAEvent evt);                           //wrapper to run simulation
    void handle_commute(int year);                               // generate commuter network and assign
    void remove_agent(Agent *agt);                                   //remove dead people from population
    void radt_model(char m);                                    //radiation model for daily trips (work/school)
    //void hndl_migrt(int day);                                //TODO long term migration between groups (to help avoid groups that have died out)
    void renew_pop(int year, int day, int dt);
    void handle_birth(int year, int day, int dt);                         // handle new births
    void calc_risk();         //find prevalence in each village
    void update_epi_status(int year, int day, int dt);                  //update agent's epi status
    void seed_lf();                                             //seed LF in population
    double mf_functional_form(char form, double worm_strength);            //converts worm strength to mf load

    void implement_mda(int year, MDAEvent evt);              //MDA!
    
    bool pop_reload();
    void read_groups();                                 //read input data
    void bld_groups();                                  //build the model groups 
    void bld_region_population();//build the population of the region
    void read_parameters();

    void reset_population();
    void reset_prev();
    void output_epidemics(int year, int day, MDAEvent evt);          //output outbreak data
    int factorial(int n);
    int n_worms();
    void prob_worms(double agg_param_init, double worm_mean);
};
#endif /* network_hpp */
