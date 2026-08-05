#ifndef mda_h
#define mda_h
#include "params.h"
#include <map>
#include <string>
#include <vector>

using namespace std;

// ── Drug effects ─────────────────────────────────────────────────────────────
class Drugs {
public:
    double kill_prob;
    double full_ster_prob;
    double part_ster_prob;
    double ster_dur;
    double part_ster_magnitude;

    Drugs(double K=0, double FSP=0, double PSP=0, double SD=0, double PSM=0)
        : kill_prob(K), full_ster_prob(FSP), part_ster_prob(PSP),
          ster_dur(SD), part_ster_magnitude(PSM) {}

    void print_drugs(ostream& out = cout) const {
        out << "kill_prob: "           << kill_prob           << endl;
        out << "full_ster_prob: "      << full_ster_prob      << endl;
        out << "part_ster_prob: "      << part_ster_prob      << endl;
        out << "ster_dur: "            << ster_dur            << endl;
        out << "part_ster_magnitude: " << part_ster_magnitude << endl;
    }
};

// DA (DEC + Albendazole) parameters are placeholders — verify against
// calibrated literature values before running WSM simulations.
inline const map<string, Drugs> DRUG_DICT = {
    {"MoxDA", Drugs(0.457, 0.514, 0, 100, 0.5)},
    {"IDA",   Drugs(0.403, 0.246, 0, 100, 0.5)},
    {"DA",    Drugs(0.403, 0.246, 0, 100, 0.5)}
};

// Combines the Hardy et al IDA and DA data
// To use a clearance rate of (27 + 45) / (39 + 72) = 64.9% for both instead?

// ── One historical MDA round ──────────────────────────────────────────────────
struct MDARound {
    int    year;
    int    month    = 0;   // 0 = not specified → default to day 28
    string drug;
    double coverage = 0.0;
    int    min_age  = 2;
};

// ── Country-level simulation configuration ────────────────────────────────────
struct CountryConfig {
    int    init_year;
    int    end_year;
    // Seeding
    double init_ant_prev;      // target ANT prevalence for initial seeding (proportion 0–1)
    double init_ant_prev_min;  // acceptance lower bound for seeded ANT prevalence (proportion 0–1)
    double init_ant_prev_max;  // acceptance upper bound for seeded ANT prevalence (proportion 0–1)
    // MF seeding mode set by which JSON keys are present in "seeding":
    //   init_mf_prev_*      → use_mf_prev_mode = true  (bounds on MF prevalence = inf/rpop)
    //   init_mf_ant_ratio_* → use_mf_prev_mode = false (bounds on MF:ANT ratio  = inf/ant_pos)
    bool   use_mf_prev_mode;
    double init_mf_prev_min        = 0;  // set when use_mf_prev_mode == true
    double init_mf_prev_max        = 0;
    double init_mf_ant_ratio_min   = 0;  // set when use_mf_prev_mode == false
    double init_mf_ant_ratio_max   = 0;
    vector<MDARound> mda_rounds;

    int sim_years() const { return end_year - init_year; }
};

// ── Per-year MDA event passed to Region::sim() ────────────────────────────────
// apply=false means no MDA this year; all other fields are ignored.
struct MDAEvent {
    bool   apply       = false;
    int    day         = 28;   // day-of-year to apply MDA
    string drug_name   = "";
    Drugs  drug_params;
    double coverage    = 0.0;
    int    min_age     = 0;
};

CountryConfig load_country_config(const string& path);
MDAEvent      make_mda_event(const CountryConfig& cfg, int calendar_year);

#endif /* mda_h */
