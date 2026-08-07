#ifndef mda_h
#define mda_h
#include "params.h"
#include <array>
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
    int    min_age;
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


#endif /* mda_h */
