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

inline const map<string, Drugs> FIXED_DRUG_DICT = {
    {"Old",   Drugs(0.55,  0.45,  0, 100, 0.5)},
    {"MoxDA", Drugs(0.457, 0.514, 0, 100, 0.5)}
};

// DA/IDA (Hardy et al: clearance rate (27 + 45) / (39 + 72) = 64.9%) share
// this combined kill+sterilise total, but the literature doesn't pin down
// how it splits between the two — Region::ster_to_kill (network.h) re-splits
// it and is fit via ABC instead of trusting the old placeholder split
// (kill_prob=0.403, full_ster_prob=0.246, i.e. ster_to_kill=0.246/0.649).
// DA_IDA_TOTAL_EFFECT is just the literature-placeholder default for
// Region::mda_total_effect (network.h) — the total itself is also fit via
// ABC now, distinct from ster_to_kill: total_effect controls how many worms
// escape MDA completely unaffected (still fertile) vs. either killed or
// sterilised, while ster_to_kill only re-splits the affected fraction.
inline constexpr double DA_IDA_TOTAL_EFFECT = 0.649;

inline Drugs da_ida_drugs(double ster_to_kill, double total_effect) {
    return Drugs(
        (1.0 - ster_to_kill) * total_effect, // kill_prob
        ster_to_kill * total_effect,         // full_ster_prob
        0, 100, 0.5
    );
}

inline Drugs lookup_drug(const string& name, double ster_to_kill, double total_effect) {
    if (name == "DA" || name == "IDA") return da_ida_drugs(ster_to_kill, total_effect);
    return FIXED_DRUG_DICT.at(name);
}

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
