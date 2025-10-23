#include "strategy.h"
#include "agent.h"   // for agent::is_antigen_positive()
#include <numeric>
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>
#include <unordered_set>

static double compute_antigen_prev_all_ages(region& rgn) {
    long long ant_pos = 0;
    long long total   = 0;
    for (auto &gkv : rgn.groups) {
        group* grp = gkv.second;
        for (auto &akv : grp->group_pop) {
            agent* a = akv.second;
            ++total;
            if (a->is_antigen_positive()) ++ant_pos;
        }
    }
    if (total == 0) return 0.0;
    return 100.0 * (double)ant_pos / (double)total;  // percent
}

void run_elimination_strategy(region& rgn,
                              mda_strat base_mda,
                              const StrategyInputs& in,
                              StrategyOutputs& out)
{
    const int START_YEAR    = start_year;   // e.g., 2010
    const int END_YEAR      = 2040;         // cap at 2040
    const int MDA_START_ABS = 2016;         // no MDA before 2016

    // Reset per-run tallies/arrays
    for (int y = 0; y < in.max_years; ++y) {
        rgn.number_treated[y]    = 0;
        rgn.achieved_coverage[y] = 0.0;
    }

    out.total_infected_people      = 0;
    out.total_years_with_infection = 0;
    out.total_mda_doses            = 0;
    out.mda_years                  = 0;
    out.tas_surveys                = 0;
    out.mf_zero_at_end             = 0;
    out.final_year                 = -1;
    out.mf_prev_2010_end           = -1.0;
    out.ant_prev_2010_end_allages  = -1.0;

    std::vector<double> annual_prev;
    annual_prev.reserve((END_YEAR - START_YEAR) + 1);

    // Track UNIQUE persons ever MF-positive across the run
    std::unordered_set<int> ever_infected_ids;

    bool threshold_reached = false;   // flips once regional MF < p_MDA
    int  tas_start_abs     = -1;      // ABS year when TAS schedule begins
    std::map<int,bool> on_targeted_mda; // gid -> needs targeted MDA next year?
    int  earliest_elim_abs = -1;      // first ABS year meeting elimination rule (record-only)

    auto lin_slope = [](const std::vector<double>& y){
        const int n = (int)y.size();
        if (n < 2) return 0.0;
        double sx=0, sy=0, sxx=0, sxy=0;
        for (int i=0; i<n; ++i) {
            const double x = i;
            sx  += x;  sy  += y[i];
            sxx += x*x; sxy += x*y[i];
        }
        const double denom = n*sxx - sx*sx;
        return (denom == 0.0) ? 0.0 : (n*sxy - sx*sy) / denom;
    };

    for (int year = 0; year < in.max_years; ++year) {
        const int ABS_YEAR = START_YEAR + year;
        if (ABS_YEAR > END_YEAR) break;

        rgn.number_treated[year]    = 0;
        rgn.achieved_coverage[year] = 0.0;
        bool applied_mda_this_year  = false;

        const bool in_mda_window = (ABS_YEAR >= MDA_START_ABS);

        // (1) Dosing decision this year
        if (in_mda_window && !threshold_reached) {
            // universal annual MDA until under threshold
            rgn.implement_MDA(year, base_mda);
            applied_mda_this_year = true;
        } else if (in_mda_window && threshold_reached) {
            // targeted MDA only for areas that previously failed TAS
            std::vector<int> resume_all;
            for (auto &kv : on_targeted_mda) if (kv.second) resume_all.push_back(kv.first);
            if (!resume_all.empty()) {
                rgn.implement_MDA_subset(year, base_mda, resume_all);
                applied_mda_this_year = true;
            }
        }

        // (2) Run one epidemiological year
        rgn.sim(year, base_mda);

        // (3) Record outcomes for the end of this ABS_YEAR
        const double prev_mf_region = rgn.mf_prev_region();   // % at end of year
        annual_prev.push_back(prev_mf_region);

        // If this was 2010, save the end-of-2010 prevalences
        if (ABS_YEAR == 2010) {
            out.mf_prev_2010_end = prev_mf_region;
            out.ant_prev_2010_end_allages = compute_antigen_prev_all_ages(rgn);
        }

        // Person-years: year-end MF snapshot
        const long long inf_count = (long long)rgn.inf_indiv.size();
        out.total_years_with_infection += inf_count;

        // Unique MF-positive persons ever seen
        for (const auto &kv : rgn.inf_indiv) {
            ever_infected_ids.insert(kv.first);
        }

        // Dose totals / counts
        out.total_mda_doses += (long long)rgn.number_treated[year];
        if (applied_mda_this_year && rgn.number_treated[year] > 0) {
            out.mda_years += 1;
        }

        out.mf_zero_at_end = (prev_mf_region == 0.0) ? 1 : 0;

        // MF==0 → record final_year and stop this run
        if (prev_mf_region == 0.0) {
            if (earliest_elim_abs < 0) earliest_elim_abs = ABS_YEAR;
            out.final_year = ABS_YEAR;    // first year MF hits exactly 0
            break;
        }

        // (4) Unlock TAS when we first cross the threshold (only from 2016)
        if (in_mda_window && !threshold_reached && (prev_mf_region < in.start_threshold_p)) {
            threshold_reached = true;
            tas_start_abs = ABS_YEAR;
            on_targeted_mda.clear();
            for (auto &grp : rgn.groups) on_targeted_mda[grp.first] = false;
        }

        // (5) TAS cadence after threshold reached
        if (threshold_reached && tas_start_abs >= 0) {
            const int delta = ABS_YEAR - tas_start_abs;
            if (delta > 0 && (delta % in.tas_every_Y_years == 0)) {
                // Build sets by area type
                std::vector<std::vector<int>> sets;
                if (in.area_type == AreaType::AS_ALL) {
                    std::vector<int> all; all.reserve(rgn.groups.size());
                    for (auto &kv : rgn.groups) all.push_back(kv.first);
                    sets.push_back(std::move(all));
                } else {
                    // No village aggregation yet; treat each group independently.
                    sets.reserve(rgn.groups.size());
                    for (auto &kv : rgn.groups) sets.push_back(std::vector<int>{kv.first});
                }

                on_targeted_mda.clear();
                for (auto &grp : rgn.groups) on_targeted_mda[grp.first] = false;

                // Run TAS for each set; mark failing sets for targeted MDA
                for (auto &set : sets) {
                    const double tas_prev =
                        rgn.tas_estimate_antigen_15plus(set, in.tas_sample_size, /*current_year=*/year);
                    if (tas_prev > 1.0) {
                        for (int gid : set) on_targeted_mda[gid] = true;
                    }
                }

                out.tas_surveys += 1;
            }
        }

        // (6) Elimination rule (record-only; never stop early because of this)
        const bool under_threshold_now = (prev_mf_region < in.start_threshold_p);
        bool nonpos_trend = false;
        {
            const int n = (int)annual_prev.size();
            if (n >= 3) {
                const int window = std::min(5, n);
                std::vector<double> tail(annual_prev.end() - window, annual_prev.end());
                nonpos_trend = (lin_slope(tail) <= 0.0);
            }
        }
        if (earliest_elim_abs < 0 && threshold_reached && under_threshold_now && nonpos_trend) {
            earliest_elim_abs = ABS_YEAR; // record only
        }

        // If we continue, update final_year to current ABS_YEAR
        out.final_year = ABS_YEAR;
    }

    // Fill summary outputs
    if (out.final_year < 0) out.final_year = END_YEAR;
    out.total_infected_people = (long long)ever_infected_ids.size();

    // years_to_elimination remains a record-only field, populated by earliest_elim_abs
    // (You can map earliest_elim_abs exactly if you want; left as-is for compatibility.)
}
