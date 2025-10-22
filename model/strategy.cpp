#include "strategy.h"
#include <numeric>
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

void run_elimination_strategy(region& rgn,
                              mda_strat base_mda,
                              const StrategyInputs& in,
                              StrategyOutputs& out)
{
    const int START_YEAR    = start_year;   // e.g., 2010 from params.h
    const int END_YEAR      = 2050;         // always run through 2030 (unless MF==0)
    const int MDA_START_ABS = 2016;         // no MDA before 2016

    // ---- Per-run resets (arrays are raw pointers, so set values explicitly) ----
    for (int y = 0; y < in.max_years; ++y) {
        rgn.number_treated[y]    = 0;
        rgn.achieved_coverage[y] = 0.0;
    }

    // ---- Per-run tallies ----
    out.total_infected_people      = 0;
    out.total_years_with_infection = 0;   // PERSON-YEARS: sum of I each year
    out.total_mda_doses            = 0;
    out.mda_years                  = 0;
    out.tas_surveys                = 0;
    out.mf_zero_at_end             = 0;

    std::vector<double> annual_prev;
    annual_prev.reserve((END_YEAR - START_YEAR) + 1);

    bool threshold_reached   = false;  // flips once region MF < p%
    int  tas_start_abs       = -1;     // ABS year TAS begins (first under-threshold year)
    std::map<int,bool> on_targeted_mda; // if true for gid, resume annual MDA there
    int  earliest_elim_abs   = -1;     // record-only (we do NOT stop early except MF==0)

    auto lin_slope = [](const std::vector<double>& y){
        int n = (int)y.size();
        if (n < 2) return 0.0;
        double sx=0, sy=0, sxx=0, sxy=0;
        for (int i=0; i<n; ++i) {
            double x = i;
            sx  += x;
            sy  += y[i];
            sxx += x*x;
            sxy += x*y[i];
        }
        double denom = n*sxx - sx*sx;
        return (denom == 0.0) ? 0.0 : (n*sxy - sx*sy) / denom;
    };

    for (int year = 0; year < in.max_years; ++year) {
        const int ABS_YEAR = START_YEAR + year;
        if (ABS_YEAR > END_YEAR) break;

        // ensure this year's slot is clean
        rgn.number_treated[year]    = 0;
        rgn.achieved_coverage[year] = 0.0;
        bool applied_mda_this_year  = false;

        const bool in_mda_window = (ABS_YEAR >= MDA_START_ABS);

        // --- 1) Decide dosing for this year ---
        if (in_mda_window && !threshold_reached) {
            // universal annual MDA until under threshold
            rgn.implement_MDA(year, base_mda);
            applied_mda_this_year = true;
        } else if (in_mda_window && threshold_reached) {
            // TAS phase: targeted MDA only for areas that previously failed TAS
            std::vector<int> resume_all;
            for (auto &kv : on_targeted_mda) {
                if (kv.second) resume_all.push_back(kv.first);
            }
            if (!resume_all.empty()) {
                rgn.implement_MDA_subset(year, base_mda, resume_all);
                applied_mda_this_year = true;
            }
        }
        // else: ABS_YEAR < 2016 → no MDA

        // --- 2) Run one epidemiological year ---
        rgn.sim(year, base_mda);

        // --- 3) Record outcomes ---
        double prev_mf_region = rgn.mf_prev_region();   // % at the end of the year
        annual_prev.push_back(prev_mf_region);

        // PERSON-YEARS infected + “total infected people” (year-end I)
        const long long inf_count = (long long)rgn.inf_indiv.size();
        out.total_infected_people      += inf_count;
        out.total_years_with_infection += inf_count;

        // Doses and MDA-year accounting
        out.total_mda_doses += (long long)rgn.number_treated[year];
        if (applied_mda_this_year && rgn.number_treated[year] > 0) {
            out.mda_years += 1;
        }

        // Update final-year zero flag as we go (last year overwrites)
        out.mf_zero_at_end = (prev_mf_region == 0.0) ? 1 : 0;

        // early stop for THIS RUN if MF==0 (move to next repeat)
        if (prev_mf_region == 0.0) {
            if (earliest_elim_abs < 0) earliest_elim_abs = ABS_YEAR;
            break;
        }

        // --- 4) Threshold crossing unlocks TAS (only from 2016 onwards) ---
        if (in_mda_window && !threshold_reached && (prev_mf_region < in.start_threshold_p)) {
            threshold_reached = true;
            tas_start_abs = ABS_YEAR;
            on_targeted_mda.clear();
            for (auto &grp : rgn.groups) {
                on_targeted_mda[grp.first] = false;
            }
        }

        // --- 5) TAS (15+ antigen) every k years AFTER threshold reached ---
        if (threshold_reached && tas_start_abs >= 0) {
            int delta = ABS_YEAR - tas_start_abs;
            if (delta > 0 && (delta % in.tas_every_k_years == 0)) {
                // Build area sets
                std::vector<std::vector<int>> sets;
                if (in.area_type == AreaType::AS_ALL) {
                    std::vector<int> all;
                    all.reserve(rgn.groups.size());
                    for (auto &kv : rgn.groups) all.push_back(kv.first);
                    sets.push_back(std::move(all));
                } else {
                    // NOTE: until you add a village-aggregation map, both village and
                    // subvillage are single-group units (i.e., identical granularity).
                    sets.reserve(rgn.groups.size());
                    for (auto &kv : rgn.groups) {
                        sets.push_back(std::vector<int>{kv.first});
                    }
                }

                // Reset flags for this TAS round
                on_targeted_mda.clear();
                for (auto &grp : rgn.groups) {
                    on_targeted_mda[grp.first] = false;
                }

                // Run TAS for each set, mark failing sets for targeted MDA
                for (auto &set : sets) {
                    double tas_prev = rgn.tas_estimate_antigen_15plus(set, in.tas_sample_size, /*current_year=*/year);
                    if (tas_prev > 1.0) {
                        for (int gid : set) {
                            on_targeted_mda[gid] = true;
                        }
                    }
                }

                // Count one survey round per TAS cadence year
                out.tas_surveys += 1;
            }
        }

        // --- 6) Elimination rule (record-only; NEVER stop early because of this) ---
        bool under_threshold_now = (prev_mf_region < in.start_threshold_p);
        bool nonpos_trend = false;
        {
            int n = (int)annual_prev.size();
            if (n >= 3) {
                int window = std::min(5, n);
                std::vector<double> tail(annual_prev.end() - window, annual_prev.end());
                nonpos_trend = (lin_slope(tail) <= 0.0);
            }
        }
        if (earliest_elim_abs < 0 && threshold_reached && under_threshold_now && nonpos_trend) {
            earliest_elim_abs = ABS_YEAR; // record but keep simming to 2030 unless MF==0
        }
    }

    // Fill output field (years to elimination relative to start_year), -1 if never
    if (earliest_elim_abs >= 0) {
        out.years_to_elimination = (earliest_elim_abs - START_YEAR + 1);
    } else {
        out.years_to_elimination = -1;
    }
}
