#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include "main.h"
#include "mda.h"
#include "strategy.h"
#include "rng.h"

using namespace std;

int SimulationNumber = 0;
string prv_out_loc;
static volatile sig_atomic_t abort_requested = 0;

static void handle_sigint(int) {
    abort_requested = 1;
}

// Remove cached population so next boot rebuilds from CSVs
static void reset_population_cache() {
    // NOTE: paths contain a literal "$" per your current params.h
    system("rm -f '../$config/ABC.init' 2>/dev/null");
    system("rm -rf '../$config/pop' 2>/dev/null");
}

int main(int argc, const char * argv[]) {
    // SIGINT handler to allow safe stop + cache reset
    std::signal(SIGINT, handle_sigint);

    // Keep legacy per-year outputs working
    prv_out_loc = (argc >= 2) ? argv[1] : (std::string)"out.csv";

    // Fixed settings
    const int    TAS_N     = 1500;    // TAS sample size (15+ antigen)
    const int    REPEATS   = 10;      // 10 repeats per combo
    const int    END_YEAR  = 2040;    // kept here for clarity; used in strategy.cpp

    // Proactive: start with a fresh population cache each program run
    reset_population_cache();

    // Build region
    region *rgn = new region(region_id, region_name);
    cout << "BOOT groups=" << rgn->groups.size()
         << " rpop=" << rgn->rpop << endl;

    // Scenarios (drug/coverage) file path
    string mda_data = datadir; mda_data = mda_data + MDA_params;

    string summary_path = outdir + string("strategy_results.csv");
    {
        bool write_header = false;
        std::ifstream test(summary_path.c_str(), std::ios::ate);
        if (!test.good() || test.tellg() == 0) write_header = true;
        test.close();

        if (write_header) {
            ofstream hdr(summary_path.c_str());
            hdr << "scenario,rep,seed,groups,pop,area_type,"
                "p_MDA_percent,p_covered,Y_TAS,tas_sample,"
                "mf_prev_2010_end,ant_prev_2010_end_allages,"
                "years_to_elimination,total_infected_people,total_years_with_infection,"
                "total_mda_doses,mda_years,tas_surveys,mf_zero_at_end,final_year\n";
        }
    }

    // Count MDA scenarios (rows in MDAParams.csv)
    int MDAScenario_count = count_mda_scenarios(mda_data);
    cout << "There are " << MDAScenario_count << " scenarios" << endl;

    // Area types to test: AS-wide and Subvillage only
    struct AreaOpt { const char* label; AreaType type; int offset; };
    AreaOpt AREAS[2] = {
        {"AS",        AreaType::AS_ALL,        0},
        {"subvillage",AreaType::SUBVILLAGE,  200}
    };

    // Grids (per your update)
    const double P_MDA_LIST[]   = {0.1, 0.5, 1}; // percent thresholds
    const double COVERED_LIST[] = {0.50, 0.65, 0.80};   // MDA coverage (0–1)
    const int    Y_TAS_LIST[]   = {2, 3};               // TAS cadence (years)

    for (int scenario_count = 0; scenario_count < MDAScenario_count; ++scenario_count) {
        // Read drug/etc. for this scenario and disable any CSV-scheduled rounds
        mda_strat base = get_mda_strat(mda_data, scenario_count + 1);
        base.NumRounds = 0;
        base.YearsBetweenRounds = 0;
        base.MDAYears.clear();

        for (const auto& A : AREAS) {
            for (double p_MDA : P_MDA_LIST) {
                for (double p_cov : COVERED_LIST) {
                    for (int Y_TAS : Y_TAS_LIST) {

                        // Early-abort check (before starting a block)
                        if (abort_requested) {
                            cerr << "\n[ABORT] Stopping on user request. Resetting caches...\n";
                            reset_population_cache();
                            return 130;
                        }

                        // Set coverage for this combo
                        mda_strat strat = base;
                        strat.Coverage = p_cov;  // override to desired p_covered

                        // Run REPEATS repeats
                        for (int rep = 1; rep <= REPEATS; ++rep) {
                            if (abort_requested) {
                                cerr << "\n[ABORT] Stopping on user request. Resetting caches...\n";
                                reset_population_cache();
                                return 130;
                            }

                            // Distinct seed per (scenario, area, p_MDA, p_cov, Y_TAS, rep)
                            const unsigned int BASE_SEED = 123456u;
                            unsigned int seed = BASE_SEED
                                              + (unsigned int)(scenario_count * 1000000)
                                              + (unsigned int)(A.offset)
                                              + (unsigned int)llround(p_MDA * 1000.0) * 1000u
                                              + (unsigned int)llround(p_cov * 100.0) * 10u
                                              + (unsigned int)(Y_TAS * 100000u)
                                              + (unsigned int)(rep);
                            reseed_rng(seed);

                            // Fresh population per run
                            rgn->reset_population();

                            // Strategy inputs
                            StrategyInputs in;
                            in.tas_every_Y_years = Y_TAS;
                            in.tas_sample_size   = TAS_N;
                            in.start_threshold_p = p_MDA;  // p_MDA in percent
                            in.max_years         = (END_YEAR - start_year + 1);
                            in.area_type         = A.type;

                            // Run controller (no early stop except when MF==0)
                            StrategyOutputs out;
                            run_elimination_strategy(*rgn, strat, in, out);

                            // Append a row to summary (single CSV for everything)
                            ofstream f(summary_path.c_str(), ios::app);
                            f << (scenario_count+1) << ","
                            << rep << ","
                            << seed << ","
                            << rgn->groups.size() << ","
                            << rgn->rpop << ","
                            << A.label << ","
                            << in.start_threshold_p << ","
                            << strat.Coverage << ","
                            << in.tas_every_Y_years << ","
                            << in.tas_sample_size << ","
                            << out.mf_prev_2010_end << ","                 // NEW
                            << out.ant_prev_2010_end_allages << ","        // NEW
                            << out.years_to_elimination << ","
                            << out.total_infected_people << ","
                            << out.total_years_with_infection << ","
                            << out.total_mda_doses << ","
                            << out.mda_years << ","
                            << out.tas_surveys << ","
                            << out.mf_zero_at_end << ","
                            << out.final_year << "\n";

                            // f closes → data flushed per row
                        }

                        // Flush after each finished block (area, p_MDA, p_cov, Y_TAS)
                        {
                            ofstream ff(summary_path.c_str(), ios::app);
                            ff << std::flush;
                        }
                        cout << "Completed: area=" << A.label
                             << " p_MDA=" << p_MDA
                             << "% p_covered=" << p_cov
                             << " Y_TAS=" << Y_TAS
                             << " (10 repeats)" << endl;
                    }
                }
            }
        }
    }

    return 0;
}
