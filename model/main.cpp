#include <iostream>
#include <fstream>
#include <string>
#include "main.h"
#include "mda.h"
#include "strategy.h"
#include "rng.h"

using namespace std;

int SimulationNumber = 0;
string prv_out_loc;

int main(int argc, const char * argv[]) {
    // Keep your legacy per-year outputs working
    prv_out_loc = (argc >= 2) ? argv[1] : (std::string)"out.csv";

    // Fixed settings per your request
    const int    K_TAS     = 2;       // TAS every 2 years
    const double P_THRESH  = 1.0;     // elimination threshold (%)
    const int    TAS_N     = 1500;    // TAS sample size (15+ antigen)
    const int    REPEATS   = 50;       // ***** 5 repeats per setting *****
    const int    END_YEAR  = 2050;    // always run through 2030 (unless MF==0)

    // Build region
    region *rgn = new region(region_id, region_name);
    cout << "BOOT groups=" << rgn->groups.size()
         << " rpop=" << rgn->rpop << endl;

    // Scenarios (drug/coverage) file path
    string mda_data = datadir; mda_data = mda_data + MDA_params;

    // Prepare summary CSV (overwrite header once)
    string summary_path = outdir + string("strategy_results.csv");
    {
        ofstream hdr(summary_path.c_str());
        hdr << "scenario,rep,seed,groups,pop,k,area_type,p_threshold,tas_sample,"
               "years_to_elimination,total_infected_people,total_years_with_infection,total_mda_doses,"
               "mda_years,tas_surveys,mf_zero_at_end\n";
    }

    // Count MDA scenarios (rows in MDAParams.csv)
    int MDAScenario_count = count_mda_scenarios(mda_data);
    cout << "There are " << MDAScenario_count << " scenarios" << endl;

    // Area types to test
    struct AreaOpt { const char* label; AreaType type; int offset; };
    AreaOpt AREAS[3] = {
        {"AS",        AreaType::AS_ALL,        0},
        {"village",   AreaType::VILLAGE,     100},
        {"subvillage",AreaType::SUBVILLAGE,  200}
    };

    for (int scenario_count = 0; scenario_count < MDAScenario_count; ++scenario_count) {
        // Read drug/coverage for this scenario
        mda_strat base = get_mda_strat(mda_data, scenario_count + 1);

        // OPTION A: disable CSV-scheduled rounds so sim.cpp never auto-fires MDA
        base.NumRounds = 0;
        base.YearsBetweenRounds = 0;
        base.MDAYears.clear();

        for (const auto& A : AREAS) {
            for (int rep = 1; rep <= REPEATS; ++rep) {

                // Distinct seed per (scenario, area, rep)
                const unsigned int BASE_SEED = 123456u;
                unsigned int seed = BASE_SEED
                                  + (unsigned int)(scenario_count * 10000)
                                  + (unsigned int)(A.offset)
                                  + (unsigned int)(rep);
                reseed_rng(seed);

                // Fresh population per run
                rgn->reset_population();

                // Strategy inputs
                StrategyInputs in;
                in.tas_every_k_years = K_TAS;
                in.tas_sample_size   = TAS_N;
                in.start_threshold_p = P_THRESH;
                in.max_years         = (END_YEAR - start_year + 1); // 2010..2030 inclusive
                in.area_type         = A.type;

                // Run controller (no early stop except when MF==0)
                StrategyOutputs out;
                run_elimination_strategy(*rgn, base, in, out);

                // Append a row to summary
                ofstream f(summary_path.c_str(), ios::app);
                f << (scenario_count+1) << ","
                  << rep << ","
                  << seed << ","
                  << rgn->groups.size() << ","
                  << rgn->rpop << ","
                  << in.tas_every_k_years << ","
                  << A.label << ","
                  << in.start_threshold_p << ","
                  << in.tas_sample_size << ","
                  << out.years_to_elimination << ","
                  << out.total_infected_people << ","
                  << out.total_years_with_infection << ","
                  << out.total_mda_doses << ","
                  << out.mda_years << ","
                  << out.tas_surveys << ","
                  << out.mf_zero_at_end << "\n";

                SimulationNumber += 1;
            }
        }
    }

    return 0;
}
