#include <iostream>
#include <ctime>
#include <unistd.h>

#include "main.h"
#include "mda.h"
#include "write_netfil_log.h"

using namespace std;

int sim_i = 0;

string prv_out_loc;

int main(int argc, const char * argv[]){
    time_t start_time = time(nullptr);
    prv_out_loc = argv[1];
   
    Region *rgn = new Region(region_id, region_name);
    
    string mda_data = string(DATADIR) + MDA_PARAMS; // Both are #define macros

    //Counting the number of different simulations we will perform
    int mda_scenario_count = count_mda_scenarios(mda_data);
    // cout << "There are " << mda_scenario_count << " scenarios" << endl;
    
    //now looping over scenarios
    for (int scenario_count = 0; scenario_count < mda_scenario_count; ++scenario_count){

        //generating mda strategy!
        MDAStrat strategy = get_mda_strat(mda_data, scenario_count + 1);

        //Now looping over simulations
         for (int i = 0; i < strategy.n_sims; ++i){
            
            //resetting the populations from previous simulation
            rgn->reset_population();
           
            //run run the simulation year by year
            for(int year = 0; year < SIM_YEARS; ++year){
                
                rgn->sim(year, strategy);
              
            }

            sim_i += 1;
        }

    }

    time_t end_time = time(nullptr);

    string filename = string(OUTDIR) + prv_out_loc;
#if !ABC_FITTING
    write_netfil(
        filename,
        start_time,
        end_time,
        rgn,
        mda_data
    );
#endif

    return 0;
} 