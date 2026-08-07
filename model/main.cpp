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
    if (argc < 2) { cerr << "Usage: ./main <output_name> [n_sims]\n"; return 1; }

    time_t start_time = time(nullptr);
    clock_t cpu_start = clock();
    prv_out_loc = argv[1];
    int n_sims = (argc >= 3) ? atoi(argv[2]) : 1;

    Region *rgn = new Region(region_id, region_name);
    rgn->start_year = rgn->init_year;
    rgn->achieved_coverage.assign(rgn->sim_years, 0.0);
    rgn->number_treated.assign(rgn->sim_years, 0);

    for (int i = 0; i < n_sims; ++i){

        rgn->reset_population();

        for (int year = 0; year < rgn->sim_years; ++year){
            cout << "Sim " << i+1 << "/" << n_sims << "  year " << rgn->init_year + year << endl;
            MDAEvent evt = rgn->make_mda_event(rgn->init_year + year);
            rgn->sim(year, evt);
        }

        sim_i += 1;
    }
    rgn->n_sims_run = n_sims;

    time_t end_time = time(nullptr);
    clock_t cpu_end = clock();

    string filename = string(OUTDIR) + prv_out_loc;
    #if !ABC_FITTING
        write_netfil(
            filename,
            start_time,
            end_time,
            cpu_start,
            cpu_end,
            rgn
        );
    #endif

    return 0;
}
