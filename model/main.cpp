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

    string config_path = string(DATADIR) + COUNTRY_CONFIG;
    CountryConfig cfg = load_country_config(config_path);

    Region *rgn = new Region(region_id, region_name);
    rgn->start_year     = cfg.init_year;
    rgn->sim_years      = cfg.sim_years();
    rgn->init_ant_prev          = cfg.init_ant_prev;
    rgn->init_ant_prev_min      = cfg.init_ant_prev_min;
    rgn->init_ant_prev_max      = cfg.init_ant_prev_max;
    rgn->use_mf_prev_mode       = cfg.use_mf_prev_mode;
    rgn->init_mf_prev_min       = cfg.init_mf_prev_min;
    rgn->init_mf_prev_max       = cfg.init_mf_prev_max;
    rgn->init_mf_ant_ratio_min  = cfg.init_mf_ant_ratio_min;
    rgn->init_mf_ant_ratio_max  = cfg.init_mf_ant_ratio_max;
    rgn->achieved_coverage.assign(rgn->sim_years, 0.0);
    rgn->number_treated.assign(rgn->sim_years, 0);

    for (int i = 0; i < n_sims; ++i){

        rgn->reset_population();

        for (int year = 0; year < rgn->sim_years; ++year){
            cout << "Sim " << i+1 << "/" << n_sims << "  year " << cfg.init_year + year << endl;
            MDAEvent evt = make_mda_event(cfg, cfg.init_year + year);
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
            rgn,
            cfg
        );
    #endif

    return 0;
}
