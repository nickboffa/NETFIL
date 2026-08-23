#include <iostream>
#include <ctime>
#include <cmath>
#include <unistd.h>

#include "main.h"
#include "mda.h"
#include "write_netfil_log.h"

using namespace std;

int sim_i = 0;

string prv_out_loc;

int main(int argc, const char * argv[]){
    if (argc < 2) { cerr << "Usage: ./main <output_name> [n_sims] [theta1 theta2 k w2n] [sigma_group beta_0]\n"; return 1; }

    time_t start_time = time(nullptr);
    clock_t cpu_start = clock();
    prv_out_loc = argv[1];
    int n_sims = (argc >= 3) ? atoi(argv[2]) : 1;

    // Transmission params optionally passed via CLI for parallel ABC runs.
    // If absent, read_parameters() loads them from tran_params.csv as normal.
    bool cli_tran = (argc >= 7);
    double cli_t1 = 0, cli_t2 = 0, cli_k = 0, cli_w = 0;
    if (cli_tran) {
        cli_t1 = atof(argv[3]);
        cli_t2 = atof(argv[4]);
        cli_k  = atof(argv[5]);
        cli_w  = atof(argv[6]);
    }

    // Clustering params (2010 seeding ICC, expressed as sigma_group/beta_0)
    // optionally passed via CLI so ABC can fit them per-particle, same reason
    // as cli_tran above — a per-particle clustering_params.csv would race
    // under parallel workers. If absent, read_parameters() loads them from
    // clustering_params.csv as normal.
    bool cli_cluster = (argc >= 9);
    double cli_sigma_group = 0, cli_beta_0 = 0;
    if (cli_cluster) {
        cli_sigma_group = atof(argv[7]);
        cli_beta_0      = atof(argv[8]);
    }

    Region *rgn = new Region(region_id, region_name);
    rgn->start_year = rgn->init_year;
    rgn->achieved_coverage.assign(rgn->sim_years, 0.0);
    rgn->number_treated.assign(rgn->sim_years, 0);

    for (int i = 0; i < n_sims; ++i){

        rgn->reset_population();

        if (cli_tran) {
            rgn->theta1    = cli_t1;
            rgn->theta2    = cli_t2;
            rgn->theta3    = 1.0 / (1.0 - exp(-cli_t2));
            rgn->agg_param = cli_k;
            rgn->agg_scale = 1.0 / cli_k;
            rgn->worktonot = cli_w;
        }

        if (cli_cluster) {
            rgn->sigma_group = cli_sigma_group;
            rgn->beta_0      = cli_beta_0;
        }

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
