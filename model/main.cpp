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
    if (argc < 2) { cerr << "Usage: ./main <output_name> [n_sims] [theta1 theta2 k w2n] [sigma_group beta_0] [ant_loss_halflife_days] [p_mda] [ster_to_kill] [mda_total_effect]\n"; return 1; }

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

    // Optional AntLoss override for ABC fitting (e.g. WSM), expressed as an
    // antigen half-life in days rather than the raw daily retention
    // probability — converted once here via daily_retention_from_halflife()
    // (network.h).
    bool cli_ant_loss = (argc >= 10);
    double cli_daily_prob_lose_ant = 0;
    if (cli_ant_loss) {
        cli_daily_prob_lose_ant = daily_retention_from_halflife(atof(argv[9]));
    }

    // Optional p_mda override for ABC fitting, same reasoning as cli_ant_loss
    // above. Position 10 is only meaningful because a run that wants p_mda
    // (or anything after it) must also send AntLoss at position 9 — the R
    // driver's particle loop always sends all three trailing values
    // together (ABC-fit or fixed-default) to guarantee this.
    bool cli_p_mda = (argc >= 11);
    double cli_p_mda_val = 0;
    if (cli_p_mda) {
        cli_p_mda_val = atof(argv[10]);
    }

    // Optional ster_to_kill override for ABC fitting, same reasoning — needs
    // positions 9 and 10 present too (see above).
    bool cli_ster_to_kill = (argc >= 12);
    double cli_ster_to_kill_val = 0;
    if (cli_ster_to_kill) {
        cli_ster_to_kill_val = atof(argv[11]);
    }

    // Optional mda_total_effect override for ABC fitting, same reasoning —
    // needs positions 9-11 present too (see above).
    bool cli_mda_total_effect = (argc >= 13);
    double cli_mda_total_effect_val = 0;
    if (cli_mda_total_effect) {
        cli_mda_total_effect_val = atof(argv[12]);
    }

    Region *rgn = new Region(region_id, region_name);
    rgn->start_year = rgn->init_year;
    rgn->achieved_coverage.assign(rgn->sim_years, 0.0);
    rgn->number_treated.assign(rgn->sim_years, 0);

    if (cli_ant_loss) {
        rgn->daily_prob_lose_ant = cli_daily_prob_lose_ant;
    }

    if (cli_p_mda) {
        rgn->p_mda = cli_p_mda_val;
    }

    if (cli_ster_to_kill) {
        rgn->ster_to_kill = cli_ster_to_kill_val;
    }

    if (cli_mda_total_effect) {
        rgn->mda_total_effect = cli_mda_total_effect_val;
    }

    for (int i = 0; i < n_sims; ++i){

        // Apply CLI overrides BEFORE reset_population(): pop_reload() (called by
        // reset_population()) constructs every Agent using the current agg_param
        // (Agent::Agent() bakes bite_scale = gamma(agg_param, 1/agg_param) in
        // permanently at construction). If reset_population() ran first, the very
        // first rep would build its whole population off whatever agg_param was
        // last set to (the tran_params.csv file default, from the Region
        // constructor's read_parameters() call) instead of this particle's k —
        // silently corrupting rep 0's transmission dynamics for its entire run.
        // Later reps "worked" only by accident, inheriting the previous rep's
        // override. Re-applied again below since reset_population()'s own
        // read_parameters() call overwrites these from file afterward.
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

        rgn->reset_population();

        // ant_loss/p_mda/ster_to_kill need the same reapplication as
        // theta/cluster above: reset_population()'s read_parameters() call
        // re-reads tran_params.csv every rep, including these as optional
        // columns 5-7 — a staged file that carries them (e.g. this
        // country/scale's own previous Fitted output) would silently
        // clobber the CLI value back to the file's on every rep after the
        // first if these weren't reapplied here too.
        if (cli_ant_loss) {
            rgn->daily_prob_lose_ant = cli_daily_prob_lose_ant;
        }

        if (cli_p_mda) {
            rgn->p_mda = cli_p_mda_val;
        }

        if (cli_ster_to_kill) {
            rgn->ster_to_kill = cli_ster_to_kill_val;
        }

        if (cli_mda_total_effect) {
            rgn->mda_total_effect = cli_mda_total_effect_val;
        }

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
