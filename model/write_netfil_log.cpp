#include "write_netfil_log.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <filesystem>

using namespace std;

void write_section(ostream& out, const string& section_name) {
    out << "\n" << section_name << endl;
    out << string(section_name.length(), '-') << endl;
}

template <typename T>
void write_value(ostream& out, const string& value_name, const T& value) {
    out << value_name << ": " << value << endl;
}

void write_netfil(
    const string& filename,
    time_t start_time,
    time_t end_time,
    clock_t cpu_start,
    clock_t cpu_end,
    Region *rgn
) {
    string basename = filename;
    if (basename.size() >= 4 && basename.substr(basename.size() - 4) == ".csv") {
        basename = basename.substr(0, basename.size() - 4);
    }
    string netfil_name = basename + ".netfil";
    ofstream netfil(netfil_name);

    if (!netfil) { // Didn't open
        std::cerr << "Could not open .netfil file to write in";
        return;
    }

    write_section(netfil, "OUTPUT FILES");
    filesystem::path csv_path = filesystem::absolute(filename).lexically_normal();
    filesystem::path netfil_path = filesystem::absolute(netfil_name).lexically_normal();
    write_value(netfil, "Output folder", csv_path.parent_path().string());
    write_value(netfil, "CSV file", csv_path.filename().string());
    write_value(netfil, "Netfil log file", netfil_path.filename().string());

    write_section(netfil, "Parameters of last simulation number");
    
    ifstream state_file(string(DATADIR) + "current_state.txt");
    string state_country, state_scale, state_theta2, line;
    while (getline(state_file, line)) {
        auto eq = line.find('=');
        if (eq == string::npos) continue;
        string key = line.substr(0, eq);
        string val = line.substr(eq + 1);
        if      (key == "country") state_country = val;
        else if (key == "scale")   state_scale   = val;
        else if (key == "theta2")  state_theta2  = val;
    }
    state_file.close();
    write_value(netfil, "Country", state_country);
    write_value(netfil, "Scale",   state_scale);
    write_value(netfil, "Theta2",  state_theta2);
    write_value(netfil, "No. of groups", rgn->group_blocks);

    write_value(netfil, "Region name", rgn->rname);
    write_value(netfil, "Theta1", rgn->theta1);
    write_value(netfil, "Theta2", rgn->theta2);
    write_value(netfil, "Theta3", rgn->theta3);
    write_value(netfil, "Work:Home biting ratio", rgn->worktonot);
    write_value(netfil, "Agg param", rgn->agg_param);
    write_value(netfil, "Agg scale", rgn->agg_scale);
    write_value(netfil, "Init beta b", rgn->init_beta_b);
    write_value(netfil, "Init poisson", rgn->init_poisson);
    

    write_section(netfil, "MDA rounds");
    for (size_t i = 0; i < rgn->mda_rounds.size(); ++i) {
        const MDARound& r = rgn->mda_rounds[i];
        write_value(netfil, "Round",    i + 1);
        write_value(netfil, "Year",     r.year);
        write_value(netfil, "Drug",     r.drug);
        write_value(netfil, "Coverage", r.coverage);
        write_value(netfil, "Min age",  r.min_age);
        netfil << endl;
    }

    write_section(netfil, "Drug parameters");

    for (const auto& [name, drug] : FIXED_DRUG_DICT) {
        write_value(netfil, "Drug", name);
        drug.print_drugs(netfil);
        netfil << endl;
    }
    // DA/IDA aren't in FIXED_DRUG_DICT — their kill/sterilise split depends
    // on rgn->ster_to_kill (network.h), so compute them with this run's
    // actual value rather than logging a fixed placeholder.
    write_value(netfil, "Drug", "DA/IDA");
    write_value(netfil, "ster_to_kill", rgn->ster_to_kill);
    write_value(netfil, "mda_total_effect", rgn->mda_total_effect);
    lookup_drug("DA", rgn->ster_to_kill, rgn->mda_total_effect).print_drugs(netfil);
    netfil << endl;

    write_section(netfil, "Year parameters");
    write_value(netfil, "Starting year of simulation", rgn->init_year);
    write_value(netfil, "Ending year of simulation",   rgn->end_year);
    write_value(netfil, "No. years simulated",         rgn->sim_years);
    
    write_section(netfil, "Worm parameters");
    write_value(netfil, "Prop worms that are male", rgn->proportion_male_worm);
    write_value(netfil, "Immature period mean", rgn->immature_period_mean);
    write_value(netfil, "Immature period stdev", rgn->immature_period_mean_std);
    write_value(netfil, "Mature period mean", rgn->mature_period_mean);
    write_value(netfil, "Mature period stdev", rgn->mature_period_mean_std);

    write_section(netfil, "Human parameters");
    write_value(netfil, "Prop humans that are male", rgn->proportion_male_agent);
    write_value(netfil, "Proportion of group >5yo that commute", rgn->commuting_prop);
    write_value(netfil, "Number of age groups", N_AGE_GROUPS);
    write_value(netfil, "Maximum age age upon init", N_AGE_GROUPS * WIDTH_AGE_GROUPS - 1);

    write_section(netfil, "Disease parameters");

    write_value(netfil, "Initial antigen prevalence mean", rgn->init_ant_prev);
    string init_ant_prev_range = to_string(rgn->init_ant_prev_min) + "—" + to_string(rgn->init_ant_prev_max);
    write_value(netfil, "Allowed initial ANT prevalence range", init_ant_prev_range);

    if (rgn->use_mf_prev_mode) {
        string mf_prev_range = to_string(rgn->init_mf_prev_min) + "—" + to_string(rgn->init_mf_prev_max);
        write_value(netfil, "Allowed initial MF prevalence range", mf_prev_range);
    } else {
        string ratio_range = to_string(rgn->init_mf_ant_ratio_min) + "—" + to_string(rgn->init_mf_ant_ratio_max);
        write_value(netfil, "Allowed initial MF:ANT ratio range", ratio_range);
    }

    write_value(netfil, "Probability you lose antigen each day", rgn->daily_prob_lose_ant);

    double seed_success_pct = 100.0 * rgn->n_sims_run / rgn->total_seed_attempts;
    write_value(netfil, "Seeding success rate (%)", seed_success_pct);
    write_value(netfil, "Total seeding attempts", rgn->total_seed_attempts);

    write_section(netfil, "Miscellaneous parameters");
    write_value(netfil, "Sigma_group (group-level stdev)", rgn->sigma_group);
    write_value(netfil, "Beta_0", rgn->beta_0);
    write_value(netfil, "Distance type (e euclidean, r road)", DISTANCE_TYPE);
    write_value(netfil, "Number of years till road network re-estimated", RECALC_YEARS);
       
    write_section(netfil, "TIMESTAMP");
    int wall_duration = (int)difftime(end_time, start_time);
    int wall_h = wall_duration / 3600;
    int wall_min = (wall_duration % 3600) / 60;
    int wall_sec = wall_duration % 60;

    int cpu_duration = (int)((cpu_end - cpu_start) / CLOCKS_PER_SEC);
    int cpu_h = cpu_duration / 3600;
    int cpu_min = (cpu_duration % 3600) / 60;
    int cpu_sec = cpu_duration % 60;

    netfil << "Started: " << ctime(&start_time);
    netfil << "Ended: " << ctime(&end_time);
    netfil << "Wall time (h:m:s): " << wall_h << ":" << wall_min << ":" << wall_sec << endl;
    netfil << "CPU time  (h:m:s): " << cpu_h << ":" << cpu_min << ":" << cpu_sec << endl;

    netfil.close();
}
