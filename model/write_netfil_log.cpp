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
    Region *rgn,
    string mda_data
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
    
    // Outputs the scale of whichever make_<>.sh file was run most recently
    ifstream scale_file(string(datadir) + "current_scale.txt");
    string current_scale;
    getline(scale_file, current_scale);
    scale_file.close();
    write_value(netfil, "Scale", current_scale);
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
    

    write_section(netfil, "MDA parameters");
    int num_scenarios = count_mda_scenarios(mda_data);
    for (int i = 1; i <= num_scenarios; i++) {
        MDAStrat strategy = get_mda_strat(mda_data, i + 1);
        write_value(netfil, "Strategy number", i);
        strategy.print_mda_strat(netfil);
        netfil << endl;
    }

    write_section(netfil, "Year parameters");
    write_value(netfil, "Starting year of simulation", start_year);
    write_value(netfil, "Ending year of simulation", start_year+sim_years);
    write_value(netfil, "No. years simulated", sim_years);
    
    write_section(netfil, "Worm parameters");
    write_value(netfil, "Prop worms that are male", proportion_male_worm);
    write_value(netfil, "Immature period mean", immature_period_mean);
    write_value(netfil, "Immature period stdev", immature_period_mean_std);
    write_value(netfil, "Mature period mean", mature_period_mean);
    write_value(netfil, "Mature period stdev", mature_period_mean_std);

    write_section(netfil, "Human parameters");
    write_value(netfil, "Prop humans that are male", proportion_male_agent);
    write_value(netfil, "Proportion of group >5yo that commute", commuting_prop);
    write_value(netfil, "Number of age groups", n_age_groups);
    write_value(netfil, "Maximum age age upon init", max_init_age);

    write_section(netfil, "Disease parameters");

    write_value(netfil, "Initial antigen prevalence mean", ant_0);
    string init_prev_range = to_string(init_prev_min) + "—" + to_string(init_prev_max);
    write_value(netfil, "Allowed initial antigen prevalence range", init_prev_range);

    string init_ratio_range = to_string(init_ratio_min) + "—" + to_string(init_ratio_max);
    write_value(netfil, "Initial mf:ant ratio range(?)", init_ratio_range);

    write_value(netfil, "Probability you lose antigen each day", DailyProbLoseAntigen);

    write_section(netfil, "Miscellaneous parameters");
    write_value(netfil, "Sigma_g (household stdev)", sigma_g);
    write_value(netfil, "Beta_0", beta_0);
    write_value(netfil, "Distance type (e euclidean, r road)", distance_type);
    write_value(netfil, "Number of years till road network re-estimated", recalc_years);
       
    write_section(netfil, "TIMESTAMP");
    int duration = (int)difftime(end_time, start_time);
    int duration_h = duration / 3600;
    int duration_min = (duration % 3600) / 60;
    int duration_sec = duration % 60;
    netfil << "Started: " << ctime(&start_time);
    netfil << "Ended: " << ctime(&end_time);
    netfil << "Time taken (h:m:s): " << duration_h << ":" <<
        duration_min << ":" << duration_sec << endl;

    netfil.close();
}
