#include "network.h"
#include "rng.h"
#include <cstring>

extern string prv_out_loc;
extern int sim_i;

void Region::output_epidemics(int year, int day, MDAStrat strategy){
    
    //total pop
    double pop_total = 0;
    double inf_total = 0;
    double ant_total = 0;

    //worm burdens
    double immature_worm_only = 0;
    double non_mated_adult = 0;
    double one_mated_adult = 0;
    double two_mated_adult = 0;
    double three_mated_adult = 0;
    double four_mated_adult = 0;
    double five_mated_adult = 0;
    double six_mated_adult = 0;
    double seven_mated_adult = 0;
    double eight_mated_adult = 0;
    double nine_mated_adult = 0;
    double tenplus_mated_adult = 0;
    
    vector<double> inf_groups;
    inf_groups.resize(groups.size());

    vector<double> antigen_pos_groups;
    antigen_pos_groups.resize(groups.size());

    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //going through groups
        Group *grp = j->second;
        pop_total += grp->group_pop.size();

        //now over people
        for(map<int,Agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
            Agent *a = k->second;
            int age = int(a->age/365);

            if(a->status == 'I'){//person is infectious
                ++inf_groups[j->first - 1];
                double ws = a->worm_strength;
                ++inf_total;
                if (ws <= 1) ++one_mated_adult;
                if (ws > 1 && ws <= 2) ++two_mated_adult;
                if (ws > 2 && ws <= 3) ++three_mated_adult;
                if (ws > 3 && ws <= 4) ++four_mated_adult;
                if (ws > 4 && ws <= 5) ++five_mated_adult;
                if (ws > 5 && ws <= 6) ++six_mated_adult;
                if (ws > 6 && ws <= 7) ++seven_mated_adult;
                if (ws > 7 && ws <= 8) ++eight_mated_adult;
                if (ws > 8 && ws <= 9) ++nine_mated_adult;
                if (ws > 9 ) ++tenplus_mated_adult;

            }
            if(a->status == 'I' || a->status == 'U'|| random_real() < pow(DAILY_PROB_LOSE_ANT, (year*365 +day) - a->last_mworm_time) ){ //all people infected with any number of mature worms or who still have lingering antibodies are counted
                
                ++antigen_pos_groups[j->first - 1];
                ++ant_total;
            }
            if (a->status == 'U') ++non_mated_adult;
            if (a->status == 'E') ++immature_worm_only;
        }

    }
    if (day == 0){
    cout << endl;
    
    cout << year+START_YEAR << ": " << "prepatent = " << pre_indiv.size() << " uninfectious = " << uninf_indiv.size() << " infectious = " << inf_indiv.size() << " antigen positive = " << ant_total << endl;
    cout << "overall mf prevalence = " << fixed << setprecision(2) << inf_indiv.size()/(double)rpop*100 << "%" << endl;
    cout<< "overall ant prevalence = " << fixed << setprecision(2) << ant_total/(double)rpop*100 << "%" << endl;
    cout<< "overall ratio prevalence = " << fixed << setprecision(2) << ant_total/inf_total << endl;
    }
    string prv_dat = OUTDIR;    prv_dat = prv_dat + prv_out_loc; 
    ofstream out;   ifstream in;
    in.open(prv_dat.c_str()); // try opening the target for output
    if(!in){ // if it doesn't exist write a heading
        out.open(prv_dat.c_str());
        out << "sim_i,";
        out << "year,";
        out << "day,";
        out << "agg_param,";
        out << "theta1,";
        out << "theta2,";
        out << "worktonot,";
        out << "immature_and_ant,";
        out << "immature_to_antigen,";
        out << "coverage,";
        out << "kill_prob,";
        out << "full_ster_prob,";
        out << "part_ster_prob,";
        out << "ster_dur,";
        out << "part_ster_magnitude,";
        out << "mda_start_year,";
        out << "n_mda_rounds,";
        out << "years_between_rounds,";
        out << "achieved_coverage,";
        out << "sim_years,";
        out << "pop_total,";
        out << "inf_total,";
        out << "ant_total,";
        out << "number_treated,";
        out << "immature_worm_only,";
        out << "non_mated_adult,";
        out << "one_mated_adult,";
        out << "two_mated_adult,";
        out << "three_mated_adult,";
        out << "four_mated_adult,";
        out << "five_mated_adult,";
        out << "six_mated_adult,";
        out << "seven_mated_adult,";
        out << "eight_mated_adult,";
        out << "nine_mated_adult,";
        out << "tenplus_mated_adult,";
        for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            out << "pop_" << group_numbers[j -> second -> gid] << ","; 
        }
        for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            out << "mf_" << group_numbers[j -> second -> gid] << ","; 
        }
        out << endl;
        out.close();
    }
    else in.close();

    //write the prevalence for whole populations, by gender, by age group and for each village
    out.open(prv_dat.c_str(), ios::app);
    
    out << sim_i << ",";
    out << year + START_YEAR << ",";
    out << day << ",";
    out << agg_param << ",";
    out << theta1 << ",";
    out << theta2 << ",";
    out << worktonot << ",";
    out << immature_and_ant  << ",";
    out << immature_to_antigen << ",";
    out << strategy.coverage << ",";
    out << strategy.drug.kill_prob << ",";
    out << strategy.drug.full_ster_prob << ",";
    out << strategy.drug.part_ster_prob << ",";
    out << strategy.drug.ster_dur << ",";
    out << strategy.drug.part_ster_magnitude << ",";
    out << strategy.mda_start_year << ",";
    out << strategy.n_mda_rounds << ",";
    out << strategy.years_between_rounds << ",";
    out << achieved_coverage[year] << ",";
    out << SIM_YEARS << ",";
    out << pop_total << ",";
    out << inf_total << ",";
    out << ant_total << ",";
    out << number_treated[year] << ",";
    out << immature_worm_only << ",";
    out << non_mated_adult << ",";
    out << one_mated_adult << ",";
    out << two_mated_adult<< ",";
    out << three_mated_adult << ",";
    out << four_mated_adult << ",";
    out << five_mated_adult<< ",";
    out << six_mated_adult << ",";
    out << seven_mated_adult<< ",";
    out << eight_mated_adult << ",";
    out << nine_mated_adult << ",";
    out << tenplus_mated_adult<< ",";
    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        double n_village = (j -> second -> group_pop).size();
        if(n_village==0) out << "NA,"; // there's a chance that populations in small villages might drop to zero - this is to avoid crashes in that situation
        else out << n_village << ",";
    }
    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        double n_village = (j -> second -> group_pop).size();
        if(n_village==0) out << "NA,"; // there's a chance that populations in small villages might drop to zero - this is to avoid crashes in that situation
        else out <<  inf_groups[j -> first - 1] << ",";
    }
    out << endl;
    out.close();
}
