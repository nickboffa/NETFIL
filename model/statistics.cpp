#include "network.h"
#include "rng.h"
#include <array>
#include <cmath>
#include <cstring>

extern string prv_out_loc;
extern int sim_i;

void Region::output_epidemics(int year, int day, MDAEvent evt){
    
    //total pop
    double pop_total = 0;
    double mf_total = 0;
    double ant_total = 0;

    //worm burdens
    double immature_worm_only = 0;
    double non_mated_adult = 0;
    // Breakdown of *why* an antigen-positive person isn't Mf-positive:
    // mda_sterilised_pair + natural_single_sex == non_mated_adult, and
    // mf_total + non_mated_adult + residual_ant_only == ant_total.
    double mda_sterilised_pair = 0;  // status 'U', but mature worms of BOTH sexes are present — MDA zeroed one side's contribution (Agent::mda_sterile)
    double natural_single_sex  = 0;  // status 'U', only ever had mature worm(s) of one sex — no MDA involvement
    double residual_ant_only   = 0;  // no live mature worms at all (status S/E) — Ag+ only from daily_prob_lose_ant decay since last_mworm_time
    array<double, 10> mated_adult{};

    const array<const char*, 10> MATED_ADULT_LABELS = {
        "one_mated_adult", "two_mated_adult", "three_mated_adult",
        "four_mated_adult", "five_mated_adult", "six_mated_adult",
        "seven_mated_adult", "eight_mated_adult", "nine_mated_adult",
        "tenplus_mated_adult"
    };
    
    vector<double> inf_groups;
    inf_groups.resize(groups.size());

    vector<double> antigen_pos_groups;
    antigen_pos_groups.resize(groups.size());

    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //going through groups
        Group *grp = j->second;
        pop_total += grp->group_pop.size();

        //now over people
        for(map<int,Agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
            Agent *agt = k->second;
            int age = int(agt->age/365);

            if(agt->status == 'I'){//person is infectious
                ++inf_groups[j->first - 1];
                ++mf_total;
                // Count number of mated adult worms
                int bucket = min(static_cast<int>(ceil(agt->worm_strength)) - 1, 9);
                ++mated_adult[bucket];
            }
            //all people infected with any number of mature worms or who still have lingering antibodies are counted
            if(agt->status == 'I' || agt->status == 'U'){
                ++antigen_pos_groups[j->first - 1];
                ++ant_total;
            } else if (random_real() < pow(daily_prob_lose_ant, (year*365 +day) - agt->last_mworm_time)) {
                ++antigen_pos_groups[j->first - 1];
                ++ant_total;
                ++residual_ant_only;
            }
            if (agt->status == 'U') {
                ++non_mated_adult;
                // Both sexes present as live mature worms, but overall still
                // 'U' (not 'I') — the only way that happens is MDA zeroing
                // out one sex's contribution via mda_sterile (see
                // Agent::update() in agent.cpp). Otherwise this agent simply
                // never had an opposite-sex mature worm to begin with.
                bool has_mature_male = false, has_mature_female = false;
                for (Worm* w : agt->wvec) {
                    if (w->status == 'M') {
                        if (w->sex == 'M') has_mature_male = true;
                        else has_mature_female = true;
                    }
                }
                if (has_mature_male && has_mature_female) ++mda_sterilised_pair;
                else ++natural_single_sex;
            }
            if (agt->status == 'E') ++immature_worm_only;
        }

    }
    if (day == 0){
    cout << endl;
    
    cout << "Start of " << year+start_year << ": " << "prepatent = " << pre_indiv.size() << " uninfectious = " << uninf_indiv.size() << " infectious = " << inf_indiv.size() << " antigen positive = " << ant_total << endl;
    cout << "overall mf prevalence = " << fixed << setprecision(2) << inf_indiv.size()/(double)rpop*100 << "%" << endl;
    cout<< "overall ant prevalence = " << fixed << setprecision(2) << ant_total/(double)rpop*100 << "%" << endl;
    cout<< "overall ratio prevalence = " << fixed << setprecision(2) << ant_total/mf_total << endl;
    }
    string prv_dat = OUTDIR;    prv_dat = prv_dat + prv_out_loc; 
    ofstream out;   ifstream in;
    in.open(prv_dat.c_str()); // try opening the target for output
    if(!in){ // if it doesn't exist write a heading
        out.open(prv_dat.c_str());
        out << "name,";
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
        out << "mf_total,";
        out << "ant_total,";
        out << "number_treated,";
        out << "immature_worm_only,";
        out << "non_mated_adult,";
        out << "mda_sterilised_pair,";
        out << "natural_single_sex,";
        out << "residual_ant_only,";
        for (const auto& label : MATED_ADULT_LABELS) {
            out << label << ",";
        }
        for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            out << "pop_" << group_numbers[j -> second -> gid] << ","; 
        }
        for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            if (j != groups.begin()) out << ",";
            out << "mf_" << group_numbers[j->second->gid];
        }
        out << "\n";
        out.close();
    }
    else in.close();

    //write the prevalence for whole populations, by gender, by age group and for each village
    out.open(prv_dat.c_str(), ios::app);
    
    out << evt.drug_name << ",";
    out << sim_i << ",";
    out << year + start_year << ",";
    out << day << ",";
    out << agg_param << ",";
    out << theta1 << ",";
    out << theta2 << ",";
    out << worktonot << ",";
    out << immature_and_ant  << ",";
    out << immature_to_antigen << ",";
    out << evt.coverage << ",";
    out << evt.drug_params.kill_prob << ",";
    out << evt.drug_params.full_ster_prob << ",";
    out << evt.drug_params.part_ster_prob << ",";
    out << evt.drug_params.ster_dur << ",";
    out << evt.drug_params.part_ster_magnitude << ",";
    out << 0 << ",";
    out << 0 << ",";
    out << 0 << ",";
    out << achieved_coverage[year] << ",";
    out << sim_years << ",";
    out << pop_total << ",";
    out << mf_total << ",";
    out << ant_total << ",";
    out << number_treated[year] << ",";
    out << immature_worm_only << ",";
    out << non_mated_adult << ",";
    out << mda_sterilised_pair << ",";
    out << natural_single_sex << ",";
    out << residual_ant_only << ",";
    for (double count : mated_adult) {
        out << count << ",";
    }
    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        double n_village = (j -> second -> group_pop).size();
        if(n_village==0) out << "NA,"; // there's a chance that populations in small villages might drop to zero - this is to avoid crashes in that situation
        else out << n_village << ",";
    }
    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        if (j != groups.begin()) out << ",";
        double n_village = (j->second->group_pop).size();
        if (n_village == 0) out << "NA";
        else out << inf_groups[j->first - 1];
    }
    out << "\n";
    out.close();
}
