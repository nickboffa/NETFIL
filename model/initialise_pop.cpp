#include "network.h"
#include "rng.h"
#include <stdio.h>
#include <string.h>
#include <cstring>

using namespace std;

//constructer of Region
Region::Region(int rid, string rname){
    this->rid = rid; //region id
    this->rname = rname;
    load_config();

    //Trackers for IDs
    rpop = 0;
    next_aid = 1;
    next_gid = 1;
    group_blocks = 0;

    //Distances
    euclid_dst = NULL;
    road_dst = NULL;

    //Clearing counts region level infection status
    pre_indiv.clear();
    uninf_indiv.clear();
    inf_indiv.clear();

    //recreate population when there multiple simulations
    init = pop_reload();
    
    if (!init){ //first simulation and we havent built population before
        read_groups();
       
        bld_groups();
       
        bld_region_population();

        //now saving all our config files (to use on other runs!)

        //Saving group meta data
        string file = CONFIG;   file = file + rname;    file = file + ".init";
        ofstream out(file.c_str());

        out << std::setprecision(2) << std::setiosflags(std::ios::fixed);
        
        out << "TRUE" << endl;
        out << rpop << endl; //region population
        out << next_aid << endl; //next agent id
        out << next_gid << endl; //next group id
        out << group_blocks << endl; //number of groups

        //saving village names/numbers and locations
        for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            out << j->first << "," << group_numbers[j->first] << "," << j->second->lon << "," << j->second->lat << endl;
        }
        out.close();

        //now saving group data, group by group! 
        for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
            Group *grp = j->second; //individual group data

            string grp_str = group_numbers[grp->gid]; //group name

            //file name 
            file = CONFIG_POP;  file = file + grp_str;  file = file + "_pop.init";

            out.open(file.c_str());
            out << "ID,age" << endl;

            //iterating over agents!
            for(map<int, Agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
                Agent *agt = k->second;
                out << agt->aid << "," << agt->age << endl;
            }
            out.close();
        }
    }

    read_parameters();

}

bool Region::pop_reload(){

    //checking if we have already generated the input!
    string file = CONFIG;   file = file + rname;    file = file + ".init";
    ifstream in(file.c_str());
    
    if (!in) {
        cerr << "pop_reload: could not open " << file << ", building population fresh" << endl;
        return false;
    }

    string line;
    getline(in, line);
    if(line != "TRUE") {
        cerr << "pop_reload: first line was \"" << line << "\" not TRUE" << endl;
        return false;
    }
    
    getline(in, line);      rpop = atoi(line.c_str());
    getline(in, line);      next_aid = atoi(line.c_str());
    getline(in, line);      next_gid = atoi(line.c_str());
    getline(in, line);      group_blocks = atoi(line.c_str());
    
    while(getline(in, line)){
        char *str = new char[line.size()+1];
        std::strcpy(str, line.c_str());
        
        char *p = std::strtok(str, ",");        int id = atoi(p);
        p = std::strtok(NULL, ",");             string grp = p;
        p = std::strtok(NULL, ",");             double lon = atof(p);
        p = std::strtok(NULL, ",");             double lat = atof(p);
        
        group_names.insert(pair<string, int>(grp, id));
        group_numbers.insert(pair<int, string>(id, grp));
        groups.insert(pair<int, Group*>(id, new Group(id, this, lat, lon)));
        
        delete []str;
    }
    in.close();
    
    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        string grp_str = group_numbers[j->first];
        
        file = CONFIG_POP;   file = file + grp_str;   file = file + "_pop.init";
        in.open(file.c_str());
        
        string line;
        getline(in, line);
        while(getline(in, line)){
            char *str = new char[line.size()+1];
            std::strcpy(str, line.c_str());
            
            char *p = std::strtok(str, ",");        int id = atoi(p);
            p = std::strtok(NULL, ",");             int age = atoi(p);
            
            Agent *agt = new Agent(id, agg_param, age);

            j->second->add_member(agt);
            
            delete []str;
        }
        in.close();
    }

    return true;
}

void Region::read_groups(){
    //function to read in group info!

    ifstream in;
    string line, file;
   
    file = DATADIR; file = file + GROUP_DATA;
    in.open(file.c_str());

    getline(in, line);          //skip header

    while(getline(in, line)){
        char* str = new char[line.size()+1];
        std::strcpy(str, line.c_str());
        
        char* p = std::strtok(str, ",");      string name = p;
        p = std::strtok(NULL, ",");     int pop = atoi(p);
        p = std::strtok(NULL, ",");     double lat = atof(p);
        p = std::strtok(NULL, ",");     double log = atof(p);
       
        group_names.insert(pair<string, int>(name, next_gid));
        group_numbers.insert(pair<int, string>(next_gid, name));

        if (pop == 0){
            cout << "Warning!!!!!! 0 population" << endl;
        }

        group_pops.insert(pair<int, int>(next_gid, pop));

        double *r = new double[2];
        r[0] = lat;     r[1] = log;
        group_coords.insert(pair<int, double*>(next_gid++, r));
        delete []str;
    }
    in.close();
    
    group_blocks = (int)group_names.size();

    //calculate cpop
    for(map<int, int>::iterator j = group_pops.begin(); j != group_pops.end(); ++j){
        int gid = j->first;
        rpop += group_pops[gid];
    }
}

void Region::bld_groups(){
    for(map<int, double*>::iterator j = group_coords.begin(); j != group_coords.end(); ++j){
        int gid = j->first;
        double lat = j->second[0], log = j->second[1];
        groups.insert(pair<int, Group*>(gid, new Group(gid, this, lat, log)));
    }
}

void Region::bld_region_population(){

    //going village by village
    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        Group *grp = j->second;
        grp->bld_group_pop();
    }
}

void Region::read_parameters(){

    ifstream in;
    string line, file;

    // Demographic/biological parameters are accessed via cfg.field (set in constructor).

    // Load the distance array selected by DISTANCE_TYPE ('r' = road, 'e' = euclidean).
    // Only one array is ever used; the other stays null.
    double*& active_dst     = (DISTANCE_TYPE == 'r') ? road_dst   : euclid_dst;
    const char* bin_file    = (DISTANCE_TYPE == 'r') ? CAR_DISTANCE_BIN  : CROW_DISTANCE_BIN;
    const char* csv_file    = (DISTANCE_TYPE == 'r') ? CAR_DISTANCE      : CROW_DISTANCE;

    if (active_dst == nullptr && group_blocks > 1) {
        int len = group_blocks*(group_blocks-1)/2;
        active_dst = new double[len];   memset(active_dst, 0, sizeof(double)*len);

        string bin_path = string(DATADIR) + bin_file;
        FILE* bf = fopen(bin_path.c_str(), "rb");

        if (bf) {
            cout << "BUILDING ROADS" << endl;
            fread(active_dst, sizeof(double), len, bf);
            fclose(bf);
            cout << "Roads loaded (" << len << " pairs)" << endl;
        } else {
            cout << "BUILDING ROADS (binary not found, reading CSV — run R to generate .bin)" << endl;

            file = DATADIR;  file += csv_file;
            in.open(file.c_str());
            getline(in, line);
            vector<string> grp_vec;
            {
                char *str = new char[line.size()+1];
                std::strcpy(str, line.c_str());
                char *p = std::strtok(str, ",");
                while (p != NULL) { grp_vec.push_back(p); p = std::strtok(NULL, ","); }
                delete []str;

                while (getline(in, line)) {
                    str = new char[line.size()+1];
                    std::strcpy(str, line.c_str());
                    p = std::strtok(str, ",");
                    int src_id = group_names[string(p)];
                    int index = 0;
                    p = std::strtok(NULL, ",");
                    while (p != NULL) {
                        int tag_id = group_names[grp_vec[index++]];
                        if (tag_id > src_id) {
                            int ii = (src_id-1)*(group_blocks*2-src_id)/2 + tag_id-src_id - 1;
                            active_dst[ii] = atof(p);
                        }
                        p = std::strtok(NULL, ",");
                    }
                    delete []str;
                }
            }
            in.close();
        }
    }
    
    file = DATADIR; file = file + TRAN_PARAM;
    in.open(file.c_str());
    getline(in,line);
    getline(in,line);
    char *str = new char[line.size()+1];
    strcpy(str, line.c_str());
    char *p = NULL;

    p = strtok(str, ",");      double theta_1 = atof(p);
    p = strtok(NULL, ",");     double theta_2 = atof(p);
    p = strtok(NULL, ",");     double k = atof(p);
    p = strtok(NULL, ",");     double w2n = atof(p);

    // Optional 5th column: a fitted AntLoss (e.g. WSM), expressed as an
    // antigen half-life in days rather than the raw daily retention
    // probability — converted via daily_retention_from_halflife() (network.h).
    // Absent in older/other countries' tran_params.csv — leaves the value
    // set by load_config() (default_params.json) untouched.
    p = strtok(NULL, ",");
    if (p != NULL) daily_prob_lose_ant = daily_retention_from_halflife(atof(p));

    // Optional 6th column: a fitted p_mda (effective-coverage discount).
    // Absent in countries that don't fit it — leaves p_mda at its default
    // (1.0, i.e. no discount) set in network.h.
    p = strtok(NULL, ",");
    if (p != NULL) p_mda = atof(p);

    // Optional 7th column: a fitted ster_to_kill (DA/IDA kill-vs-sterilise
    // split). Absent in countries that don't fit it — leaves ster_to_kill at
    // its literature-placeholder default set in network.h.
    p = strtok(NULL, ",");
    if (p != NULL) ster_to_kill = atof(p);

    // Optional 8th column: a fitted mda_total_effect (combined DA/IDA
    // kill+sterilise probability). Absent in countries that don't fit it —
    // leaves mda_total_effect at its literature-placeholder default set in
    // network.h.
    p = strtok(NULL, ",");
    if (p != NULL) mda_total_effect = atof(p);

    delete []str;
    in.close();

    theta1 = theta_1;
    theta2 = theta_2;
    theta3 = 1 / (1 - exp(-theta2));
    agg_param = k;
    agg_scale = 1 / k;
    worktonot  = w2n;

    file = DATADIR; file = file + CLUSTERING_PARAMS;
    in.open(file.c_str());
    getline(in, line);  // skip header
    getline(in, line);
    char *sp_str = new char[line.size()+1];
    strcpy(sp_str, line.c_str());
    char *sp = NULL;
    sp = strtok(sp_str, ",");  sigma_group = atof(sp);
    sp = strtok(NULL, ",");    beta_0 = atof(sp);
    delete []sp_str;
    in.close();

    // ant_0, init_prev_min/max, init_ratio_min/max are set from country.json in main()
}

void Region::reset_population(){
   
    //resetting population
    pre_indiv.clear();
    inf_indiv.clear();
    uninf_indiv.clear();
    no_worms_indiv.clear();
    for(map<int, Group*>::iterator j = groups.begin();  j != groups.end(); ++j){ //iterating through groups
        delete j->second;
    }
    groups.clear();

    for(map<int, double*>::iterator j = group_coords.begin();  j != group_coords.end(); ++j){ //iterating through groups
        delete [] j->second;
    }
    group_coords.clear();

    for(int i = 0; i < N_AGE_GROUPS; ++i){
        pvec[i].clear();
    }
    group_names.clear();
    group_numbers.clear();

    group_pops.clear();
    
    // delete[] road_dst;
    // delete[] euclid_dst;
    // euclid_dst = nullptr;
    // road_dst = nullptr;

    rpop = 0;
    next_aid = 1;
    next_gid = 1;
    group_blocks = 0;

    if(!pop_reload()){
        cout << "reload pop err" << endl;
        exit(1);
    }
    
    read_parameters();
    

}

void Region::reset_prev(){
    //now need to clera worms from people
    for(map<int, Agent*>::iterator j = inf_indiv.begin(); j != inf_indiv.end(); ++j){
        Agent *agt =j->second;
        agt->status = 'S';
        agt->worm_strength = 0;
        for(int i = 0; i < agt->wvec.size(); ++i){
            delete agt->wvec[i];
        }
        agt->wvec.clear();
    }

    for(map<int, Agent*>::iterator j = pre_indiv.begin(); j != pre_indiv.end(); ++j){
        Agent *agt =j->second;
        agt->status = 'S';
        for(int i = 0; i < agt->wvec.size(); ++i){
            delete agt->wvec[i];
        }
        agt->wvec.clear();
    }

    for(map<int, Agent*>::iterator j = uninf_indiv.begin(); j != uninf_indiv.end(); ++j){
        Agent *agt =j->second;
        agt->status = 'S';
        for(int i = 0; i < agt->wvec.size(); ++i){
            delete agt->wvec[i];
        }
        agt->wvec.clear();
    }

    //resetting population
    pre_indiv.clear();
    inf_indiv.clear();
    uninf_indiv.clear();
    no_worms_indiv.clear();
}
//constructer of groups
Group::Group(int gid, Region *rgn, double lat, double lon){
    this->gid = gid;
    this->rgn = rgn;
    this->lat = lat;
    this->lon = lon;

    this->sum_mf = 0;
}

Group::~Group(){
    rgn = NULL;

    for(map<int, Agent*>::iterator j = group_pop.begin(); j != group_pop.end(); ++j)
        delete j->second;
    group_pop.clear();
    for (auto node : commuting_dist) {
        delete node;
    }

    // Clear the vector to remove all pointer elements (now dangling pointers)
    commuting_dist.clear();
    commuting_dist.clear();
    day_population.clear();
    commuting_pop.clear();
    day_population.clear();
    commuting_cumsum.clear();

}

//build individual group populations!
void Group::bld_group_pop(){

    int group_population = rgn->group_pops[gid];
    const auto& age_dist = rgn->age_dist;
    
    while(group_population > 0 ){
        double age_p = random_real();
        for(int i = 0; i < N_AGE_GROUPS; ++i){
            int lower_bound = WIDTH_AGE_GROUPS*i; //lower bound of our age bracket
            int upper_bound = WIDTH_AGE_GROUPS*(i+1) - 1; // upper bound of age bracket
            
            double ll = age_dist[i];
            double uu = age_dist[i+1];
            
            if ((age_p >= ll) && (age_p < uu)){
                int id = rgn->next_aid++;
                int age = 365*(lower_bound + (upper_bound - lower_bound)*random_real()); // age 
                Agent *agt = new Agent(id,rgn->agg_param,age); //creating new agent of correct age!
                add_member(agt);
                break;
            }
        }
        group_population--;
    }
}

void Group::add_member(Agent *agt){
    group_pop.insert(pair<int, Agent*>(agt->aid, agt));
}
