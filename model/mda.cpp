#include <iostream>
#include <fstream>
#include <cstring>
#include "mda.h"
using namespace std;

int count_mda_scenarios(string filename){
    //counting the number of different MDA scenerios that we are testing!

    ifstream in;
    string line;
    
    in.open(filename.c_str());
    if(!in){
        cout << "open " << filename << " failed" << endl;
        exit(1);
    }
    int mda_scenario_num {-1}; //Count the number of lines (not including the first, which is just the column titles)

    while (getline(in, line)){
        ++mda_scenario_num;
    }
    in.close();
    return mda_scenario_num;
}

MDAStrat get_mda_strat(string filename, int N)
{
    ifstream in;
    
    in.open(filename.c_str());
    if(!in){
        cout << "open " << filename << " failed" << endl;
        exit(1);
    }
    
    string line;

    //skip N lines
    for(int i = 0; i < N; ++i){
        getline(in, line);
    }
    
    getline(in,line);

    char *str = new char[line.size()+1];
    strcpy(str, line.c_str());
    char *p = NULL;

    p = strtok(str, ",");       double coverage = atof(p);
    p = strtok(NULL, ",");      double kill_prob = atof(p);
    p = strtok(NULL, ",");      double full_ster_prob = atof(p);
    p = strtok(NULL, ",");      double part_ster_prob = atof(p);
    p = strtok(NULL, ",");      double ster_dur = atof(p);
    p = strtok(NULL, ",");      double part_ster_magnitude = atof(p);
    p = strtok(NULL, ",");      int min_age = atoi(p);
    p = strtok(NULL, ",");      int mda_start_year = atoi(p);
    p = strtok(NULL, ",");      int mda_num_round = atoi(p);
    p = strtok(NULL, ",");      int mda_years_between_rounds = atoi(p);
    p = strtok(NULL, ",");      int num_sims = atoi(p);
   
    delete []str;

    Drugs drug {kill_prob, full_ster_prob, part_ster_prob, ster_dur, part_ster_magnitude};
    MDAStrat strat {coverage, drug, min_age, mda_start_year, mda_num_round, mda_years_between_rounds, num_sims};

    return strat;
}